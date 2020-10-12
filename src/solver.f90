module mod_solver
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public init_matrix,create_solver,setup_solver,add_to_diagonal, &
         solve_helmholtz,finalize_matrix,finalize_solver, &
         hypre_solver, &
         HYPRESolverSMG,HYPRESolverPFMG,HYPRESolverGMRES,HYPRESolverBiCGSTAB
  integer, parameter :: HYPRESolverSMG      = 1, &
                        HYPRESolverPFMG     = 2, &
                        HYPRESolverGMRES    = 3, &
                        HYPRESolverBiCGSTAB = 4
  type hypre_solver 
    integer(8) :: grid,stencil,precond,solver,mat,rhs,sol
    integer    :: stype,comm_hypre
  end type hypre_solver 
  contains
  subroutine init_matrix(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,ng, &
                         dx1,dx2,dy1,dy2,dz1,dz2,rhsx,rhsy,rhsz,asolver)
    !
    ! description
    !
    implicit none
    integer, parameter :: nstencil = 7
    character(len=1)  , intent(in ), dimension(0:1,3) :: cbc
    real(rp)          , intent(in ), dimension(0:1,3) ::  bc
    real(rp)          , intent(in ), dimension(0:1,3) ::  dl
    logical           , intent(in ), dimension(0:1,3) ::  is_bound
    logical           , intent(in )                   ::  is_uniform_grid
    logical           , intent(in ), dimension(    3) ::  is_centered
    integer           , intent(in ), dimension(    3) :: lo,hi,ng
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dx1,dx2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dy1,dy2
    real(rp)          , intent(in ), target, dimension(lo(3)-1:) :: dz1,dz2
    type(hypre_solver), intent(out)                              :: asolver
    real(rp)          , intent(out), dimension(lo(2):,lo(3):,0:)    :: rhsx
    real(rp)          , intent(out), dimension(lo(1):,lo(3):,0:)    :: rhsy
    real(rp)          , intent(out), dimension(lo(1):,lo(2):,0:)    :: rhsz
    integer, dimension(3         ) :: periods,qqq
    integer, dimension(3,nstencil) :: offsets
    real(rp), dimension(product(hi(:)-lo(:)+1)*nstencil) :: matvalues
    real(rp), dimension(0:1,3) :: factor,sgn
    integer(8) :: grid,stencil,mat,rhs,sol
    integer :: i,j,k,q,qq
    real(rp) :: cc,cxm,cxp,cym,cyp,czm,czp
    integer            :: comm_hypre
    !
    comm_hypre = MPI_COMM_WORLD
    periods(:) = 0
    where(cbc(0,:)//cbc(1,:) == 'PP') periods(:) = ng(:)
    qqq(:) = 0
    where(.not.is_centered(:)) qqq(:) = 1
    factor(:,:) = 0._rp
    sgn(   :,:) = 0._rp
    do q=1,3
      do qq=0,1
        if(is_bound(qq,q)) then
          select case(cbc(qq,q))
          case('N')
            factor(qq,q) = 1._rp*dl(qq,q)*bc(qq,q)
            if(qq == 1) factor(qq,q) = -factor(qq,q)
            sgn(   qq,q) = 1._rp
          case('D')
            if(is_centered(q)) then
              factor(qq,q) = -2._rp*bc(qq,q)
              sgn(   qq,q) = -1._rp
            else
              factor(qq,q) = -1._rp*bc(qq,q)
              sgn(   qq,q) =  0._rp
            endif
          end select
        endif
      enddo
    enddo
    !
    ! create 3D grid object
    !
    call HYPRE_StructGridCreate(comm_hypre,3,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
    call HYPRE_StructGridSetExtents(grid,lo,hi,ierr)
    call HYPRE_StructGridAssemble(grid,ierr)
    !
    ! setup the finite-difference stencil
    !
    call HYPRE_StructStencilCreate(3,nstencil,stencil,ierr)
    offsets = reshape([ 0, 0, 0, &
                       -1, 0, 0, &
                        1, 0, 0, &
                        0,-1, 0, &
                        0, 1, 0, &
                        0, 0,-1, &
                        0, 0, 1 ],shape(offsets))
    do q=1,nstencil
      call HYPRE_StructStencilSetElement(stencil,q-1,offsets(:,q),ierr)
    enddo
    !
    ! create coefficient matrix, and solution & right-hand-side vectors
    !
    call HYPRE_StructMatrixCreate(comm_hypre,grid,stencil,mat,ierr)
    if(is_uniform_grid) call HYPRE_StructMatrixSetSymmetric(mat,1,ierr)
    call HYPRE_StructMatrixInitialize(mat,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,sol,ierr)
    call HYPRE_StructVectorInitialize(sol,ierr)
    call HYPRE_StructVectorCreate(comm_hypre,grid,rhs,ierr)
    call HYPRE_StructVectorInitialize(rhs,ierr)
    q = 0
    rhsx(:,:,:) = 0._rp
    rhsy(:,:,:) = 0._rp
    rhsz(:,:,:) = 0._rp
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          cxm = 1._rp/(dx1(i-1+qqq(1))*dx2(i))
          cxp = 1._rp/(dx1(i  +qqq(1))*dx2(i))
          cym = 1._rp/(dy1(j-1+qqq(2))*dy2(j))
          cyp = 1._rp/(dy1(j  +qqq(2))*dy2(j))
          czm = 1._rp/(dz1(k-1+qqq(3))*dz2(k))
          czp = 1._rp/(dz1(k  +qqq(3))*dz2(k))
          cc  = -(cxm+cxp+cym+cyp+czm+czp)
          if(periods(1) == 0) then
            if(is_bound(0,1).and.i == lo(1)) then
              rhsx(j,k,0) = rhsx(j,k,0) + cxm*factor(0,1)
              cc = cc + sgn(0,1)*cxm
              cxm = 0._rp
            endif
            if(is_bound(1,1).and.i == hi(1)) then
              rhsx(j,k,1) = rhsx(j,k,1) + cxp*factor(1,1)
              cc = cc + sgn(1,1)*cxp
              cxp = 0._rp
            endif
          endif
          if(periods(2) == 0) then
            if(is_bound(0,2).and.j == lo(2)) then
              rhsy(i,k,0) = rhsy(i,k,0) + cym*factor(0,2)
              cc = cc + sgn(0,2)*cym
              cym = 0._rp
            endif
            if(is_bound(1,2).and.j == hi(2)) then
              rhsy(i,k,1) = rhsy(i,k,1) + cyp*factor(1,2)
              cc = cc + sgn(1,2)*cyp
              cyp = 0._rp
            endif
          endif
          if(periods(3) == 0) then
            if(is_bound(0,3).and.k == lo(3)) then
              rhsz(i,j,0) = rhsz(i,j,0) + czm*factor(0,3)
              cc = cc + sgn(0,3)*czm
              czm = 0._rp
            endif
            if(is_bound(1,3).and.k == hi(3)) then
              rhsz(i,j,1) = rhsz(i,j,1) + czp*factor(1,3)
              cc = cc + sgn(1,3)*czp
              czp = 0._rp
            endif
          endif
          qq = (q-1)*nstencil
          matvalues(qq+1) = cc
          matvalues(qq+2) = cxm
          matvalues(qq+3) = cxp
          matvalues(qq+4) = cym
          matvalues(qq+5) = cyp
          matvalues(qq+6) = czm
          matvalues(qq+7) = czp
          !
        enddo
      enddo
    enddo
    call HYPRE_StructMatrixSetBoxValues(mat,lo,hi,nstencil, &
                                        [0,1,2,3,4,5,6],matvalues,ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    asolver%grid       = grid
    asolver%stencil    = stencil
    asolver%mat        = mat
    asolver%rhs        = rhs
    asolver%sol        = sol
    asolver%comm_hypre = comm_hypre
  end subroutine init_matrix
  subroutine create_solver(maxiter,maxerror,stype,asolver)
    implicit none
    integer           ,         intent(   in) :: maxiter
    real(rp)          ,         intent(   in) :: maxerror
    integer           ,         intent(   in) :: stype
    type(hypre_solver), target, intent(inout) :: asolver
    integer(8) :: solver,precond
    integer :: precond_id
    !
    ! setup solver
    !
    ! note: this part was taken from the Paris Simulator code
    !       freely available under a GPL license
    !       http://www.ida.upmc.fr/~zaleski/paris
    !
    if     ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGCreate(asolver%comm_hypre,solver,ierr)
      call HYPRE_StructSMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructSMGSetTol(solver,maxerror,ierr)
      call hypre_structSMGsetLogging(solver,1,ierr)
      call HYPRE_StructSMGSetPrintLevel(solver,1,ierr)
    elseif ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGCreate(asolver%comm_hypre,solver,ierr)
      call HYPRE_StructPFMGSetMaxIter(solver,maxiter,ierr)
      call HYPRE_StructPFMGSetTol(solver,maxerror,ierr)
      call HYPRE_structPFMGsetLogging(solver,1,ierr)
      call HYPRE_StructPFMGSetPrintLevel(solver,1,ierr)
      call HYPRE_StructPFMGSetRelChange(solver,1,ierr)
      ! Relaxiation Method: 2 is the fastest if symm matrix
      ! 0: Jacobi
      ! 1: Weighted Jacobi (default)
      ! 2: Red/Black Gauss-Seidel (symmetric: RB pre- and post-relaxation)
      ! 3: Red/Black Gauss-Seidel (nonsymmetric: RB pre- and post-relaxation)
      call HYPRE_StructPFMGSetRelaxType(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPreRelax(solver,1,ierr)
      call HYPRE_StructPFMGSetNumPostRelax(solver,1,ierr)
    elseif ( stype == HYPRESolverGMRES .or. &
             stype == HYPRESolverBiCGSTAB   ) then
      if     (stype == HYPRESolverGMRES) then
        call HYPRE_StructGMRESCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructGMRESSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructGMRESSetTol(solver,maxerror,ierr)
        !call HYPRE_StructGMRESSetLogging(solver, 1 ,ierr)
      elseif (stype == HYPRESolverBiCGSTAB) then
        call HYPRE_StructBiCGSTABCreate(asolver%comm_hypre,solver,ierr)
        call HYPRE_StructBiCGSTABSetMaxIter(solver,maxiter,ierr)
        call HYPRE_StructBiCGSTABSetTol(solver,maxerror,ierr)
      endif
      ! Use PFMG as preconditioner
      call HYPRE_StructPFMGCreate(asolver%comm_hypre,precond,ierr)
      call HYPRE_StructPFMGSetMaxIter(precond,10,ierr)
      call HYPRE_StructPFMGSetTol(precond,0._rp,ierr)
      call HYPRE_StructPFMGSetZeroGuess(precond,ierr)
      call HYPRE_StructPFMGSetRelChange(precond,1,ierr)
      call HYPRE_StructPFMGSetRelaxType(precond,2,ierr)
      precond_id = 1   ! Set PFMG as preconditioner
      if     (stype == HYPRESolverGMRES) then
        call HYPRE_StructGMRESSetPrecond(solver,precond_id,precond,ierr)
      elseif (stype == HYPRESolverBiCGSTAB) then
        call HYPRE_StructBiCGSTABSetPrecond(solver,precond_id,precond,ierr)
      endif
      asolver%precond = precond
    endif
    asolver%solver  = solver
    asolver%stype   = stype
  end subroutine create_solver
  subroutine setup_solver(asolver)
    implicit none
    type(hypre_solver), target, intent(inout) :: asolver
    integer(8), pointer :: solver,mat,rhs,sol
    integer   , pointer :: stype
    !
    solver => asolver%solver
    stype  => asolver%stype
    mat    => asolver%mat
    rhs    => asolver%rhs
    sol    => asolver%sol
    !
    ! setup solver
    !
    ! note: this part was taken from the Paris Simulator code
    !       freely available under a GPL license
    !       http://www.ida.upmc.fr/~zaleski/paris
    !
    if     ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGSetup(solver,mat,rhs,sol,ierr)
    elseif ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGSetup(solver,mat,rhs,sol,ierr)
    elseif ( stype == HYPRESolverGMRES .or. &
             stype == HYPRESolverBiCGSTAB   ) then
      if     (stype == HYPRESolverGMRES) then
        call HYPRE_StructGMRESSetup(solver,mat,rhs,sol,ierr)
      elseif (stype == HYPRESolverBiCGSTAB) then
        call HYPRE_StructBiCGSTABSetup(solver,mat,rhs,sol,ierr)
      endif
    endif
  end subroutine setup_solver
  !
  subroutine add_to_diagonal(lo,hi,alpha,mat)
    implicit none
    integer   , intent(in   ), dimension(3) :: lo,hi
    real(rp)  , intent(in   ) :: alpha
    integer(8), intent(inout) :: mat
    real(rp), dimension(product(hi(:)-lo(:)+1)) :: matvalues
    integer :: i,j,k,q
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q=q+1
          matvalues(q) = alpha
        enddo
      enddo
    enddo
    call HYPRE_StructMatrixAddToBoxValues(mat,lo,hi,1,[0],matvalues,ierr)
  end subroutine add_to_diagonal
  subroutine solve_helmholtz(asolver,lo,hi,p,po)
    implicit none
    type(hypre_solver), target, intent(in   )               :: asolver
    integer           ,         intent(in   ), dimension(3) :: lo,hi
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p,po
    integer(8), pointer :: solver,mat,rhs,sol
    integer   , pointer :: stype
    real(rp), dimension(product(hi(:)-lo(:)+1)) :: solvalues,rhsvalues
    integer :: i,j,k,q
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          rhsvalues(q) = p( i,j,k)
          solvalues(q) = po(i,j,k)
        enddo
      enddo
    enddo
    !
    !
    ! create rhs and solution vector
    !
    call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, &
                                        rhsvalues,ierr)
    call HYPRE_StructVectorAssemble(rhs,ierr)
    !
    call HYPRE_StructVectorSetBoxValues(sol,lo,hi, &
                                        solvalues,ierr)
    call HYPRE_StructVectorAssemble(sol,ierr)
    !
    ! setup solver, and solve
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if ( stype == HYPRESolverSMG ) then 
      call HYPRE_StructSMGSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructSMGGetNumIterations(solver,num_iterations,ierr)
    elseif ( stype == HYPRESolverPFMG ) then  
      call HYPRE_StructPFMGSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructPFMGGetNumIterations(solver,num_iterations,ierr)
    elseif (stype == HYPRESolverGMRES) then 
      call HYPRE_StructGMRESSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructGMRESGetNumIterations(solver, num_iterations,ierr)
    elseif (stype == HYPRESolverBiCGSTAB) then 
      call HYPRE_StructBiCGSTABSolve(solver,mat,rhs,sol,ierr)
      !call HYPRE_StructBiCGSTABGetNumIterations(solver, num_iterations,ierr)
    endif ! stype
    !
    ! end of part based on the Paris Simulator code
    !
    ! retrieve solution
    !
    call HYPRE_StructVectorGetBoxValues(sol,lo,hi,solvalues,ierr)
    q = 0
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          q = q + 1
          p( i,j,k) = solvalues(q)
          po(i,j,k) = p(i,j,k)
        enddo
      enddo
    enddo
  end subroutine solve_helmholtz
  !
  subroutine finalize_solver(asolver)
    implicit none
    type(hypre_solver), target, intent(in) :: asolver
    integer(8), pointer :: precond,solver
    integer   , pointer :: stype
    !
    precond => asolver%precond
    solver  => asolver%solver
    stype   => asolver%stype
    !
    ! note: this part was based on the the Paris Simulator code
    !       freely available under a GPL license; see:
    !       http://www.ida.upmc.fr/~zaleski/paris/
    !
    if     ( stype == HYPRESolverSMG ) then
      call HYPRE_StructSMGDestroy(solver,ierr)
    elseif ( stype == HYPRESolverPFMG ) then
      call HYPRE_StructPFMGDestroy(solver,ierr)
    elseif ( stype == HYPRESolverGMRES ) then
      call HYPRE_StructGMRESDestroy(solver,ierr)
      call HYPRE_StructPFMGDestroy(precond,ierr)
    elseif ( stype == HYPRESolverBiCGSTAB ) then
      call HYPRE_StructBiCGSTABDestroy(solver,ierr)
      call HYPRE_StructPFMGDestroy(precond,ierr)
    endif
    !
  end subroutine finalize_solver
  subroutine finalize_matrix(asolver)
    implicit none
    type(hypre_solver), target, intent(in) :: asolver
    integer(8), pointer :: grid,stencil,mat,rhs,sol
    !
    grid    => asolver%grid
    stencil => asolver%stencil
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    !
    call HYPRE_StructGridDestroy(grid,ierr)
    call HYPRE_StructStencilDestroy(stencil,ierr)
    call HYPRE_StructMatrixDestroy(mat,ierr)
    call HYPRE_StructVectorDestroy(rhs,ierr)
    call HYPRE_StructVectorDestroy(sol,ierr)
    !
  end subroutine finalize_matrix
end module mod_solver
