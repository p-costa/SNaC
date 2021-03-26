module mod_solver
  use mpi_f08
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public init_bc_rhs,init_matrix_3d,create_solver,setup_solver, &
         add_constant_to_diagonal,solve_helmholtz,finalize_matrix,finalize_solver, &
         hypre_solver, &
         HYPRESolverSMG,HYPRESolverPFMG,HYPRESolverGMRES,HYPRESolverBiCGSTAB
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  public init_fft_reduction,init_n_2d_matrices,create_n_solvers,setup_n_solvers, &
         add_constant_to_n_diagonals,solve_n_helmholtz_2d, &
         finalize_n_matrices,finalize_n_solvers
  type alltoallw
    integer            :: counts
    integer            :: disps
    type(MPI_DATATYPE) :: types
  end type alltoallw
#endif
!  ! parameterized derived types like below would make more sense but
!  ! they are very well supported yet
!  type alltoallw(n) ! Parameterized derived types are not working so well yet!
!    integer, len :: n
!    integer           , dimension(n) :: counts
!    integer           , dimension(n) :: disps
!    type(MPI_DATATYPE), dimension(n) :: types
!  end type alltoallw
  integer, parameter :: HYPRESolverSMG      = 1, &
                        HYPRESolverPFMG     = 2, &
                        HYPRESolverGMRES    = 3, &
                        HYPRESolverBiCGSTAB = 4
  type hypre_solver 
    integer(8) :: grid,stencil,precond,solver,mat,rhs,sol
    integer    :: stype,comm_hypre
  end type hypre_solver 
  contains
  subroutine init_bc_rhs(cbc,bc,dl,is_bound,is_centered,lo,hi,periods, &
                         dx1,dx2,dy1,dy2,dz1,dz2,rhsx,rhsy,rhsz)
    !
    ! description
    !
    implicit none
    character(len=1)  , intent(in ), dimension(0:1,3) :: cbc
    real(rp)          , intent(in ), dimension(0:1,3) ::  bc
    real(rp)          , intent(in ), dimension(0:1,3) ::  dl
    logical           , intent(in ), dimension(0:1,3) ::  is_bound
    logical           , intent(in ), dimension(    3) ::  is_centered
    integer           , intent(in ), dimension(    3) :: lo,hi,periods
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dx1,dx2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dy1,dy2
    real(rp)          , intent(in ), target, dimension(lo(3)-1:) :: dz1,dz2
    real(rp)          , intent(out), dimension(lo(2):,lo(3):,0:)    :: rhsx
    real(rp)          , intent(out), dimension(lo(1):,lo(3):,0:)    :: rhsy
    real(rp)          , intent(out), dimension(lo(1):,lo(2):,0:)    :: rhsz
    integer, dimension(3) :: qqq
    real(rp), dimension(0:1,3) :: factor,sgn
    integer :: i,j,k,q,qq,idir,ib
    integer, dimension(3,3) :: eye
    !
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
    eye(:,:) = 0
    do idir=1,3
      eye(idir,idir) = 1
    enddo
    rhsx(:,:,:) = 0._rp
    rhsy(:,:,:) = 0._rp
    rhsz(:,:,:) = 0._rp
    do idir=1,3
      do k=lo(3),hi(3),max((hi(3)-lo(3))*eye(idir,3),1)
        do j=lo(2),hi(2),max((hi(2)-lo(1))*eye(idir,2),1)
          do i=lo(1),hi(1),max((hi(1)-lo(1))*eye(idir,1),1)
            if(periods(idir) == 0) then
              select case(idir)
              case(1)
                if(    i == lo(idir)) then
                  ib = 0
                elseif(i == hi(idir)) then
                  ib = 1
                endif
                  if(is_bound(ib,idir)) &
                  rhsx(j,k,ib) = rhsx(j,k,ib) + factor(ib,idir)/(dx1(i-(1-ib)+qqq(idir))*dx2(i))
              case(2)
                if(    j == lo(idir)) then
                  ib = 0
                elseif(j == hi(idir)) then
                  ib = 1
                endif
                  if(is_bound(ib,idir)) &
                  rhsy(i,k,ib) = rhsy(i,k,ib) + factor(ib,idir)/(dy1(j-(1-ib)+qqq(idir))*dy2(j))
              case(3)
                if(    k == lo(idir)) then
                  ib = 0
                elseif(k == hi(idir)) then
                  ib = 1
                endif
                  if(is_bound(ib,idir)) &
                  rhsz(i,j,ib) = rhsz(i,j,ib) + factor(ib,idir)/(dz1(k-(1-ib)+qqq(idir))*dz2(k))
              end select
            endif
          enddo
        enddo
      enddo
    enddo
  end subroutine init_bc_rhs
  subroutine init_matrix_3d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                            dx1,dx2,dy1,dy2,dz1,dz2,asolver,lambda)
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
    integer           , intent(in ), dimension(    3) :: lo,hi,periods
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dx1,dx2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dy1,dy2
    real(rp)          , intent(in ), target, dimension(lo(3)-1:) :: dz1,dz2
    type(hypre_solver), intent(out)                              :: asolver
    real(rp)          , intent(in ), optional, dimension(:)      :: lambda
    integer, dimension(3         ) :: qqq
    integer, dimension(3,nstencil) :: offsets
    real(rp), dimension(product(hi(:)-lo(:)+1)*nstencil) :: matvalues
    real(rp), dimension(0:1,3) :: factor,sgn
    integer(8) :: grid,stencil,mat,rhs,sol
    integer :: i,j,k,q,qq
    real(rp) :: cc,cxm,cxp,cym,cyp,czm,czp
    integer            :: comm_hypre
    !
    comm_hypre = MPI_COMM_WORLD%MPI_VAL
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
#ifdef _FFT_X
          cxm = 0._rp
          cxp = 0._rp
          qq  = i - lo(1) + 1
#elif  _FFT_Y
          cym = 0._rp
          cyp = 0._rp
          qq  = j - lo(2) + 1
#elif  _FFT_Z
          czm = 0._rp
          czp = 0._rp
          qq  = k - lo(3) + 1
#endif
          cc  = -(cxm+cxp+cym+cyp+czm+czp)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
          cc  = cc + lambda(qq)
#endif
          if(periods(1) == 0) then
            if(is_bound(0,1).and.i == lo(1)) then
              cc = cc + sgn(0,1)*cxm
              cxm = 0._rp
            endif
            if(is_bound(1,1).and.i == hi(1)) then
              cc = cc + sgn(1,1)*cxp
              cxp = 0._rp
            endif
          endif
          if(periods(2) == 0) then
            if(is_bound(0,2).and.j == lo(2)) then
              cc = cc + sgn(0,2)*cym
              cym = 0._rp
            endif
            if(is_bound(1,2).and.j == hi(2)) then
              cc = cc + sgn(1,2)*cyp
              cyp = 0._rp
            endif
          endif
          if(periods(3) == 0) then
            if(is_bound(0,3).and.k == lo(3)) then
              cc = cc + sgn(0,3)*czm
              czm = 0._rp
            endif
            if(is_bound(1,3).and.k == hi(3)) then
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
  end subroutine init_matrix_3d
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
  subroutine add_constant_to_diagonal(lo,hi,alpha,mat)
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
  end subroutine add_constant_to_diagonal
  subroutine solve_helmholtz(asolver,lo,hi,p,po)
    implicit none
    type(hypre_solver), target, intent(in   )               :: asolver
    integer           ,         intent(in   ), dimension(3) :: lo,hi
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p,po
    integer(8), pointer :: solver,mat,rhs,sol
    integer   , pointer :: stype
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    !
    call HYPRE_StructVectorSetBoxValues(rhs,lo,hi,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    call HYPRE_StructVectorAssemble(rhs,ierr)
    !
    ! create soluction vector
    !
    call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
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
    ! fecth results
    !
    call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),ierr)
    po(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
  end subroutine solve_helmholtz
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
  end subroutine finalize_matrix
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  subroutine init_matrix_2d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                            dl1_1,dl1_2,dl2_1,dl2_2,asolver)
    !
    ! description
    !
    implicit none
    integer, parameter :: nstencil = 5
    character(len=1)  , intent(in ), dimension(0:1,2) :: cbc
    real(rp)          , intent(in ), dimension(0:1,2) ::  bc
    real(rp)          , intent(in ), dimension(0:1,2) ::  dl
    logical           , intent(in )                   ::  is_uniform_grid
    logical           , intent(in ), dimension(0:1,2) ::  is_bound
    logical           , intent(in ), dimension(    2) ::  is_centered
    integer           , intent(in ), dimension(    2) :: lo,hi,periods
    real(rp)          , intent(in ), target, dimension(lo(1)-1:) :: dl1_1,dl1_2
    real(rp)          , intent(in ), target, dimension(lo(2)-1:) :: dl2_1,dl2_2
    type(hypre_solver), intent(out)                              :: asolver
    integer, dimension(2         ) :: qqq
    integer, dimension(2,nstencil) :: offsets
    real(rp), dimension(product(hi(:)-lo(:)+1)*nstencil) :: matvalues
    real(rp), dimension(0:1,2) :: factor,sgn
    integer(8) :: grid,stencil,mat,rhs,sol
    integer :: i1,i2,q,qq
    real(rp) :: cc,c1m,c1p,c2m,c2p
    integer            :: comm_hypre
    !
    comm_hypre = MPI_COMM_WORLD%MPI_VAL
    !
    qqq(:) = 0
    where(.not.is_centered(:)) qqq(:) = 1
    factor(:,:) = 0._rp
    sgn(   :,:) = 0._rp
    do q=1,2
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
    ! create 2D grid object
    !
    call HYPRE_StructGridCreate(comm_hypre,2,grid,ierr)
    call HYPRE_StructGridSetPeriodic(grid,periods,ierr)
    call HYPRE_StructGridSetExtents(grid,lo,hi,ierr)
    call HYPRE_StructGridAssemble(grid,ierr)
    !
    ! setup the finite-difference stencil
    !
    call HYPRE_StructStencilCreate(2,nstencil,stencil,ierr)
    offsets = reshape([ 0, 0, &
                       -1, 0, &
                        1, 0, &
                        0,-1, &
                        0, 1 ],shape(offsets))
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
    do i2=lo(2),hi(2)
      do i1=lo(1),hi(1)
        c1m = 1._rp/(dl1_1(i1-1+qqq(1))*dl1_2(i1))
        c1p = 1._rp/(dl1_1(i1  +qqq(1))*dl1_2(i1))
        c2m = 1._rp/(dl2_1(i2-1+qqq(2))*dl2_2(i2))
        c2p = 1._rp/(dl2_1(i2  +qqq(2))*dl2_2(i2))
        cc  = -(c1m+c1p+c2m+c2p)
        if(periods(1) == 0) then
          if(is_bound(0,1).and.i1 == lo(1)) then
            cc = cc + sgn(0,1)*c1m
            c1m = 0._rp
          endif
          if(is_bound(1,1).and.i1 == hi(1)) then
            cc = cc + sgn(1,1)*c1p
            c1p = 0._rp
          endif
        endif
        if(periods(2) == 0) then
          if(is_bound(0,2).and.i2 == lo(2)) then
            cc = cc + sgn(0,2)*c2m
            c2m = 0._rp
          endif
          if(is_bound(1,2).and.i2 == hi(2)) then
            cc = cc + sgn(1,2)*c2p
            c2p = 0._rp
          endif
        endif
        q  = q + 1
        qq = (q-1)*nstencil
        matvalues(qq+1) = cc
        matvalues(qq+2) = c1m
        matvalues(qq+3) = c1p
        matvalues(qq+4) = c2m
        matvalues(qq+5) = c2p
      enddo
    enddo
    call HYPRE_StructMatrixSetBoxValues(mat,lo,hi,nstencil, &
                                        [0,1,2,3,4],matvalues,ierr)
    call HYPRE_StructMatrixAssemble(mat,ierr)
    asolver%grid       = grid
    asolver%stencil    = stencil
    asolver%mat        = mat
    asolver%rhs        = rhs
    asolver%sol        = sol
    asolver%comm_hypre = comm_hypre
  end subroutine init_matrix_2d
  subroutine solve_n_helmholtz_2d(asolver,lo_out,hi_out,lo,hi,p,po)
    implicit none
    type(hypre_solver), target, intent(inout), dimension(:) :: asolver
    integer           ,         intent(in   )               :: lo_out,hi_out
    integer           ,         intent(in   ), dimension(2) :: lo,hi
#ifdef _FFT_Z
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo_out-1:) :: p,po
#elif  _FFT_Y
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo_out-1:,lo(2)-1:) :: p,po
#elif  _FFT_X
    real(rp)          ,         intent(inout), dimension(lo_out-1:,lo(1)-1:,lo(2)-1:) :: p,po
#endif
    integer(8), pointer :: solver,mat,rhs,sol
    integer   , pointer :: stype
    integer :: i_out,ii,i1,i2,n
    do i_out=lo_out,hi_out
      n = i_out-lo_out+1 
      solver  => asolver(n)%solver
      mat     => asolver(n)%mat
      rhs     => asolver(n)%rhs
      sol     => asolver(n)%sol
      stype   => asolver(n)%stype
      !
      ! setup soluction and rhs vectors
      !
#ifdef _FFT_Z
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, p(lo(1):hi(1),lo(2):hi(2),i_out),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(lo(1):hi(1),lo(2):hi(2),i_out),ierr)
#elif  _FFT_Y
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, p(lo(1):hi(1),i_out,lo(2):hi(2)),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(lo(1):hi(1),i_out,lo(2):hi(2)),ierr)
#elif  _FFT_X
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, p(i_out,lo(1):hi(1),lo(2):hi(2)),ierr)
      call HYPRE_StructVectorSetBoxValues(sol,lo,hi,po(i_out,lo(1):hi(1),lo(2):hi(2)),ierr)
#endif
      call HYPRE_StructVectorAssemble(rhs,ierr)
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
      ! fecth results
      !
#ifdef _FFT_Z
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),lo(2):hi(2),i_out),ierr)
      po(lo(1):hi(1),lo(2):hi(2),i_out) = p(lo(1):hi(1),lo(2):hi(2),i_out)
#elif  _FFT_Y
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(lo(1):hi(1),i_out,lo(2):hi(2)),ierr)
      po(lo(1):hi(1),i_out,lo(2):hi(2)) = p(lo(1):hi(1),i_out,lo(2):hi(2))
#elif  _FFT_X
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,p(i_out,lo(1):hi(1),lo(2):hi(2)),ierr)
      po(i_out,lo(1):hi(1),lo(2):hi(2)) = p(i_out,lo(1):hi(1),lo(2):hi(2))
#endif
    enddo
  end subroutine solve_n_helmholtz_2d
  subroutine init_n_2d_matrices(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo_out,hi_out,lo,hi,periods, &
                                dl1_1,dl1_2,dl2_1,dl2_2,lambda,asolver)
    character(len=1)  , intent(in   ), dimension(0:1,2) :: cbc
    real(rp)          , intent(in   ), dimension(0:1,2) ::  bc
    real(rp)          , intent(in   ), dimension(0:1,2) ::  dl
    logical           , intent(in   )                   ::  is_uniform_grid
    logical           , intent(in   ), dimension(0:1,2) ::  is_bound
    logical           , intent(in   ), dimension(    2) ::  is_centered
    integer           , intent(in   )                   :: lo_out,hi_out
    integer           , intent(in   ), dimension(    2) :: lo,hi,periods
    real(rp)          , intent(in   ), target, dimension(lo(1)-1:) :: dl1_1,dl1_2
    real(rp)          , intent(in   ), target, dimension(lo(2)-1:) :: dl2_1,dl2_2
    real(rp)          , intent(in   ), dimension(:) :: lambda
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    type(hypre_solver) :: asolver_aux
    integer :: i_out,q
    do i_out=lo_out,hi_out
      q = i_out-lo_out+1
      asolver_aux = asolver(q)
      call init_matrix_2d(cbc,bc,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                          dl1_1,dl1_2,dl2_1,dl2_2,asolver_aux)
      call add_constant_to_diagonal([lo(1),lo(2),1],[hi(1),hi(2),1],lambda(q),asolver_aux%mat)
      asolver(q) = asolver_aux
    enddo
  end subroutine init_n_2d_matrices
  subroutine create_n_solvers(n,maxiter,maxerror,stype,asolver)
    integer           , intent(   in) :: n
    integer           , intent(   in) :: maxiter
    real(rp)          , intent(   in) :: maxerror
    integer           , intent(   in) :: stype
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call create_solver(maxiter,maxerror,stype,asolver(q))
    enddo
  end subroutine create_n_solvers
  subroutine setup_n_solvers(n,asolver)
    integer           , intent(   in) :: n
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call setup_solver(asolver(q))
    enddo
  end subroutine setup_n_solvers
  subroutine add_constant_to_n_diagonals(n,lo,hi,alpha,mat)
    integer   , intent(in   )               :: n
    integer   , intent(in   ), dimension(2) :: lo,hi
    real(rp)  , intent(in   )               :: alpha
    integer(8), intent(inout), dimension(:) :: mat
    integer :: q
    do q=1,n
      call add_constant_to_diagonal([lo(1),lo(2),1],[hi(1),hi(2),1],alpha,mat(q))
    enddo
  end subroutine add_constant_to_n_diagonals
  subroutine finalize_n_matrices(n,asolver)
    integer           , intent(   in) :: n
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call finalize_matrix(asolver(q))
    enddo
  end subroutine finalize_n_matrices
  subroutine finalize_n_solvers(n,asolver)
    integer           , intent(   in) :: n
    type(hypre_solver), intent(inout), dimension(:) :: asolver
    integer :: q
    do q=1,n
      call finalize_solver(asolver(q))
    enddo
  end subroutine finalize_n_solvers
  subroutine init_fft_reduction(idir,n,bc,is_centered,dl,arrplan,normfft,lambda)
    use iso_c_binding , only: C_PTR
    use mod_fft       , only: fftini,eigenvalues
    implicit none
    integer         , intent(in )                 :: idir
    integer         , intent(in ), dimension(3  ) :: n
    character(len=1), intent(in ), dimension(0:1) :: bc
    logical         , intent(in )                 :: is_centered
    real(rp)        , intent(in )                 :: dl       
    type(C_PTR)     , intent(out), dimension(2  ) :: arrplan
    real(rp)        , intent(out)                 :: normfft
    real(rp)        , intent(out), dimension(:  ) :: lambda
    call fftini(idir,n,bc,is_centered,arrplan,normfft)
    call eigenvalues(n(idir),bc,is_centered,lambda)
    lambda(:) = lambda(:)/dl**2
  end subroutine init_fft_reduction
  subroutine solve_n_helmholtz_2d_old(asolver,maxiter,maxerror,lo_out,hi_out,lo,hi,lambda,alpha_old,p,po)
    implicit none
    type(hypre_solver), target, intent(inout)               :: asolver
    integer           ,         intent(in   )               :: maxiter
    real(rp)          ,         intent(in   )               :: maxerror
    integer           ,         intent(in   )               :: lo_out,hi_out
    integer           ,         intent(in   ), dimension(2) :: lo,hi
    real(rp)          ,         intent(inout)               :: alpha_old
    real(rp)          ,         intent(in   ), dimension(:) :: lambda
#ifdef _FFT_Z
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo_out-1:) :: p,po
#elif  _FFT_Y
    real(rp)          ,         intent(inout), dimension(lo(1)-1:,lo_out-1:,lo(2)-1:) :: p,po
#elif  _FFT_X
    real(rp)          ,         intent(inout), dimension(lo_out-1:,lo(1)-1:,lo(2)-1:) :: p,po
#endif
    integer(8), pointer :: solver,mat,rhs,sol
    integer   , pointer :: stype
    real(rp), dimension(product(hi(:)-lo(:)+1)) :: solvalues,rhsvalues
    integer :: i_out,ii,i1,i2,q
    real(rp) :: alpha
    solver  => asolver%solver
    mat     => asolver%mat
    rhs     => asolver%rhs
    sol     => asolver%sol
    stype   => asolver%stype
    do i_out=lo_out,hi_out
      q = 0
      do i2=lo(2),hi(2)
        do i1=lo(1),hi(1)
          q = q + 1
#ifdef _FFT_Z
          rhsvalues(q) = p( i1,i2,i_out)
          solvalues(q) = po(i1,i2,i_out)
#elif  _FFT_Y
          rhsvalues(q) = p( i1,i_out,i2)
          solvalues(q) = po(i1,i_out,i2)
#elif  _FFT_X
          rhsvalues(q) = p( i_out,i1,i2)
          solvalues(q) = po(i_out,i1,i2)
#endif
        enddo
      enddo
      ii = i_out-lo_out+1
      alpha = lambda(ii) - alpha_old
      call add_constant_to_diagonal([lo(1),lo(2),1],[hi(1),hi(2),1],alpha,mat)
      alpha_old = lambda(ii)
      !
      call HYPRE_StructVectorSetBoxValues(rhs,lo,hi, &
                                          rhsvalues,ierr)
      call HYPRE_StructVectorAssemble(rhs,ierr)
      !
      ! create soluction vector
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
      call create_solver(maxiter,maxerror,stype,asolver)
      call setup_solver(asolver)
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
      ! fecth results
      !
      call HYPRE_StructVectorGetBoxValues(sol,lo,hi,solvalues,ierr)
      q = 0
      do i2=lo(2),hi(2)
        do i1=lo(1),hi(1)
          q = q + 1
#ifdef _FFT_Z
          p( i1,i2,i_out) = solvalues(q)
          po(i1,i2,i_out) = p(i1,i2,i_out)
#elif  _FFT_Y
          p( i1,i_out,i2) = solvalues(q)
          po(i1,i_out,i2) = p(i1,i_out,i2)
#elif  _FFT_X
          p( i_out,i1,i2) = solvalues(q)
          po(i_out,i1,i2) = p(i_out,i1,i2)
#endif
        enddo
      enddo
      call finalize_solver(asolver)
    enddo
  end subroutine solve_n_helmholtz_2d_old
  subroutine init_transpose(lo_idir,idir,dims,n_p,n_s,transpose_params,comm_slab)
    use mod_common_mpi, only: myid
    implicit none
    integer, intent(in)                            :: lo_idir,idir
    integer, intent(in), dimension(3)              :: dims,n_p,n_s
    type(alltoallw), dimension(product(dims),2), intent(out) :: transpose_params
    type(MPI_COMM), intent(out)                    :: comm_slab
    integer, dimension(3) :: n_i
    type(MPI_DATATYPE) :: type_sub_p,type_sub_s
    integer(MPI_ADDRESS_KIND) :: lb, ext_rp
    integer, dimension(3,3) :: eye
    integer :: i,j,k,ii,jj,kk,irank
    !
    eye(:,:) = 0
    do i=1,3
      eye(i,i) = 1
    enddo
    select case(idir)
      case(1)
        n_i(:) = [n_s(1),n_p(2),n_p(3)]
      case(2)
        n_i(:) = [n_p(1),n_s(2),n_p(3)]
      case(3)
        n_i(:) = [n_p(1),n_p(2),n_s(3)]
    end select
    transpose_params(:,:)%counts = 1
    call MPI_TYPE_GET_EXTENT(MPI_REAL_RP,lb,ext_rp)
    irank = 0
    do k=0,dims(3)-1
      do j=0,dims(2)-1
        do i=0,dims(1)-1
          ii = irank*n_s(1)*eye(1,idir) + 1
          jj = irank*n_s(2)*eye(2,idir) + 1
          kk = irank*n_s(3)*eye(3,idir) + 1
          transpose_params(irank+1,1)%disps = extent_f(ii,jj,kk,n_p)*ext_rp
          ii = i*n_p(1)*(1-eye(1,idir)) + 1
          jj = j*n_p(2)*(1-eye(2,idir)) + 1
          kk = k*n_p(3)*(1-eye(3,idir)) + 1
          transpose_params(irank+1,2)%disps = extent_f(ii,jj,kk,n_s)*ext_rp
          irank = irank + 1
        enddo
      enddo
    enddo
    call MPI_TYPE_CREATE_SUBARRAY(3,n_p,n_i,[0,0,0],MPI_ORDER_FORTRAN,MPI_REAL_RP,type_sub_p)
    call MPI_TYPE_CREATE_SUBARRAY(3,n_s,n_i,[0,0,0],MPI_ORDER_FORTRAN,MPI_REAL_RP,type_sub_s)
    call MPI_TYPE_COMMIT(type_sub_p)
    call MPI_TYPE_COMMIT(type_sub_s)
    transpose_params(:,1)%types = type_sub_p
    transpose_params(:,2)%types = type_sub_s
    !
    ! determine communicator pertaining to a slab (share the same value of lower bound in the slab-decomposed direction)
    !
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,lo_idir,myid,comm_slab)
  end subroutine init_transpose
  pure integer function extent_f(i,j,k,n) ! memory position of array with point i,j,k assuming Fortran ordering
     implicit none
     integer, intent(in) :: i,j,k,n(3)
     extent_f = (i-1) + (j-1)*n(1) + (k-1)*n(1)*n(2)
  end function extent_f
  subroutine transpose_slab(idir,nhi,nho,nrank,t_param,comm_block,arr_in,arr_out)
    implicit none
    integer, intent(in)                            :: idir,nhi,nho,nrank
    type(alltoallw), dimension(nrank,2), intent(in) :: t_param
    type(MPI_COMM), intent(in)                    :: comm_block
    real(rp), intent(in ), dimension(1-nhi:,1-nhi:,1-nhi:) :: arr_in
    real(rp), intent(out), dimension(1-nho:,1-nho:,1-nho:) :: arr_out
    integer :: is,ir
    !
    if(idir == 1) then
      is = 1
      ir = 2
    else
      is = 2
      ir = 1
    endif
    call MPI_ALLTOALLW(arr_in( 1,1,1),t_param(:,is)%counts,t_param(:,is)%disps,t_param(:,is)%types, &
                       arr_out(1,1,1),t_param(:,ir)%counts,t_param(:,ir)%disps,t_param(:,ir)%types, &
                       comm_block)
  end subroutine transpose_slab
#endif
end module mod_solver
