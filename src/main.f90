!    SSSSSSSSSSSSSSS NNNNNNNN        NNNNNNNN                         CCCCCCCCCCCCC
!  SS:::::::::::::::SN:::::::N       N::::::N                      CCC::::::::::::C
! S:::::SSSSSS::::::SN::::::::N      N::::::N                    CC:::::::::::::::C
! S:::::S     SSSSSSSN:::::::::N     N::::::N                   C:::::CCCCCCCC::::C
! S:::::S            N::::::::::N    N::::::N  aaaaaaaaaaaaa   C:::::C       CCCCCC
! S:::::S            N:::::::::::N   N::::::N  a::::::::::::a C:::::C              
!  S::::SSSS         N:::::::N::::N  N::::::N  aaaaaaaaa:::::aC:::::C              
!   SS::::::SSSSS    N::::::N N::::N N::::::N           a::::aC:::::C              
!     SSS::::::::SS  N::::::N  N::::N:::::::N    aaaaaaa:::::aC:::::C              
!        SSSSSS::::S N::::::N   N:::::::::::N  aa::::::::::::aC:::::C              
!             S:::::SN::::::N    N::::::::::N a::::aaaa::::::aC:::::C              
!             S:::::SN::::::N     N:::::::::Na::::a    a:::::a C:::::C       CCCCCC
! SSSSSSS     S:::::SN::::::N      N::::::::Na::::a    a:::::a  C:::::CCCCCCCC::::C
! S::::::SSSSSS:::::SN::::::N       N:::::::Na:::::aaaa::::::a   CC:::::::::::::::C
! S:::::::::::::::SS N::::::N        N::::::N a::::::::::aa:::a    CCC::::::::::::C
!  SSSSSSSSSSSSSSS   NNNNNNNN         NNNNNNN  aaaaaaaaaa  aaaa       CCCCCCCCCCCCC
!-----------------------------------------------------------------------------------
! **SNaC**, also known as *Slow CaNS*
! Pedro Costa (p.simoes.costa@gmail.com)
!-----------------------------------------------------------------------------------
program snac
  use mpi
  use mod_bound          , only: bounduvw,boundp,updt_rhs
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,myid_block,ierr,comm_block
  use mod_correc         , only: correc
  use mod_initflow       , only: initflow
  use mod_initgrid       , only: initgrid,distribute_grid,bound_grid,save_grid
  use mod_initmpi        , only: initmpi
  use mod_fillps         , only: fillps
  use mod_load           , only: load
  use mod_output         , only: out0d,out1d,write_visu_3d
  use mod_param          , only: read_input, &
                                 datadir,    &
                                 small,      &
                                 rkcoeff,    &
                                 cfl,dtmin,uref,lref,rey,visc,                 &
                                 nstep,time_max,tw_max,stop_type,              &
                                 restart,is_overwrite_save,                    &
                                 nthreadsmax,                                  &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave,     &
                                 dims,lo,hi,lmin,lmax,                         &
                                 gt,gr,                                        &
                                 cbcvel,bcvel,cbcpre,bcpre,                    &
                                 bforce,periods,inivel,                        &
                                 vol_all,my_block,id_first,nblocks,nrank,      &
                                 is_periodic,l_periodic,                       &
                                 hypre_tol,hypre_maxiter
  use mod_updt_pressure  , only: updt_pressure
  use mod_rk             , only: rk_mom
  use mod_sanity         , only: test_sanity
  use mod_solver         , only: init_matrix,create_solver,setup_solver,add_to_diagonal, &
                                 solve_helmholtz,finalize_solver,finalize_matrix,        &
                                 hypre_solver,HYPRESolverPFMG
  use mod_types
  !$ use omp_lib
  implicit none
  integer , dimension(0:1,3) :: nb
  logical , dimension(0:1,3) :: is_bound
  integer , dimension(    3) :: halos
  integer , dimension(    3) :: ng,lo_g,hi_g,lo_1,hi_1
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,up,vp,wp,pp,po
#ifdef _IMPDIFF
  real(rp), allocatable, dimension(:,:,:) :: uo,vo,wo
#endif
  real(rp), allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsp
  real(rp) :: alpha
#ifdef _IMPDIFF
  type(rhs_bound) :: rhsu,rhsv,rhsw
#endif
  real(rp), dimension(0:1,3) :: dl
#ifdef _IMPDIFF
  integer , dimension(    3) :: hiu,hiv,hiw
#endif
  type(hypre_solver) :: psolver
#ifdef _IMPDIFF
  type(hypre_solver) :: usolver,vsolver,wsolver
  real(rp)           :: alphai,alphaoi
#endif
  !
  real(rp) :: dt,dtmax,time,dtrk,divtot,divmax
  integer  :: irk,istep
  real(rp), allocatable, dimension(:) :: dxc  ,dxf  ,xc  ,xf  , &
                                         dyc  ,dyf  ,yc  ,yf  , &
                                         dzc  ,dzf  ,zc  ,zf  , &
                                         dxc_g,dxf_g,xc_g,xf_g, &
                                         dyc_g,dyf_g,yc_g,yf_g, &
                                         dzc_g,dzf_g,zc_g,zf_g
  !
  real(rp), dimension(100) :: var
  character(len=3  ) :: cblock
  character(len=7  ) :: fldnum
  character(len=100) :: filename
  integer :: iunit
  !
  real(rp) :: twi,tw
  logical  :: is_done,kill
#ifdef _TIMING
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  twi = MPI_WTIME()
  !
  ! read parameter file
  !
  call read_input()
  lo_g(:) = lo(:)
  hi_g(:) = hi(:)
  !
  ! check sanity of input file
  !
  call test_sanity(gr,stop_type,cbcvel,cbcpre)
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(my_block,id_first,dims,cbcpre,bcpre,periods,lo,hi,ng,nb,is_bound,halos)
  lo_1(:) = lo(:) - lo_g(:) + 1 ! lo(:) with 1 as first index in the begining of each block
  hi_1(:) = hi(:) - lo_g(:) + 1 ! hi(:) with 1 as first index in the begining of each block
  !
  ! allocate variables
  !
  allocate(u( lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           v( lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           w( lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           p( lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           up(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           vp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           wp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           pp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           po(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
#ifdef _IMPDIFF
  allocate(uo(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           vo(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
           wo(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
#endif
  allocate(dudtrko(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
           dvdtrko(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
           dwdtrko(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
  allocate(dxc(lo(1)-1:hi(1)+1), &
           dxf(lo(1)-1:hi(1)+1), &
            xc(lo(1)-1:hi(1)+1), &
            xf(lo(1)-1:hi(1)+1), &
           dyc(lo(2)-1:hi(2)+1), &
           dyf(lo(2)-1:hi(2)+1), &
            yc(lo(2)-1:hi(2)+1), &
            yf(lo(2)-1:hi(2)+1), &
           dzc(lo(3)-1:hi(3)+1), &
           dzf(lo(3)-1:hi(3)+1), &
            zc(lo(3)-1:hi(3)+1), &
            zf(lo(3)-1:hi(3)+1))
  allocate(dxc_g(lo_g(1)-1:hi_g(1)+1), &
           dxf_g(lo_g(1)-1:hi_g(1)+1), &
            xc_g(lo_g(1)-1:hi_g(1)+1), &
            xf_g(lo_g(1)-1:hi_g(1)+1), &
           dyc_g(lo_g(2)-1:hi_g(2)+1), &
           dyf_g(lo_g(2)-1:hi_g(2)+1), &
            yc_g(lo_g(2)-1:hi_g(2)+1), &
            yf_g(lo_g(2)-1:hi_g(2)+1), &
           dzc_g(lo_g(3)-1:hi_g(3)+1), &
           dzf_g(lo_g(3)-1:hi_g(3)+1), &
            zc_g(lo_g(3)-1:hi_g(3)+1), &
            zf_g(lo_g(3)-1:hi_g(3)+1))
  allocate(rhsp%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsp%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsp%z(lo(1):hi(1),lo(2):hi(2),0:1))
#ifdef _IMPDIFF
  allocate(rhsu%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsu%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsu%z(lo(1):hi(1),lo(2):hi(2),0:1), &
           rhsv%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsv%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsv%z(lo(1):hi(1),lo(2):hi(2),0:1), &
           rhsw%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsw%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsw%z(lo(1):hi(1),lo(2):hi(2),0:1))
#endif
  !
  if(myid == 0) then
    write(stdout,*) '*******************************'
    write(stdout,*) '*** Beginning of simulation ***'
    write(stdout,*) '*******************************'
    write(stdout,*) ''
  endif
  !
  ! generate grid
  !
  call initgrid(lo_g(1),hi_g(1),gt(1),gr(1),lmin(1),lmax(1),dxc_g,dxf_g,xc_g,xf_g)
  call initgrid(lo_g(2),hi_g(2),gt(2),gr(2),lmin(2),lmax(2),dyc_g,dyf_g,yc_g,yf_g)
  call initgrid(lo_g(3),hi_g(3),gt(3),gr(3),lmin(3),lmax(3),dzc_g,dzf_g,zc_g,zf_g)
  write(cblock,'(i3.3)') my_block
  if(myid_block == 0) then
    call save_grid(trim(datadir)//'grid_x_b_'//cblock,lo_g(1),hi_g(1),xf_g,xc_g,dxf_g,dxc_g)
    call save_grid(trim(datadir)//'grid_y_b_'//cblock,lo_g(2),hi_g(2),yf_g,yc_g,dyf_g,dyc_g)
    call save_grid(trim(datadir)//'grid_z_b_'//cblock,lo_g(3),hi_g(3),zf_g,zc_g,dzf_g,dzc_g)
    open(newunit=iunit,status='replace',file=trim(datadir)//'geometry_b_'//cblock//'.out')
      write(iunit,*) lo_g(1),lo_g(2),lo_g(3) 
      write(iunit,*) hi_g(1),hi_g(2),hi_g(3) 
      write(iunit,*) lmin(1),lmin(2),lmin(3) 
      write(iunit,*) lmax(1),lmax(2),lmax(3) 
    close(iunit)
  endif
  call distribute_grid(lo_g(1),lo(1),hi(1),dxc_g,dxc)
  call distribute_grid(lo_g(1),lo(1),hi(1),dxf_g,dxf)
  call distribute_grid(lo_g(1),lo(1),hi(1), xc_g, xc)
  call distribute_grid(lo_g(1),lo(1),hi(1), xf_g, xf)
  call distribute_grid(lo_g(2),lo(2),hi(2),dyc_g,dyc)
  call distribute_grid(lo_g(2),lo(2),hi(2),dyf_g,dyf)
  call distribute_grid(lo_g(2),lo(2),hi(2), yc_g, yc)
  call distribute_grid(lo_g(2),lo(2),hi(2), yf_g, yf)
  call distribute_grid(lo_g(3),lo(3),hi(3),dzc_g,dzc)
  call distribute_grid(lo_g(3),lo(3),hi(3),dzf_g,dzf)
  call distribute_grid(lo_g(3),lo(3),hi(3), zc_g, zc)
  call distribute_grid(lo_g(3),lo(3),hi(3), zf_g, zf)
  call bound_grid(lo_g(1),hi_g(1),lo(1),hi(1),nb(0:1,1),dxc)
  call bound_grid(lo_g(1),hi_g(1),lo(1),hi(1),nb(0:1,1),dxf)
  call bound_grid(lo_g(2),hi_g(2),lo(2),hi(2),nb(0:1,2),dyc)
  call bound_grid(lo_g(2),hi_g(2),lo(2),hi(2),nb(0:1,2),dyf)
  call bound_grid(lo_g(3),hi_g(3),lo(3),hi(3),nb(0:1,3),dzc)
  call bound_grid(lo_g(3),hi_g(3),lo(3),hi(3),nb(0:1,3),dzf)
  !
  ! initialization of the flow fields
  !
  u(:,:,:)      = 0._rp
  v(:,:,:)      = 0._rp
  w(:,:,:)      = 0._rp
  p(:,:,:)      = 0._rp
  if(.not.restart) then
    istep = 0
    time = 0._rp
    call initflow(inivel,.false.,lo,hi,lo_g,hi_g,lmax-lmin,uref,lref,visc,bforce(1), &
                  xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,p)
    if(myid == 0) write(stdout,*) '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)//'fld_b_'//cblock//'.bin',comm_block,ng,[1,1,1],lo_1,hi_1,u,v,w,p,time,istep)
    if(myid == 0) write(stdout,*) '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  endif
  call bounduvw(cbcvel,lo,hi,bcvel,.false.,halos,is_bound,nb, &
                dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
  call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
  up(:,:,:)      = 0._rp
  vp(:,:,:)      = 0._rp
  wp(:,:,:)      = 0._rp
  pp(:,:,:)      = 0._rp
  po(:,:,:)      = 0._rp
  dudtrko(:,:,:) = 0._rp
  dvdtrko(:,:,:) = 0._rp
  dwdtrko(:,:,:) = 0._rp
#ifdef _IMPDIFF
  uo(:,:,:) = 0._rp
  vo(:,:,:) = 0._rp
  wo(:,:,:) = 0._rp
  alphaoi = 0._rp
#endif
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  include 'out1d.h90'
  include 'out3d.h90'
  !
  ! determine time step
  !
  call chkdt(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,visc,u,v,w,dtmax)
  dt = min(cfl*dtmax,dtmin)
  if(myid == 0) write(stdout,*) 'dtmax = ', dtmax, 'dt = ',dt
  !
  ! initialize Poisson solver
  !
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  call init_matrix(cbcpre,bcpre,dl,is_bound,[.true.,.true.,.true.],lo,hi,periods, &
                   dxc,dxf,dyc,dyf,dzc,dzf,rhsp%x,rhsp%y,rhsp%z,psolver)
  call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,psolver)
  call setup_solver(psolver)
#ifdef _IMPDIFF
  dl = reshape([dxf_g(lo_g(1)-0),dxf_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  hiu(:) = hi(:)
  if(is_bound(1,1)) hiu(:) = hiu(:)-[1,0,0]
  call init_matrix(cbcvel(:,:,1),bcvel(:,:,1),dl,is_bound,[.false.,.true.,.true.],lo,hiu,periods, &
                   dxf,dxc,dyc,dyf,dzc,dzf,rhsu%x,rhsu%y,rhsu%z,usolver)
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyf_g(lo_g(2)-0),dyf_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  hiv(:) = hi(:)
  if(is_bound(1,2)) hiv(:) = hiv(:)-[0,1,0]
  call init_matrix(cbcvel(:,:,2),bcvel(:,:,2),dl,is_bound,[.true.,.false.,.true.],lo,hiv,periods, &
                   dxc,dxf,dyf,dyc,dzc,dzf,rhsv%x,rhsv%y,rhsv%z,vsolver)
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzf_g(lo_g(3)-0),dzf_g(hi_g(3))],shape(dl))
  hiw(:) = hi(:)
  if(is_bound(1,3)) hiw(:) = hiw(:)-[0,0,1]
  call init_matrix(cbcvel(:,:,3),bcvel(:,:,3),dl,is_bound,[.true.,.true.,.false.],lo,hiw,periods, &
                   dxc,dxf,dyc,dyf,dzf,dzc,rhsw%x,rhsw%y,rhsw%z,wsolver)
#endif
  !
  ! main loop
  !
  if(myid == 0) write(stdout,*) '*** Calculation loop starts now ***'
  kill    = .false.
  is_done = .false.
  do while(.not.is_done)
#ifdef _TIMING
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time  = time  + dt
    if(myid == 0) write(stdout,*) 'Timestep #', istep, 'Time = ', time
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      alpha = -visc*dtrk/2._rp
      call rk_mom(rkcoeff(:,irk),lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,dt,bforce, &
                  visc,u,v,w,p,dudtrko,dvdtrko,dwdtrko,up,vp,wp)
#ifdef _IMPDIFF
      alphai = alpha**(-1)
      !
      !$OMP WORKSHARE
      up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END WORKSHARE
      call updt_rhs(lo,hiu,is_bound,rhsu%x,rhsu%y,rhsu%z,up)
      call add_to_diagonal(lo,hiu,alphai-alphaoi,usolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,usolver)
      call setup_solver(usolver)
      call solve_helmholtz(usolver,lo,hiu,up,uo)
      call finalize_solver(usolver)
      !
      !$OMP WORKSHARE
      vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END WORKSHARE
      call updt_rhs(lo,hiv,is_bound,rhsv%x,rhsv%y,rhsv%z,vp)
      call add_to_diagonal(lo,hiv,alphai-alphaoi,vsolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,vsolver)
      call setup_solver(vsolver)
      call solve_helmholtz(vsolver,lo,hiv,vp,vo)
      call finalize_solver(vsolver)
      !
      !$OMP WORKSHARE
      wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END WORKSHARE
      call updt_rhs(lo,hiw,is_bound,rhsw%x,rhsw%y,rhsw%z,wp)
      call add_to_diagonal(lo,hiw,alphai-alphaoi,wsolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,wsolver)
      call setup_solver(wsolver)
      call solve_helmholtz(wsolver,lo,hiw,wp,wo)
      call finalize_solver(wsolver)
      !
      alphaoi = alphai
#endif
      call bounduvw(cbcvel,lo,hi,bcvel,.false.,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,up,vp,wp)
#if defined(_IMPDIFF) && defined(_ONE_PRESS_CORR)
      dtrk  = dt
      if(irk < 3) then ! pressure correction only at the last RK step
        !$OMP WORKSHARE
        u(:,:,:) = up(:,:,:)
        v(:,:,:) = vp(:,:,:)
        w(:,:,:) = wp(:,:,:)
        !$OMP END WORKSHARE
        cycle
      endif
#endif
      call fillps(lo,hi,dxf,dyf,dzf,dtrk,up,vp,wp,pp)
      call updt_rhs(lo,hi,is_bound,rhsp%x,rhsp%y,rhsp%z,pp)
      call solve_helmholtz(psolver,lo,hi,pp,po)
      call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,pp)
      call correc(lo,hi,dxc,dyc,dzc,dtrk,pp,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,lo,hi,bcvel,.true.,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
      call updt_pressure(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,pp,p)
      call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
    enddo
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    endif
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    endif
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600._rp
      call MPI_ALLREDUCE(MPI_IN_PLACE,tw,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    endif
    if(mod(istep,icheck) == 0) then
      if(myid == 0) write(stdout,*) 'Checking stability and divergence...'
      !
      call chkdt(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,visc,u,v,w,dtmax)
      dt = min(cfl*dtmax,dtmin)
      if(myid == 0) write(stdout,*) 'dtmax = ', dtmax, 'dt = ',dt
      if(dtmax < small) then
        if(myid == 0) write(stderr,*) 'ERROR: timestep is too small.'
        if(myid == 0) write(stderr,*) 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
      !
      call chkdiv(lo,hi,dxf,dyf,dzf,u,v,w,vol_all,MPI_COMM_WORLD,divtot,divmax)
      if(myid == 0) write(stdout,*) 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
      if(divtot /= divtot) then!divmax > small.or.divtot /= divtot) then
        if(myid == 0) write(stderr,*) 'ERROR: maximum divergence is too large.'
        if(myid == 0) write(stderr,*) 'Aborting...'
        is_done = .true.
        kill = .true.
      endif
    endif
    !
    ! output routines below
    !
    if(mod(istep,iout0d) == 0) then
      !allocate(var(4))
      var(1) = 1._rp*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,var)
    endif
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d) == 0) then
      include 'out1d.h90'
    endif
    if(mod(istep,iout2d) == 0) then
      include 'out2d.h90'
    endif
    if(mod(istep,iout3d) == 0) then
      include 'out3d.h90'
    endif
    if(mod(istep,isave ) == 0.or.(is_done.and..not.kill)) then
      if(is_overwrite_save) then
        filename = 'fld_b_'//cblock//'.bin'
      else
        filename = 'fld_'//fldnum//'_b_'//cblock//'.bin'
      endif
      call load('w',trim(datadir)//trim(filename),comm_block,ng,[1,1,1],lo_1,hi_1,u,v,w,p,time,istep)
      if(.not.is_overwrite_save) then
        !
        ! fld.bin -> last checkpoint file (symbolic link)
        !
        if(myid_block == 0) call execute_command_line('ln -sf '//trim(filename)//' '//trim(datadir)//'fld_b_'//cblock//'.bin')
      endif
      if(myid_block == 0) write(stdout,*) '*** Checkpoint saved at time = ', time, &
                                          'time step = ', istep, 'block = ',my_block,'. ***'
    endif
#ifdef _TIMING
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(myid == 0) write(stdout,*) 'Avrg, min & max elapsed time: '
      if(myid == 0) write(stdout,*) dt12av/(1._rp*nrank),dt12min,dt12max
#endif
  enddo
  call finalize_solver(psolver)
  call finalize_matrix(psolver)
#ifdef _IMPDIFF
  call finalize_matrix(usolver)
  call finalize_matrix(vsolver)
  call finalize_matrix(wsolver)
#endif
  if(myid == 0.and.(.not.kill)) write(stdout,*) '*** Fim ***'
  call MPI_FINALIZE(ierr)
  stop
end program snac
