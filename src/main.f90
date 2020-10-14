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
  use iso_c_binding      , only: C_PTR
  use mpi_f08
  use mod_bound          , only: bounduvw,boundp,updt_rhs,inflow
  use mod_chkdiv         , only: chkdiv
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,myid_block,comm_block
  use mod_correc         , only: correc
  use mod_initflow       , only: initflow,init_inflow
  use mod_initgrid       , only: initgrid,distribute_grid,bound_grid,save_grid
  use mod_initmpi        , only: initmpi
  use mod_fillps         , only: fillps
  use mod_load           , only: load
  use mod_output         , only: out0d,out1d,write_visu_3d
  use mod_param          , only: read_input, &
                                 datadir,    &
                                 small,      &
                                 rkcoeff,    &
                                 cfl,dtmin,uref,lref,rey,visc,             &
                                 nstep,time_max,tw_max,stop_type,          &
                                 restart,is_overwrite_save,                &
                                 nthreadsmax,                              &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 dims,lo,hi,lmin,lmax,                     &
                                 gt,gr,                                    &
                                 cbcvel,bcvel,cbcpre,bcpre,                &
                                 inflow_type,                              &
                                 bforce,periods,inivel,                    &
                                 vol_all,my_block,id_first,nblocks,nrank,  &
                                 is_periodic,l_periodic,                   &
                                 lmax_max,lmin_min,lo_min,hi_max,          &
                                 hypre_tol,hypre_maxiter
  use mod_updt_pressure  , only: updt_pressure
  use mod_rk             , only: rk_mom
  use mod_sanity         , only: test_sanity
  use mod_solver         , only: init_bc_rhs,init_matrix_3d,create_solver,setup_solver, & 
                                 add_constant_to_diagonal,solve_helmholtz,finalize_solver,finalize_matrix, &
                                 hypre_solver,HYPRESolverPFMG
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  use mod_solver         , only: init_fft_reduction,init_n_2d_matrices,create_n_solvers,setup_n_solvers,solve_n_helmholtz_2d, &
                                 add_constant_to_n_diagonals,finalize_n_solvers,finalize_n_matrices
  use mod_fft            , only: fft,fftend
  use mod_sanity         , only: test_sanity_fft
#endif
  use mod_types
  !$ use omp_lib
  implicit none
  integer , dimension(0:1,3) :: nb
  logical , dimension(0:1,3) :: is_bound,is_bound_inflow
  type(MPI_DATATYPE) , dimension(    3) :: halos
  integer , dimension(    3) :: ng,lo_g,hi_g,lo_1,hi_1
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,up,vp,wp,pp,po
  real(rp), allocatable, dimension(:,:,:) :: velin_x,velin_y,velin_z
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
  logical                             :: is_uniform_grid
  logical, dimension(3)               :: is_centered
  !
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  target                              :: dxc,dxf,dyc,dyf,dzc,dzf
  real(rp), pointer, dimension(:)     :: dl1_1,dl1_2,dl2_1,dl2_2
  type(C_PTR)          , dimension(2) :: arrplan_p
  real(rp)                            :: normfft_p
  real(rp), allocatable, dimension(:) :: lambda_p
  real(rp)                            :: alpha_lambda_p
  type(hypre_solver), allocatable, dimension(:) :: psolver_fft
#ifdef _IMPDIFF
  real(rp), pointer, dimension(:)     :: dlu1_1,dlu1_2,dlu2_1,dlu2_2, &
                                         dlv1_1,dlv1_2,dlv2_1,dlv2_2, &
                                         dlw1_1,dlw1_2,dlw2_1,dlw2_2
  type(C_PTR)          , dimension(2) :: arrplan_u,arrplan_v,arrplan_w
  real(rp)                            :: normfft_u,normfft_v,normfft_w
  real(rp), allocatable, dimension(:) :: lambda_u,lambda_v,lambda_w
  real(rp)                            :: alpha_lambda_u, &
                                         alpha_lambda_v, &
                                         alpha_lambda_w
  type(hypre_solver), allocatable, dimension(:) :: usolver_fft, &
                                                   vsolver_fft, &
                                                   wsolver_fft
#endif
#endif
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
  integer :: idir,ib,il,iu,iskip
  !
  call MPI_INIT()
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid)
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
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
#if   defined(_FFT_X) && !(defined(_FFT_Y) || defined(_FFT_Z))
  idir = 1
#elif defined(_FFT_Y) && !(defined(_FFT_X) || defined(_FFT_Z))
  idir = 2
#elif defined(_FFT_Z) && !(defined(_FFT_X) || defined(_FFT_Y))
  idir = 3
#else
  if(myid == 0) write(stderr,*) 'ERROR: there can be only one FFT direction; check the pre-processor flags.'
  if(myid == 0) write(stderr,*) 'Aborting...'
  call MPI_FINALIZE(ierr)
  error stop
#endif
  call test_sanity_fft(dims(idir),lo(idir),hi(idir),lmin(idir),lmax(idir),gr(idir))
#endif
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(my_block,id_first,dims,cbcpre,bcpre,periods,lmin,lmax,gt,gr,lo,hi,ng,nb,is_bound,halos)
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
  allocate(velin_x( lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,0:1), &
           velin_y( lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1,0:1), &
           velin_z( lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,0:1))
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
  call bound_grid(lo_g(1),hi_g(1),lo(1),hi(1),nb(0:1,1),is_periodic(1),lo_min(1),hi_max(1),dxf,dxc)
  call bound_grid(lo_g(2),hi_g(2),lo(2),hi(2),nb(0:1,2),is_periodic(2),lo_min(2),hi_max(2),dyf,dyc)
  call bound_grid(lo_g(3),hi_g(3),lo(3),hi(3),nb(0:1,3),is_periodic(3),lo_min(3),hi_max(3),dzf,dzc)
  is_uniform_grid = all(dzf(:) == dzf(lo(3))) .and. &
                    all(dyf(:) == dyf(lo(2))) .and. &
                    all(dxf(:) == dxf(lo(1))) .and. &
#ifdef _FFT_X
                    dyf(lo(2)) == dzf(lo(3))
#elif  _FFT_Y
                    dxf(lo(1)) == dzf(lo(3))
#elif  _FFT_Z
                    dyf(lo(2)) == dzf(lo(3))
#else
                    dxf(lo(1)) == dyf(lo(2)) .and. dyf(lo(2)) == dzf(lo(3))
#endif
  call mpi_allreduce(MPI_IN_PLACE,is_uniform_grid,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
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
    call load('r',trim(datadir)//'fld_b_'//cblock//'.bin',comm_block,ng,[1,1,1],lo_1,hi_1,u,v,w,p,po,time,istep)
    if(myid == 0) write(stdout,*) '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  endif
  call bounduvw(cbcvel,lo,hi,bcvel,.false.,halos,is_bound,nb, &
                dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
  call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
  velin_x(:,:,:) = 0._rp
  velin_y(:,:,:) = 0._rp
  velin_z(:,:,:) = 0._rp
  do idir=1,3
    do ib=0,1
      is_bound_inflow(ib,idir) = is_bound(ib,idir).and.inflow_type(ib,idir)>0.and.cbcvel(ib,idir,idir)=='D'
      if(is_bound_inflow(ib,idir)) then
        select case(idir)
        case(1)
          il = 2;iu = 3;iskip = 1
          call init_inflow(periods(il:iu:iskip),lo(il:iu:iskip),hi(il:iu:iskip),lmin(il:iu:iskip),lmax(il:iu:iskip), &
                           yc,zc,bcvel(ib,idir,idir),velin_x(:,:,ib))
        case(2)
          il = 1;iu = 3;iskip = 2
          call init_inflow(periods(il:iu:iskip),lo(il:iu:iskip),hi(il:iu:iskip),lmin(il:iu:iskip),lmax(il:iu:iskip), &
                           xc,zc,bcvel(ib,idir,idir),velin_y(:,:,ib))
        case(3)
          il = 1;iu = 2;iskip = 1
          call init_inflow(periods(il:iu:iskip),lo(il:iu:iskip),hi(il:iu:iskip),lmin(il:iu:iskip),lmax(il:iu:iskip), &
                           xc,yc,bcvel(ib,idir,idir),velin_z(:,:,ib))
        end select
      endif
    enddo
  enddo
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
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
#ifdef _FFT_X
  idir = 1
  il = 2;iu = 3;iskip = 1
  dl1_1  => dyc;dl1_2  => dyf;dl2_1  => dzc;dl2_2  => dzf
#ifdef _IMPDIFF
  dlu1_1 => dyc;dlu1_2 => dyf;dlu2_1 => dzc;dlu2_2 => dzf
  dlv1_1 => dyf;dlv1_2 => dyc;dlv2_1 => dzc;dlv2_2 => dzf
  dlw1_1 => dyc;dlw1_2 => dyf;dlw2_1 => dzf;dlw2_2 => dzc
#endif
#elif  _FFT_Y
  idir = 2
  il = 1;iu = 3;iskip = 2
  dl1_1  => dxc;dl1_2  => dxf;dl2_1  => dzc;dl2_2  => dzf
#ifdef _IMPDIFF
  dlu1_1 => dxf;dlu1_2 => dxc;dlu2_1 => dzc;dlu2_2 => dzf
  dlv1_1 => dxc;dlv1_2 => dxf;dlv2_1 => dzc;dlv2_2 => dzf
  dlw1_1 => dxc;dlw1_2 => dxf;dlw2_1 => dzf;dlw2_2 => dzc
#endif
#elif  _FFT_Z
  idir = 3
  il = 1;iu = 2;iskip = 1
  dl1_1  => dxc;dl1_2  => dxf;dl2_1  => dyc;dl2_2  => dyf
#ifdef _IMPDIFF
  dlu1_1 => dxf;dlu1_2 => dxc;dlu2_1 => dyc;dlu2_2 => dyf
  dlv1_1 => dxc;dlv1_2 => dxf;dlv2_1 => dyf;dlv2_2 => dyc
  dlw1_1 => dxc;dlw1_2 => dxf;dlw2_1 => dyc;dlw2_2 => dyf
#endif
#endif
#endif
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  is_centered(:) = [.true.,.true.,.true.]
  call init_bc_rhs(cbcpre,bcpre,dl,is_bound,is_centered,lo,hi,periods, &
                   dxc,dxf,dyc,dyf,dzc,dzf,rhsp%x,rhsp%y,rhsp%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_p(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcpre(:,idir),.true.,dl(0,idir),arrplan_p,normfft_p,lambda_p)
  alpha_lambda_p = 0._rp
  allocate(psolver_fft(hi(idir)-lo(idir)+1))
  call init_n_2d_matrices(cbcpre(:,il:iu:iskip),bcpre(:,il:iu:iskip),dl(:,il:iu:iskip), &
                          is_uniform_grid,is_bound(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo(idir),hi(idir),lo(il:iu:iskip),hi(il:iu:iskip),periods(il:iu:iskip), &
                          dl1_1,dl1_2,dl2_1,dl2_2,lambda_p,psolver_fft)
  call create_n_solvers(hi(idir)-lo(idir)+1,hypre_maxiter,hypre_tol,HYPRESolverPFMG,psolver_fft)
  call setup_n_solvers(hi(idir)-lo(idir)+1,psolver_fft)
#else
  call init_matrix_3d(cbcpre,bcpre,dl,is_uniform_grid,is_bound,is_centered,lo,hi,periods, &
                      dxc,dxf,dyc,dyf,dzc,dzf,psolver)
  call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,psolver)
  call setup_solver(psolver)
#endif
#ifdef _IMPDIFF
  dl = reshape([dxf_g(lo_g(1)-0),dxf_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  hiu(:) = hi(:)
  if(is_bound(1,1)) hiu(:) = hiu(:)-[1,0,0]
  is_centered(:) = [.false.,.true.,.true.]
  call init_bc_rhs(cbcvel(:,:,1),bcvel(:,:,1),dl,is_bound,is_centered,lo,hiu,periods, &
                   dxf,dxc,dyc,dyf,dzc,dzf,rhsu%x,rhsu%y,rhsu%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_u(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcvel(:,idir,1),is_centered(idir),dl(0,idir),arrplan_u,normfft_u,lambda_u)
  alpha_lambda_u = 0._rp
  allocate(usolver_fft(hi(idir)-lo(idir)+1))
  call init_n_2d_matrices(cbcvel(:,il:iu:iskip,1),bcvel(:,il:iu:iskip,1),dl(:,il:iu:iskip), &
                          is_uniform_grid,is_bound(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo(idir),hiu(idir),lo(il:iu:iskip),hiu(il:iu:iskip),periods(il:iu:iskip), &
                          dlu1_1,dlu1_2,dlu2_1,dlu2_2,lambda_u,usolver_fft)
#else
  call init_matrix_3d(cbcvel(:,:,1),bcvel(:,:,1),dl,is_uniform_grid,is_bound,is_centered,lo,hiu,periods, &
                      dxf,dxc,dyc,dyf,dzc,dzf,usolver)
#endif
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyf_g(lo_g(2)-0),dyf_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  hiv(:) = hi(:)
  if(is_bound(1,2)) hiv(:) = hiv(:)-[0,1,0]
  is_centered(:) = [.true.,.false.,.true.]
  call init_bc_rhs(cbcvel(:,:,2),bcvel(:,:,2),dl,is_bound,is_centered,lo,hiv,periods, &
                   dxc,dxf,dyf,dyc,dzc,dzf,rhsv%x,rhsv%y,rhsv%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_v(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcvel(:,idir,2),is_centered(idir),dl(0,idir),arrplan_v,normfft_v,lambda_v)
  alpha_lambda_v = 0._rp
  allocate(vsolver_fft(hi(idir)-lo(idir)+1))
  call init_n_2d_matrices(cbcvel(:,il:iu:iskip,2),bcvel(:,il:iu:iskip,2),dl(:,il:iu:iskip), &
                          is_uniform_grid,is_bound(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo(idir),hiv(idir),lo(il:iu:iskip),hiv(il:iu:iskip),periods(il:iu:iskip), &
                          dlv1_1,dlv1_2,dlv2_1,dlv2_2,lambda_v,vsolver_fft)
#else
  call init_matrix_3d(cbcvel(:,:,2),bcvel(:,:,2),dl,is_uniform_grid,is_bound,is_centered,lo,hiv,periods, &
                      dxc,dxf,dyf,dyc,dzc,dzf,vsolver)
#endif
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzf_g(lo_g(3)-0),dzf_g(hi_g(3))],shape(dl))
  hiw(:) = hi(:)
  if(is_bound(1,3)) hiw(:) = hiw(:)-[0,0,1]
  is_centered(:) = [.true.,.true.,.false.]
  call init_bc_rhs(cbcvel(:,:,3),bcvel(:,:,3),dl,is_bound,is_centered,lo,hiw,periods, &
                   dxc,dxf,dyc,dyf,dzf,dzc,rhsw%x,rhsw%y,rhsw%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_w(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcvel(:,idir,3),is_centered(idir),dl(0,idir),arrplan_w,normfft_w,lambda_w)
  alpha_lambda_w = 0._rp
  allocate(wsolver_fft(hi(idir)-lo(idir)+1))
  call init_n_2d_matrices(cbcvel(:,il:iu:iskip,3),bcvel(:,il:iu:iskip,3),dl(:,il:iu:iskip), &
                          is_uniform_grid,is_bound(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo(idir),hiw(idir),lo(il:iu:iskip),hiw(il:iu:iskip),periods(il:iu:iskip), &
                          dlw1_1,dlw1_2,dlw2_1,dlw2_2,lambda_w,wsolver_fft)
#else
  call init_matrix_3d(cbcvel(:,:,3),bcvel(:,:,3),dl,is_uniform_grid,is_bound,is_centered,lo,hiw,periods, &
                      dxc,dxf,dyc,dyf,dzf,dzc,wsolver)
#endif
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
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call add_constant_to_n_diagonals(hiu(idir)-lo(idir)+1,lo(il:iu:iskip),hiu(il:iu:iskip), &
                                       alphai-alphaoi,usolver_fft(:)%mat) ! correct diagonal term
      call create_n_solvers(hiu(idir)-lo(idir)+1,hypre_maxiter,hypre_tol,HYPRESolverPFMG,usolver_fft)
      call setup_n_solvers(hiu(idir)-lo(idir)+1,usolver_fft)
      call fft(arrplan_u(1),up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      call solve_n_helmholtz_2d(usolver_fft,lo(idir),hiu(idir),lo(il:iu:iskip),hiu(il:iu:iskip),up,uo)
      call fft(arrplan_u(2),up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_u
      call finalize_n_solvers(hiu(idir)-lo(idir)+1,usolver_fft)
#else
      call add_constant_to_diagonal(lo,hiu,alphai-alphaoi,usolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,usolver)
      call setup_solver(usolver)
      call solve_helmholtz(usolver,lo,hiu,up,uo)
      call finalize_solver(usolver)
#endif
      !
      !$OMP WORKSHARE
      vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END WORKSHARE
      call updt_rhs(lo,hiv,is_bound,rhsv%x,rhsv%y,rhsv%z,vp)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call add_constant_to_n_diagonals(hiv(idir)-lo(idir)+1,lo(il:iu:iskip),hiv(il:iu:iskip), &
                                       alphai-alphaoi,vsolver_fft(:)%mat) ! correct diagonal term
      call create_n_solvers(hiv(idir)-lo(idir)+1,hypre_maxiter,hypre_tol,HYPRESolverPFMG,vsolver_fft)
      call setup_n_solvers(hiv(idir)-lo(idir)+1,vsolver_fft)
      call fft(arrplan_v(1),vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      call solve_n_helmholtz_2d(vsolver_fft,lo(idir),hiv(idir),lo(il:iu:iskip),hiv(il:iu:iskip),vp,vo)
      call fft(arrplan_v(2),vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_v
      call finalize_n_solvers(hiv(idir)-lo(idir)+1,vsolver_fft)
#else
      call add_constant_to_diagonal(lo,hiv,alphai-alphaoi,vsolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,vsolver)
      call setup_solver(vsolver)
      call solve_helmholtz(vsolver,lo,hiv,vp,vo)
      call finalize_solver(vsolver)
#endif
      !
      !$OMP WORKSHARE
      wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END WORKSHARE
      call updt_rhs(lo,hiw,is_bound,rhsw%x,rhsw%y,rhsw%z,wp)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call add_constant_to_n_diagonals(hiw(idir)-lo(idir)+1,lo(il:iu:iskip),hiw(il:iu:iskip), &
                                       alphai-alphaoi,wsolver_fft(:)%mat) ! correct diagonal term
      call create_n_solvers(hiw(idir)-lo(idir)+1,hypre_maxiter,hypre_tol,HYPRESolverPFMG,wsolver_fft)
      call setup_n_solvers(hiw(idir)-lo(idir)+1,wsolver_fft)
      call fft(arrplan_w(1),wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      call solve_n_helmholtz_2d(wsolver_fft,lo(idir),hiw(idir),lo(il:iu:iskip),hiw(il:iu:iskip),wp,wo)
      call fft(arrplan_w(2),wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_w
      call finalize_n_solvers(hiw(idir)-lo(idir)+1,wsolver_fft)
#else
      call add_constant_to_diagonal(lo,hiw,alphai-alphaoi,wsolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,HYPRESolverPFMG,wsolver)
      call setup_solver(wsolver)
      call solve_helmholtz(wsolver,lo,hiw,wp,wo)
      call finalize_solver(wsolver)
#endif
      !
      alphaoi = alphai
#endif
      call bounduvw(cbcvel,lo,hi,bcvel,.false.,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,up,vp,wp)
      call inflow(is_bound_inflow,lo,hi,velin_x,velin_y,velin_z,up,vp,wp)
#if !defined(_IMPDIFF) && defined(_ONE_PRESS_CORR)
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
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call fft(arrplan_p(1),pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      call solve_n_helmholtz_2d(psolver_fft,lo(idir),hi(idir),lo(il:iu:iskip),hi(il:iu:iskip),pp,po)
      call fft(arrplan_p(2),pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_p
#else
      call solve_helmholtz(psolver,lo,hi,pp,po)
#endif
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
      call MPI_ALLREDUCE(MPI_IN_PLACE,tw,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
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
      call load('w',trim(datadir)//trim(filename),comm_block,ng,[1,1,1],lo_1,hi_1,u,v,w,p,po,time,istep)
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
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
      if(myid == 0) write(stdout,*) 'Avrg, min & max elapsed time: '
      if(myid == 0) write(stdout,*) dt12av/(1._rp*nrank),dt12min,dt12max
#endif
  enddo
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  call finalize_n_matrices(hi(idir)-lo(idir)+1,psolver_fft)
  call fftend(arrplan_p)
  call finalize_n_solvers(hi(idir)-lo(idir)+1,psolver_fft)
#else
  call finalize_matrix(psolver)
  call finalize_solver(psolver)
#endif
#ifdef _IMPDIFF
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  call finalize_n_matrices(hiu(idir)-lo(idir)+1,usolver_fft)
  call finalize_n_matrices(hiv(idir)-lo(idir)+1,vsolver_fft)
  call finalize_n_matrices(hiw(idir)-lo(idir)+1,wsolver_fft)
  call fftend(arrplan_u)
  call fftend(arrplan_v)
  call fftend(arrplan_w)
#else
  call finalize_matrix(usolver)
  call finalize_matrix(vsolver)
  call finalize_matrix(wsolver)
#endif
#endif
  if(myid == 0.and.(.not.kill)) write(stdout,*) '*** Fim ***'
  call MPI_FINALIZE()
  stop
end program snac
