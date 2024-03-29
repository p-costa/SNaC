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
! **SNaC**
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
  use mod_output         , only: out0d,out1d,out2d,write_visu_3d
  use mod_param          , only: read_input, &
                                 datadir,    &
                                 small,      &
                                 rkcoeff,    &
                                 cfl,dtmin,uref,lref,rey,visc,             &
                                 nstep,time_max,tw_max,stop_type,          &
                                 restart,is_overwrite_save,nsaves_max,     &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                                 dims,ng,lmin,lmax,                        &
                                 gt,gr,                                    &
                                 cbcvel,bcvel,cbcpre,bcpre,                &
                                 inflow_type,                              &
                                 bforce,periods,inivel,                    &
                                 vol_all,my_block,id_first,nblocks,nrank,  &
                                 is_periodic,l_periodic,                   &
                                 lmax_max,lmin_min,                        &
                                 hypre_tol,hypre_maxiter,hypre_solver_i
  use mod_updt_pressure  , only: updt_pressure
  use mod_rk             , only: rk_mom
  use mod_post           , only: cmpt_wall_forces, updt_wall_forces
  use mod_sanity         , only: test_sanity
  use mod_solver         , only: init_bc_rhs,init_matrix_3d,create_solver,setup_solver, &
                                 add_constant_to_diagonal,solve_helmholtz,finalize_solver,finalize_matrix, &
                                 hypre_solver
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  use mod_solver         , only: init_fft_reduction,init_n_2d_matrices,create_n_solvers,setup_n_solvers,solve_n_helmholtz_2d, &
                                 init_n_3d_matrices,solve_n_helmholtz_3d, &
                                 add_constant_to_n_diagonals,finalize_n_solvers,finalize_n_matrices
#ifdef _FFT_USE_SLABS
  use mod_solver         , only: alltoallw,init_comm_slab,init_transpose_slab_uneven,transpose_slab
#endif
  use mod_fft            , only: fft,fftend
  use mod_sanity         , only: test_sanity_fft
#endif
  use mod_types
  !$ use omp_lib
  implicit none
  integer , dimension(0:1,3) :: nb
  logical , dimension(0:1,3) :: is_bound,is_bound_inflow
  type(MPI_DATATYPE) , dimension(    3) :: halos
  integer , dimension(    3) :: lo,hi,lo_g,hi_g,lo_1,hi_1
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
  type(rhs_bound) :: u_in,v_in,w_in
  real(rp) :: alpha,alpha_bc(0:1,1:3),alpha_bc_o(0:1,1:3)
#ifdef _IMPDIFF
  type(rhs_bound) :: rhsu,rhsv,rhsw,bcu,bcv,bcw
#endif
  real(rp), dimension(    0:1,3) :: dl
#ifdef _IMPDIFF
  real(rp), dimension(0:1,3) :: dlu,dlv,dlw
  integer , dimension(    3) :: hiu,hiv,hiw
#endif
  type(hypre_solver) :: psolver
  logical                             :: is_symm_matrix_p
#ifdef _IMPDIFF
  type(hypre_solver) :: usolver,vsolver,wsolver
  real(rp)           :: alphai,alphaoi
  logical                             :: is_symm_matrix_u, &
                                         is_symm_matrix_v, &
                                         is_symm_matrix_w
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
  integer , dimension(    3) :: lo_a,hi_a
  logical , dimension(0:1,3) :: is_bound_a
  target                              :: dxc,dxf,dyc,dyf,dzc,dzf
  real(rp), pointer, dimension(:)     :: dl1_1,dl1_2,dl2_1,dl2_2
  type(C_PTR)          , dimension(2) :: arrplan_p
  real(rp)                            :: normfft_p
  real(rp), allocatable, dimension(:) :: lambda_p,lambda_p_a
  real(rp)                            :: alpha_lambda_p
  type(hypre_solver), allocatable, dimension(:) :: psolver_fft
  type(MPI_COMM)    , allocatable, dimension(:) :: comms_fft
  integer                             :: npsolvers
#ifdef _IMPDIFF
  integer , dimension(    3) :: hiu_a,hiv_a,hiw_a
  real(rp), pointer, dimension(:)     :: dlu1_1,dlu1_2,dlu2_1,dlu2_2, &
                                         dlv1_1,dlv1_2,dlv2_1,dlv2_2, &
                                         dlw1_1,dlw1_2,dlw2_1,dlw2_2
  type(C_PTR)          , dimension(2) :: arrplan_u,arrplan_v,arrplan_w
  real(rp)                            :: normfft_u,normfft_v,normfft_w
  real(rp), allocatable, dimension(:) :: lambda_u,lambda_v,lambda_w
  real(rp), allocatable, dimension(:) :: lambda_u_a,lambda_v_a,lambda_w_a
  real(rp)                            :: alpha_lambda_u, &
                                         alpha_lambda_v, &
                                         alpha_lambda_w
  type(hypre_solver), allocatable, dimension(:) :: usolver_fft, &
                                                   vsolver_fft, &
                                                   wsolver_fft
#endif
#ifdef _FFT_USE_SLABS
  target :: dxc_g,dxf_g,dyc_g,dyf_g,dzc_g,dzf_g
  type(alltoallw), allocatable, dimension(:,:) :: t_params
  integer              , dimension(3)     :: n_s,n_p,lo_s,hi_s
  real(rp), pointer    , dimension(:)     :: dl1_1_g,dl1_2_g,dl2_1_g,dl2_2_g
  real(rp), allocatable, dimension(:,:,:) :: pp_s
  integer                                 :: nrank_block
#ifdef _IMPDIFF
  real(rp), allocatable, dimension(:,:,:) :: up_s,vp_s,wp_s
#endif
#endif
#ifdef _FFT_USE_SLICED_PENCILS
  integer :: nslices
  integer :: q
#ifndef _FFT_USE_SLABS
  integer             , dimension(3  ) :: n_s,lo_s,hi_s
#endif
  integer, allocatable, dimension(:,:) :: lo_sp,hi_sp
#endif
#endif
  real(rp), dimension(0:1,3) :: tau_x    ,tau_y    ,tau_z    , &
                                tau_x_o  ,tau_y_o  ,tau_z_o  , &
                                tau_x_acc,tau_y_acc,tau_z_acc
  !
  real(rp), dimension(100) :: var
  character(len=3  ) :: cblock
  character(len=7  ) :: fldnum
  character(len=4  ) :: chkptnum
  character(len=100) :: filename
  integer :: iunit
  !
  real(rp) :: twi,tw
  integer  :: savecounter,ichkptnum
  logical  :: is_done,kill
#ifdef _TIMING
  real(rp) :: dt12,dt12av,dt12min,dt12max
#endif
  integer :: idir,ib,il,iu,iskip
  !
  call MPI_INIT()
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid)
  twi = MPI_WTIME()
  savecounter = 0
  !
  ! read parameter file
  !
  call read_input()
  !
  ! check sanity of input file
  !
  call test_sanity(ng,dims,gr,stop_type,cbcvel,cbcpre,periods,inflow_type)
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
  call MPI_FINALIZE()
  error stop
#endif
  call test_sanity_fft(dims(idir),ng(idir),lmin(idir),lmax(idir),gr(idir))
#endif
  !
  ! initialize MPI/OpenMP
  !
  !$ call omp_set_num_threads(omp_get_max_threads())
  call initmpi(my_block,nblocks,id_first,dims,cbcpre,bcpre,is_periodic,periods, &
               lmin,lmax,gt,gr,lo,hi,lo_g,hi_g,ng,nb,is_bound,halos)
  lo_1(:) = lo(:) - lo_g(:) + 1 ! lo(:) with 1 as first index in the beginning of each block
  hi_1(:) = hi(:) - lo_g(:) + 1 ! hi(:) with 1 as first index in the beginning of each block
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
  allocate(u_in%x(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,0:1), &
           u_in%y(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1,0:1), &
           u_in%z(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,0:1), &
           v_in%x(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,0:1), &
           v_in%y(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1,0:1), &
           v_in%z(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,0:1), &
           w_in%x(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1,0:1), &
           w_in%y(lo(1)-1:hi(1)+1,lo(3)-1:hi(3)+1,0:1), &
           w_in%z(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,0:1))
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
  allocate( bcu%x(lo(2):hi(2),lo(3):hi(3),0:1), &
            bcu%y(lo(1):hi(1),lo(3):hi(3),0:1), &
            bcu%z(lo(1):hi(1),lo(2):hi(2),0:1), &
            bcv%x(lo(2):hi(2),lo(3):hi(3),0:1), &
            bcv%y(lo(1):hi(1),lo(3):hi(3),0:1), &
            bcv%z(lo(1):hi(1),lo(2):hi(2),0:1), &
            bcw%x(lo(2):hi(2),lo(3):hi(3),0:1), &
            bcw%y(lo(1):hi(1),lo(3):hi(3),0:1), &
            bcw%z(lo(1):hi(1),lo(2):hi(2),0:1))
#endif
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
#ifdef _FFT_USE_SLICED_PENCILS
nslices = min(16,ng(idir))
if(nslices > ng(idir)) then
  if(myid == 0) write(stderr,*) 'ERROR: number of pencil slices cannot exceed the number of grid points along the FFT direction.'
  call MPI_FINALIZE()
  error stop
end if
#if defined(_IMPDIFF) || defined(_FFT_USE_SLABS)
#if   defined(_IMPDIFF)
  if(myid == 0) write(stderr,*) 'ERROR: implicit diffusion not yet supported with "_FFT_USE_SLICED_PENCILS".'
#elif defined(_FFT_USE_SLABS)
  if(myid == 0) write(stderr,*) 'ERROR: cannot select both "_FFT_USE_SLICED_PENCILS" and "_FFT_USE_SLABS".'
#endif
  if(myid == 0) write(stderr,*) 'Aborting...'
  call MPI_FINALIZE()
  error stop
#endif
#endif
  lo_a(:) = lo(:)
  hi_a(:) = hi(:)
#ifdef _FFT_USE_SLABS
  nrank_block = product(dims(:))
  if(ng(idir) >= nrank_block) then
    n_p = hi(:)-lo(:)+1
#ifdef _FFT_X
    n_s = [ng(1)/nrank_block,ng(2),ng(3)]
#elif  _FFT_Y
    n_s = [ng(1),ng(2)/nrank_block,ng(3)]
#elif  _FFT_Z
    n_s = [ng(1),ng(2),ng(3)/nrank_block]
#endif
  else
    write(stderr,*) 'ERROR: number of tasks in block', my_block, 'exceeds the number of grid points along direction', idir,':'
    write(stderr,*) nrank_block, ' > ', ng(idir)
    write(stderr,*) 'Aborting...'
    call MPI_FINALIZE()
    stop
  end if
  !
  ! distribute slab subdomains as evenly as possible
  !
  irk = mod(ng(idir),nrank_block)
  if(myid_block+1 <= irk) n_s(idir) = n_s(idir) + 1
  lo_s(:) = lo_g(:)
  hi_s(:) = hi_g(:)
  lo_s(idir) = (myid_block*n_s(idir)   + 1) - (lo_g(idir) - 1)
  hi_s(idir) = lo_s(idir) + n_s(idir) - 1
  if(myid_block+1 > irk) then
    lo_s(idir) = lo_s(idir) + irk
    hi_s(idir) = hi_s(idir) + irk
  end if
  !
  lo_a(:) = lo_s(:)
  hi_a(:) = hi_s(:)
  allocate(t_params(nrank_block,2))
  deallocate(po)
  allocate(pp_s(lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0), &
           po(  lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0))
  pp_s(:,:,:) = 0._rp
#ifdef _IMPDIFF
  deallocate(uo,vo,wo)
  allocate(up_s(lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0), &
           vp_s(lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0), &
           wp_s(lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0), &
           uo(  lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0), &
           vo(  lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0), &
           wo(  lo_s(1)-0:hi_s(1)+0,lo_s(2)-0:hi_s(2)+0,lo_s(3)-0:hi_s(3)+0))
  up_s(:,:,:) = 0._rp
  vp_s(:,:,:) = 0._rp
  wp_s(:,:,:) = 0._rp
#endif
#endif
  npsolvers = hi_a(idir)-lo_a(idir)+1
#ifdef _FFT_USE_SLICED_PENCILS
  !
  ! distribute sliced pencils as evenly as possible
  !
  allocate(lo_sp(3,nslices),hi_sp(3,nslices))
  do q=1,nslices
    n_s(:) = ng(:)
    n_s(idir) = ng(idir)/nslices
    irk = mod(ng(idir),nslices)
    if(q <= irk) n_s(idir) = n_s(idir) + 1
    lo_s(:) = lo(:)
    hi_s(:) = hi(:)
    lo_s(idir) = lo(idir) + (q-1)*(n_s(idir))
    hi_s(idir) = lo_s(idir) + n_s(idir) - 1
    if(q > irk) then
      lo_s(idir) = lo_s(idir) + irk
      hi_s(idir) = hi_s(idir) + irk
    end if
    lo_sp(:,q) = lo_s(:)
    hi_sp(:,q) = hi_s(:)
  end do
  npsolvers = nslices
#endif
#endif
  !
  if(myid == 0) then
    write(stdout,*) '*******************************'
    write(stdout,*) '*** Beginning of simulation ***'
    write(stdout,*) '*******************************'
    write(stdout,*) ''
  end if
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
  end if
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
  call bound_grid(lo_g(1),hi_g(1),lo(1),hi(1),nb(0:1,1),is_periodic(1),dxf,dxc)
  call bound_grid(lo_g(2),hi_g(2),lo(2),hi(2),nb(0:1,2),is_periodic(2),dyf,dyc)
  call bound_grid(lo_g(3),hi_g(3),lo(3),hi(3),nb(0:1,3),is_periodic(3),dzf,dzc)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
#ifdef _FFT_USE_SLABS
  dxc_g(lo_g(1)-1) = 0._rp
  dxf_g(lo_g(1)-1) = 0._rp
  if(lo(1) == lo_g(1)) then
    dxc_g(lo_g(1)-1) = dxc(lo(1)-1)
    dxf_g(lo_g(1)-1) = dxf(lo(1)-1)
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,dxc_g(lo_g(1)-1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dxf_g(lo_g(1)-1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  dxc_g(hi_g(1)  ) = 0._rp
  dxf_g(hi_g(1)  ) = 0._rp
  dxc_g(hi_g(1)+1) = 0._rp
  dxf_g(hi_g(1)+1) = 0._rp
  if(hi(1) == hi_g(1)) then
    dxc_g(hi_g(1)  ) = dxc(hi(1)  )
    dxf_g(hi_g(1)  ) = dxf(hi(1)  )
    dxc_g(hi_g(1)+1) = dxc(hi(1)+1)
    dxf_g(hi_g(1)+1) = dxf(hi(1)+1)
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,dxc_g(hi_g(1)  ),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dxf_g(hi_g(1)  ),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dxc_g(hi_g(1)+1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dxf_g(hi_g(1)+1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  dyc_g(lo_g(2)-1) = 0._rp
  dyf_g(lo_g(2)-1) = 0._rp
  if(lo(2) == lo_g(2)) then
    dyc_g(lo_g(2)-1) = dyc(lo(2)-1)
    dyf_g(lo_g(2)-1) = dyf(lo(2)-1)
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,dyc_g(lo_g(2)-1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dyf_g(lo_g(2)-1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  dyc_g(hi_g(2)  ) = 0._rp
  dyf_g(hi_g(2)  ) = 0._rp
  dyc_g(hi_g(2)+1) = 0._rp
  dyf_g(hi_g(2)+1) = 0._rp
  if(hi(2) == hi_g(2)) then
    dyc_g(hi_g(2)  ) = dyc(hi(2)  )
    dyf_g(hi_g(2)  ) = dyf(hi(2)  )
    dyc_g(hi_g(2)+1) = dyc(hi(2)+1)
    dyf_g(hi_g(2)+1) = dyf(hi(2)+1)
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,dyc_g(hi_g(2)  ),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dyf_g(hi_g(2)  ),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dyc_g(hi_g(2)+1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dyf_g(hi_g(2)+1),1,MPI_REAL_RP,MPI_MAX,comm_block)
    dzc_g(lo_g(3)-1) = 0._rp
    dzf_g(lo_g(3)-1) = 0._rp
  if(lo(3) == lo_g(3)) then
    dzc_g(lo_g(3)-1) = dzc(lo(3)-1)
    dzf_g(lo_g(3)-1) = dzf(lo(3)-1)
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,dzc_g(lo_g(3)-1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dzf_g(lo_g(3)-1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  dzc_g(hi_g(3)  ) = 0._rp
  dzf_g(hi_g(3)  ) = 0._rp
  dzc_g(hi_g(3)+1) = 0._rp
  dzf_g(hi_g(3)+1) = 0._rp
  if(hi(3) == hi_g(3)) then
    dzc_g(hi_g(3)  ) = dzc(hi(3)  )
    dzf_g(hi_g(3)  ) = dzf(hi(3)  )
    dzc_g(hi_g(3)+1) = dzc(hi(3)+1)
    dzf_g(hi_g(3)+1) = dzf(hi(3)+1)
  end if
  call MPI_ALLREDUCE(MPI_IN_PLACE,dzc_g(hi_g(3)  ),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dzf_g(hi_g(3)  ),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dzc_g(hi_g(3)+1),1,MPI_REAL_RP,MPI_MAX,comm_block)
  call MPI_ALLREDUCE(MPI_IN_PLACE,dzf_g(hi_g(3)+1),1,MPI_REAL_RP,MPI_MAX,comm_block)
#endif
#endif
  is_uniform_grid = all(dzf(:) == dzf(lo(3))) .and. &
                    all(dyf(:) == dyf(lo(2))) .and. &
                    all(dxf(:) == dxf(lo(1))) .and. &
#ifdef _FFT_X
                    dyf(lo(2)) == dzf(lo(3))
#elif  _FFT_Y
                    dxf(lo(1)) == dzf(lo(3))
#elif  _FFT_Z
                    dxf(lo(1)) == dyf(lo(2))
#else
                    dxf(lo(1)) == dyf(lo(2)) .and. &
                    dxf(lo(1)) == dzf(lo(3)) .and. &
                    dyf(lo(2)) == dzf(lo(3))
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
  end if
  call bounduvw(cbcvel,lo,hi,bcvel,.false.,halos,is_bound,nb, &
                dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
  call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
  u_in%x(:,:,:) = 0._rp
  v_in%x(:,:,:) = 0._rp
  w_in%x(:,:,:) = 0._rp
  u_in%y(:,:,:) = 0._rp
  v_in%y(:,:,:) = 0._rp
  w_in%y(:,:,:) = 0._rp
  u_in%z(:,:,:) = 0._rp
  v_in%z(:,:,:) = 0._rp
  w_in%z(:,:,:) = 0._rp
#ifdef _IMPDIFF
  do ib=0,1
    bcu%x(:,:,ib) = bcvel(ib,1,1)
    bcu%y(:,:,ib) = bcvel(ib,2,1)
    bcu%z(:,:,ib) = bcvel(ib,3,1)
    bcv%x(:,:,ib) = bcvel(ib,1,2)
    bcv%y(:,:,ib) = bcvel(ib,2,2)
    bcv%z(:,:,ib) = bcvel(ib,3,2)
    bcw%x(:,:,ib) = bcvel(ib,1,3)
    bcw%y(:,:,ib) = bcvel(ib,2,3)
    bcw%z(:,:,ib) = bcvel(ib,3,3)
  end do
#endif
  !
  ! inflow BCs
  !
  do idir=1,3
    do ib=0,1
      is_bound_inflow(ib,idir) = is_bound(ib,idir).and.inflow_type(ib,idir)>0.and.cbcvel(ib,idir,idir)=='D'
      if(is_bound_inflow(ib,idir)) then
        select case(idir)
        case(1)
          il = 2;iu = 3;iskip = 1
          call init_inflow(inflow_type(ib,idir),periods(il:iu:iskip), &
                           lo(il:iu:iskip),hi(il:iu:iskip),lmin(il:iu:iskip),lmax(il:iu:iskip), &
                           yc,yf,zc,zf,bcvel(ib,idir,idir),lref,visc,u_in%x(:,:,ib),v_in%x(:,:,ib),w_in%x(:,:,ib))
#ifdef _IMPDIFF
  bcu%x(:,:,ib) = u_in%x(lo(2):hi(2),lo(3):hi(3),ib)
  bcv%x(:,:,ib) = v_in%x(lo(2):hi(2),lo(3):hi(3),ib)
  bcw%x(:,:,ib) = w_in%x(lo(2):hi(2),lo(3):hi(3),ib)
#endif
        case(2)
          il = 1;iu = 3;iskip = 2
          call init_inflow(inflow_type(ib,idir),periods(il:iu:iskip), &
                           lo(il:iu:iskip),hi(il:iu:iskip),lmin(il:iu:iskip),lmax(il:iu:iskip), &
                           xc,xf,zc,zf,bcvel(ib,idir,idir),lref,visc,v_in%y(:,:,ib),u_in%y(:,:,ib),w_in%y(:,:,ib))
#ifdef _IMPDIFF
  bcu%y(:,:,ib) = u_in%y(lo(1):hi(1),lo(3):hi(3),ib)
  bcv%y(:,:,ib) = v_in%y(lo(1):hi(1),lo(3):hi(3),ib)
  bcw%y(:,:,ib) = w_in%y(lo(1):hi(1),lo(3):hi(3),ib)
#endif
        case(3)
          il = 1;iu = 2;iskip = 1
          call init_inflow(inflow_type(ib,idir),periods(il:iu:iskip), &
                           lo(il:iu:iskip),hi(il:iu:iskip),lmin(il:iu:iskip),lmax(il:iu:iskip), &
                           xc,xf,yc,yf,bcvel(ib,idir,idir),lref,visc,w_in%z(:,:,ib),u_in%z(:,:,ib),v_in%z(:,:,ib))
#ifdef _IMPDIFF
  bcu%z(:,:,ib) = u_in%z(lo(1):hi(1),lo(2):hi(2),ib)
  bcv%z(:,:,ib) = v_in%z(lo(1):hi(1),lo(2):hi(2),ib)
  bcw%z(:,:,ib) = w_in%z(lo(1):hi(1),lo(2):hi(2),ib)
#endif
        end select
      end if
    end do
  end do
  call inflow(is_bound_inflow,.false.,lo,hi,u_in%x,v_in%x,w_in%x, &
                                            u_in%y,v_in%y,w_in%y, &
                                            u_in%z,v_in%z,w_in%z,u,v,w)
  alpha_bc(  :,:) = 0._rp
  alpha_bc_o(:,:) = alpha_bc(:,:)
  !
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
  tau_x_o(:,:) = 0._rp
  tau_y_o(:,:) = 0._rp
  tau_z_o(:,:) = 0._rp
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
#ifndef _FFT_USE_SLABS
  dl1_1  => dyc;dl1_2  => dyf;dl2_1  => dzc;dl2_2  => dzf
#else
  dl1_1  => dyc_g;dl1_2  => dyf_g;dl2_1  => dzc_g;dl2_2  => dzf_g
#endif
#ifdef _IMPDIFF
#ifndef _FFT_USE_SLABS
  dlu1_1 => dyc;dlu1_2 => dyf;dlu2_1 => dzc;dlu2_2 => dzf
  dlv1_1 => dyf;dlv1_2 => dyc;dlv2_1 => dzc;dlv2_2 => dzf
  dlw1_1 => dyc;dlw1_2 => dyf;dlw2_1 => dzf;dlw2_2 => dzc
#else
  dlu1_1 => dyc_g;dlu1_2 => dyf_g;dlu2_1 => dzc_g;dlu2_2 => dzf_g
  dlv1_1 => dyf_g;dlv1_2 => dyc_g;dlv2_1 => dzc_g;dlv2_2 => dzf_g
  dlw1_1 => dyc_g;dlw1_2 => dyf_g;dlw2_1 => dzf_g;dlw2_2 => dzc_g
#endif
#endif
#elif  _FFT_Y
  idir = 2
  il = 1;iu = 3;iskip = 2
#ifndef _FFT_USE_SLABS
  dl1_1  => dxc;dl1_2  => dxf;dl2_1  => dzc;dl2_2  => dzf
#else
  dl1_1  => dxc_g;dl1_2  => dxf_g;dl2_1  => dzc_g;dl2_2  => dzf_g
#endif
#ifdef _IMPDIFF
#ifndef _FFT_USE_SLABS
  dlu1_1 => dxf;dlu1_2 => dxc;dlu2_1 => dzc;dlu2_2 => dzf
  dlv1_1 => dxc;dlv1_2 => dxf;dlv2_1 => dzc;dlv2_2 => dzf
  dlw1_1 => dxc;dlw1_2 => dxf;dlw2_1 => dzf;dlw2_2 => dzc
#else
  dlu1_1 => dxf_g;dlu1_2 => dxc_g;dlu2_1 => dzc_g;dlu2_2 => dzf_g
  dlv1_1 => dxc_g;dlv1_2 => dxf_g;dlv2_1 => dzc_g;dlv2_2 => dzf_g
  dlw1_1 => dxc_g;dlw1_2 => dxf_g;dlw2_1 => dzf_g;dlw2_2 => dzc_g
#endif
#endif
#elif  _FFT_Z
  idir = 3
  il = 1;iu = 2;iskip = 1
#ifndef _FFT_USE_SLABS
  dl1_1  => dxc;dl1_2  => dxf;dl2_1  => dyc;dl2_2  => dyf
#else
  dl1_1  => dxc_g;dl1_2  => dxf_g;dl2_1  => dyc_g;dl2_2  => dyf_g
#endif
#ifdef _IMPDIFF
#ifndef _FFT_USE_SLABS
  dlu1_1 => dxf;dlu1_2 => dxc;dlu2_1 => dyc;dlu2_2 => dyf
  dlv1_1 => dxc;dlv1_2 => dxf;dlv2_1 => dyf;dlv2_2 => dyc
  dlw1_1 => dxc;dlw1_2 => dxf;dlw2_1 => dyc;dlw2_2 => dyf
#else
  dlu1_1 => dxf_g;dlu1_2 => dxc_g;dlu2_1 => dyc_g;dlu2_2 => dyf_g
  dlv1_1 => dxc_g;dlv1_2 => dxf_g;dlv2_1 => dyf_g;dlv2_2 => dyc_g
  dlw1_1 => dxc_g;dlw1_2 => dxf_g;dlw2_1 => dyc_g;dlw2_2 => dyf_g
#endif
#endif
#endif
#endif
  dl = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dl))
  is_centered(:) = [.true.,.true.,.true.]
  call init_bc_rhs(cbcpre,bcpre,dl,is_bound,is_centered,lo,hi,periods, &
                   dxc,dxf,dyc,dyf,dzc,dzf,rhsp%x,rhsp%y,rhsp%z)
  alpha         = 0._rp
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_p(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcpre(:,idir),.true.,dl(0,idir),arrplan_p,normfft_p,lambda_p)
  alpha_lambda_p = 0._rp
  allocate(psolver_fft(npsolvers))
  allocate(comms_fft(hi_a(idir)-lo_a(idir)+1))
  allocate(lambda_p_a(hi_a(idir)-lo_a(idir)+1))
  is_bound_a(:,:) = is_bound(:,:)
  is_symm_matrix_p = is_uniform_grid
#ifndef _FFT_USE_SLABS
  comms_fft(:) = MPI_COMM_WORLD
  lambda_p_a(:) = lambda_p
#else
  do ib=1,3 ! determine is_bound pertaining to the slabs
    is_bound_a(:,ib) = .false.
    where(cbcpre(0:1,ib) /= 'F') is_bound_a(0:1,ib) = .true.
  end do
  call init_comm_slab(lo(idir),hi(idir),lo_s(idir),hi_s(idir),myid,comms_fft)
  lambda_p_a(:) = lambda_p(lo_s(idir)-lo(idir)+1:hi_s(idir)-lo(idir)+1)
  call init_transpose_slab_uneven(idir,1,0,dims,lo_1-1,lo_s-1,n_p,n_s,comm_block,t_params)
#endif
#ifndef _FFT_USE_SLICED_PENCILS
  call init_n_2d_matrices(cbcpre(:,il:iu:iskip),bcpre(:,il:iu:iskip),dl(:,il:iu:iskip), &
                          is_symm_matrix_p,is_bound_a(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo_a(idir),hi_a(idir),lo_a(il:iu:iskip),hi_a(il:iu:iskip),periods(il:iu:iskip), &
                          dl1_1,dl1_2,dl2_1,dl2_2,alpha,alpha_bc,lambda_p_a,comms_fft,psolver_fft)
#else
  call init_n_3d_matrices(idir,nslices,cbcpre,bcpre,dl,is_symm_matrix_p,is_bound,is_centered,lo,periods, &
                          lo_sp,hi_sp,dxc,dxf,dyc,dyf,dzc,dzf,alpha,alpha_bc,lambda_p_a,psolver_fft)
#endif
  call create_n_solvers(npsolvers,hypre_maxiter,hypre_tol,hypre_solver_i,psolver_fft)
  call setup_n_solvers(npsolvers,psolver_fft)
#else
  call init_matrix_3d(cbcpre,bcpre,dl,is_symm_matrix_p,is_bound,is_centered,lo,hi,periods, &
                      dxc,dxf,dyc,dyf,dzc,dzf,alpha,alpha_bc,psolver)
  call create_solver(hypre_maxiter,hypre_tol,hypre_solver_i,psolver)
  call setup_solver(psolver)
#endif
#ifdef _IMPDIFF
  dlu = reshape([dxf_g(lo_g(1)-0),dxf_g(hi_g(1)), &
                 dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                 dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dlu))
  hiu(:) = hi(:)
  if(is_bound(1,1)) hiu(:) = hiu(:)-[1,0,0]
  is_centered(:) = [.false.,.true.,.true.]
  is_symm_matrix_u = is_uniform_grid
  call init_bc_rhs(cbcvel(:,:,1),bcvel(:,:,1),dlu,is_bound,is_centered,lo,hiu,periods, &
                   dxf,dxc,dyc,dyf,dzc,dzf,rhsu%x,rhsu%y,rhsu%z,bcu%x,bcu%y,bcu%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_u(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcvel(:,idir,1),is_centered(idir),dlu(0,idir),arrplan_u,normfft_u,lambda_u)
  alpha_lambda_u = 0._rp
  allocate(lambda_u_a(hi_a(idir)-lo_a(idir)+1))
  allocate(usolver_fft(hi_a(idir)-lo_a(idir)+1))
#ifndef _FFT_USE_SLABS
  lambda_u_a(:) = lambda_u(:)
  hiu_a(:) = hiu(:)
#else
  lambda_u_a(:) = lambda_u(lo_s(idir)-lo(idir)+1:hi_s(idir)-lo(idir)+1)
  hiu_a(:) = hi_s(:)
  if(periods(1) == 0) hiu_a(:) = hiu_a(:)-[1,0,0]
#endif
  call init_n_2d_matrices(cbcvel(:,il:iu:iskip,1),bcvel(:,il:iu:iskip,1),dlu(:,il:iu:iskip), &
                          is_symm_matrix_u,is_bound_a(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo_a(idir),hiu_a(idir),lo_a(il:iu:iskip),hiu_a(il:iu:iskip),periods(il:iu:iskip), &
                          dlu1_1,dlu1_2,dlu2_1,dlu2_2,alpha,alpha_bc,lambda_u_a,comms_fft,usolver_fft)
#else
  call init_matrix_3d(cbcvel(:,:,1),bcvel(:,:,1),dlu,is_symm_matrix_u,is_bound,is_centered,lo,hiu,periods, &
                      dxf,dxc,dyc,dyf,dzc,dzf,alpha,alpha_bc,usolver)
#endif
  dlv = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                 dyf_g(lo_g(2)-0),dyf_g(hi_g(2)), &
                 dzc_g(lo_g(3)-1),dzc_g(hi_g(3))],shape(dlv))
  hiv(:) = hi(:)
  if(is_bound(1,2)) hiv(:) = hiv(:)-[0,1,0]
  is_centered(:) = [.true.,.false.,.true.]
  call init_bc_rhs(cbcvel(:,:,2),bcvel(:,:,2),dlv,is_bound,is_centered,lo,hiv,periods, &
                   dxc,dxf,dyf,dyc,dzc,dzf,rhsv%x,rhsv%y,rhsv%z,bcv%x,bcv%y,bcv%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_v(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcvel(:,idir,2),is_centered(idir),dlv(0,idir),arrplan_v,normfft_v,lambda_v)
  alpha_lambda_v = 0._rp
  allocate(lambda_v_a(hi_a(idir)-lo_a(idir)+1))
  allocate(vsolver_fft(hi_a(idir)-lo_a(idir)+1))
#ifndef _FFT_USE_SLABS
  lambda_v_a(:) = lambda_v(:)
  hiv_a(:) = hiv(:)
#else
  lambda_v_a(:) = lambda_v(lo_s(idir)-lo(idir)+1:hi_s(idir)-lo(idir)+1)
  hiv_a(:) = hi_s(:)
  if(periods(2) == 0) hiv_a(:) = hiv_a(:)-[0,1,0]
  is_symm_matrix_v = is_uniform_grid
#endif
  call init_n_2d_matrices(cbcvel(:,il:iu:iskip,2),bcvel(:,il:iu:iskip,2),dlv(:,il:iu:iskip), &
                          is_symm_matrix_v,is_bound_a(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo_a(idir),hiv_a(idir),lo_a(il:iu:iskip),hiv_a(il:iu:iskip),periods(il:iu:iskip), &
                          dlv1_1,dlv1_2,dlv2_1,dlv2_2,alpha,alpha_bc,lambda_v_a,comms_fft,vsolver_fft)
#else
  call init_matrix_3d(cbcvel(:,:,2),bcvel(:,:,2),dlv,is_symm_matrix_v,is_bound,is_centered,lo,hiv,periods, &
                      dxc,dxf,dyf,dyc,dzc,dzf,alpha,alpha_bc,vsolver)
#endif
  dlw = reshape([dxc_g(lo_g(1)-1),dxc_g(hi_g(1)), &
                 dyc_g(lo_g(2)-1),dyc_g(hi_g(2)), &
                 dzf_g(lo_g(3)-0),dzf_g(hi_g(3))],shape(dlw))
  hiw(:) = hi(:)
  if(is_bound(1,3)) hiw(:) = hiw(:)-[0,0,1]
  is_centered(:) = [.true.,.true.,.false.]
  call init_bc_rhs(cbcvel(:,:,3),bcvel(:,:,3),dlw,is_bound,is_centered,lo,hiw,periods, &
                   dxc,dxf,dyc,dyf,dzf,dzc,rhsw%x,rhsw%y,rhsw%z,bcw%x,bcw%y,bcw%z)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  allocate(lambda_w(hi(idir)-lo(idir)+1))
  call init_fft_reduction(idir,hi(:)-lo(:)+1,cbcvel(:,idir,3),is_centered(idir),dlw(0,idir),arrplan_w,normfft_w,lambda_w)
  alpha_lambda_w = 0._rp
  allocate(lambda_w_a(hi_a(idir)-lo_a(idir)+1))
  allocate(wsolver_fft(hi_a(idir)-lo_a(idir)+1))
#ifndef _FFT_USE_SLABS
  lambda_w_a(:) = lambda_w(:)
  hiw_a(:) = hiw(:)
#else
  lambda_w_a(:) = lambda_w(lo_s(idir)-lo(idir)+1:hi_s(idir)-lo(idir)+1)
  hiw_a(:) = hi_s(:)
  if(periods(3) == 0) hiw_a(:) = hiw_a(:)-[0,0,1]
  is_symm_matrix_w = is_uniform_grid
#endif
  call init_n_2d_matrices(cbcvel(:,il:iu:iskip,3),bcvel(:,il:iu:iskip,3),dlw(:,il:iu:iskip), &
                          is_symm_matrix_w,is_bound_a(:,il:iu:iskip),is_centered(il:iu:iskip), &
                          lo_a(idir),hiw_a(idir),lo_a(il:iu:iskip),hiw_a(il:iu:iskip),periods(il:iu:iskip), &
                          dlw1_1,dlw1_2,dlw2_1,dlw2_2,alpha,alpha_bc,lambda_w_a,comms_fft,wsolver_fft)
#else
  call init_matrix_3d(cbcvel(:,:,3),bcvel(:,:,3),dlw,is_symm_matrix_w,is_bound,is_centered,lo,hiw,periods, &
                      dxc,dxf,dyc,dyf,dzf,dzc,alpha,alpha_bc,wsolver)
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
    tau_x_acc(:,:) = 0._rp; tau_y_acc(:,:) = 0._rp; tau_z_acc(:,:) = 0._rp
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      alpha = -visc*dtrk/2._rp
      call rk_mom(rkcoeff(:,irk),lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,dt,bforce, &
                  visc,u,v,w,p,dudtrko,dvdtrko,dwdtrko,up,vp,wp)
      call cmpt_wall_forces(hi(:)-lo(:)+1,is_bound,dxc,dxf,dyc,dyf,dzc,dzf,visc, &
                            u,v,w,p,bforce,tau_x,tau_y,tau_z)
      call updt_wall_forces(rkcoeff(:,irk),tau_x,tau_y,tau_z,tau_x_o,tau_y_o,tau_z_o, &
                            tau_x_acc,tau_y_acc,tau_z_acc)
#ifdef _IMPDIFF
      alphai = alpha**(-1)
      !
      !$OMP PARALLEL WORKSHARE
      up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END PARALLEL WORKSHARE
      call updt_rhs(lo,hiu,is_bound,rhsu%x,rhsu%y,rhsu%z,up)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call add_constant_to_n_diagonals(hiu_a(idir)-lo_a(idir)+1,lo_a(il:iu:iskip),hiu_a(il:iu:iskip), &
                                       alphai-alphaoi,usolver_fft(:)%mat) ! correct diagonal term
      call create_n_solvers(hiu_a(idir)-lo_a(idir)+1,hypre_maxiter,hypre_tol,hypre_solver_i,usolver_fft)
      call setup_n_solvers(hiu_a(idir)-lo_a(idir)+1,usolver_fft)
      call fft(arrplan_u(1),up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#ifndef _FFT_USE_SLABS
      call solve_n_helmholtz_2d(usolver_fft,lo(idir),hiu(idir),1,lo(il:iu:iskip),hiu(il:iu:iskip),up,uo)
#else
      call transpose_slab(1,0,t_params(:,1:2:1 ),comm_block,up,up_s)
      call solve_n_helmholtz_2d(usolver_fft,lo_a(idir),hiu_a(idir),0,lo_a(il:iu:iskip),hiu_a(il:iu:iskip),up_s,uo)
      call transpose_slab(0,1,t_params(:,2:1:-1),comm_block,up_s,up)
#endif
      call fft(arrplan_u(2),up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      !$OMP PARALLEL WORKSHARE
      up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_u
      !$OMP END PARALLEL WORKSHARE
      call finalize_n_solvers(hiu_a(idir)-lo_a(idir)+1,usolver_fft)
#else
      call add_constant_to_diagonal(lo,hiu,alphai-alphaoi,usolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,hypre_solver_i,usolver)
      call setup_solver(usolver)
      call solve_helmholtz(usolver,lo,hiu,up,uo)
      call finalize_solver(usolver)
#endif
      !
      !$OMP PARALLEL WORKSHARE
      vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END PARALLEL WORKSHARE
      call updt_rhs(lo,hiv,is_bound,rhsv%x,rhsv%y,rhsv%z,vp)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call add_constant_to_n_diagonals(hiv_a(idir)-lo_a(idir)+1,lo_a(il:iu:iskip),hiv_a(il:iu:iskip), &
                                       alphai-alphaoi,vsolver_fft(:)%mat) ! correct diagonal term
      call create_n_solvers(hiv_a(idir)-lo_a(idir)+1,hypre_maxiter,hypre_tol,hypre_solver_i,vsolver_fft)
      call setup_n_solvers(hiv_a(idir)-lo_a(idir)+1,vsolver_fft)
      call fft(arrplan_v(1),vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#ifndef _FFT_USE_SLABS
      call solve_n_helmholtz_2d(vsolver_fft,lo(idir),hiv(idir),1,lo(il:iu:iskip),hiv(il:iu:iskip),vp,vo)
#else
      call transpose_slab(1,0,t_params(:,1:2:1 ),comm_block,vp,vp_s)
      call solve_n_helmholtz_2d(vsolver_fft,lo_a(idir),hiv_a(idir),0,lo_a(il:iu:iskip),hiv_a(il:iu:iskip),vp_s,vo)
      call transpose_slab(0,1,t_params(:,2:1:-1),comm_block,vp_s,vp)
#endif
      call fft(arrplan_v(2),vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      !$OMP PARALLEL WORKSHARE
      vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_v
      !$OMP END PARALLEL WORKSHARE
      call finalize_n_solvers(hiv_a(idir)-lo_a(idir)+1,vsolver_fft)
#else
      call add_constant_to_diagonal(lo,hiv,alphai-alphaoi,vsolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,hypre_solver_i,vsolver)
      call setup_solver(vsolver)
      call solve_helmholtz(vsolver,lo,hiv,vp,vo)
      call finalize_solver(vsolver)
#endif
      !
      !$OMP PARALLEL WORKSHARE
      wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*alphai
      !$OMP END PARALLEL WORKSHARE
      call updt_rhs(lo,hiw,is_bound,rhsw%x,rhsw%y,rhsw%z,wp)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call add_constant_to_n_diagonals(hiw_a(idir)-lo_a(idir)+1,lo_a(il:iu:iskip),hiw_a(il:iu:iskip), &
                                       alphai-alphaoi,wsolver_fft(:)%mat) ! correct diagonal term
      call create_n_solvers(hiw_a(idir)-lo_a(idir)+1,hypre_maxiter,hypre_tol,hypre_solver_i,wsolver_fft)
      call setup_n_solvers(hiw_a(idir)-lo_a(idir)+1,wsolver_fft)
      call fft(arrplan_w(1),wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#ifndef _FFT_USE_SLABS
      call solve_n_helmholtz_2d(wsolver_fft,lo(idir),hiw(idir),1,lo(il:iu:iskip),hiw(il:iu:iskip),wp,wo)
#else
      call transpose_slab(1,0,t_params(:,1:2:1 ),comm_block,wp,wp_s)
      call solve_n_helmholtz_2d(wsolver_fft,lo_a(idir),hiw_a(idir),0,lo_a(il:iu:iskip),hiw_a(il:iu:iskip),wp_s,wo)
      call transpose_slab(0,1,t_params(:,2:1:-1),comm_block,wp_s,wp)
#endif
      call fft(arrplan_w(2),wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      !$OMP PARALLEL WORKSHARE
      wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_w
      !$OMP END PARALLEL WORKSHARE
      call finalize_n_solvers(hiw_a(idir)-lo_a(idir)+1,wsolver_fft)
#else
      call add_constant_to_diagonal(lo,hiw,alphai-alphaoi,wsolver%mat) ! correct diagonal term
      call create_solver(hypre_maxiter,hypre_tol,hypre_solver_i,wsolver)
      call setup_solver(wsolver)
      call solve_helmholtz(wsolver,lo,hiw,wp,wo)
      call finalize_solver(wsolver)
#endif
      !
      alphaoi = alphai
#endif
      call bounduvw(cbcvel,lo,hi,bcvel,.false.,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,up,vp,wp)
      call inflow(is_bound_inflow,.false.,lo,hi,u_in%x,v_in%x,w_in%x, &
                                                u_in%y,v_in%y,w_in%y, &
                                                u_in%z,v_in%z,w_in%z,up,vp,wp)
#if !defined(_IMPDIFF) && defined(_ONE_PRESS_CORR)
      dtrk  = dt
      if(irk < 3) then ! pressure correction only at the last RK step
        !$OMP PARALLEL WORKSHARE
        u(:,:,:) = up(:,:,:)
        v(:,:,:) = vp(:,:,:)
        w(:,:,:) = wp(:,:,:)
        !$OMP END PARALLEL WORKSHARE
        cycle
      end if
#endif
      call fillps(lo,hi,dxf,dyf,dzf,dtrk,up,vp,wp,pp)
      call updt_rhs(lo,hi,is_bound,rhsp%x,rhsp%y,rhsp%z,pp)
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
      call fft(arrplan_p(1),pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#ifndef _FFT_USE_SLABS
#ifdef _FFT_USE_SLICED_PENCILS
      call solve_n_helmholtz_3d(psolver_fft,nslices,lo_sp,hi_sp,1,lo,hi,pp,po)
#else
      call solve_n_helmholtz_2d(psolver_fft,lo(idir),hi(idir),1,lo(il:iu:iskip),hi(il:iu:iskip),pp,po)
#endif
#else
      call transpose_slab(1,0,t_params(:,1:2:1 ),comm_block,pp,pp_s)
      call solve_n_helmholtz_2d(psolver_fft,lo_s(idir),hi_s(idir),0,lo_s(il:iu:iskip),hi_s(il:iu:iskip),pp_s,po)
      call transpose_slab(0,1,t_params(:,2:1:-1),comm_block,pp_s,pp)
#endif
      call fft(arrplan_p(2),pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      !$OMP PARALLEL WORKSHARE
      pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = pp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft_p
      !$OMP END PARALLEL WORKSHARE
#else
      call solve_helmholtz(psolver,lo,hi,pp,po)
#endif
      alpha_bc_o(:,:) = alpha_bc(:,:)
      call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,pp)
      call correc(lo,hi,dxc,dyc,dzc,dtrk,pp,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,lo,hi,bcvel,.true.,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
      call inflow(is_bound_inflow,.true.,lo,hi,u_in%x,v_in%x,w_in%x, &
                                               u_in%y,v_in%y,w_in%y, &
                                               u_in%z,v_in%z,w_in%z,u,v,w)
      call updt_pressure(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,pp,p)
      call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
    end do
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep >= nstep   ) is_done = is_done.or..true.
    end if
    if(stop_type(2)) then ! maximum simulation time reached
      if(time  >= time_max) is_done = is_done.or..true.
    end if
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600._rp
      call MPI_ALLREDUCE(MPI_IN_PLACE,tw,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
      if(tw    >= tw_max  ) is_done = is_done.or..true.
    end if
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
      end if
      !
      call chkdiv(lo,hi,dxf,dyf,dzf,u,v,w,vol_all,MPI_COMM_WORLD,divtot,divmax)
      if(myid == 0) write(stdout,*) 'Total divergence = ', divtot, '| Maximum divergence = ', divmax
      if(divtot /= divtot) then!divmax > small.or.divtot /= divtot) then
        if(myid == 0) write(stderr,*) 'ERROR: maximum divergence is too large.'
        if(myid == 0) write(stderr,*) 'Aborting...'
        is_done = .true.
        kill = .true.
      end if
    end if
    !
    ! output routines below
    !
    if(mod(istep,iout0d) == 0) then
      !allocate(var(4))
      var(1) = 1._rp*istep
      var(2) = dt
      var(3) = time
      call out0d(trim(datadir)//'time.out',3,myid,var)
      include 'out0d.h90'
    end if
    write(fldnum,'(i7.7)') istep
    if(mod(istep,iout1d) == 0) then
      include 'out1d.h90'
    end if
    if(mod(istep,iout2d) == 0) then
      include 'out2d.h90'
    end if
    if(mod(istep,iout3d) == 0) then
      include 'out3d.h90'
    end if
    if(mod(istep,isave ) == 0.or.(is_done.and..not.kill)) then
      if(is_overwrite_save) then
        filename = 'fld_b_'//cblock//'.bin'
      else
        if(nsaves_max > 0) then
          savecounter = savecounter + 1
          ichkptnum = mod(savecounter,nsaves_max) + 1
          write(chkptnum,'(i4.4)') ichkptnum
          filename = 'fld_'//chkptnum//'_b_'//cblock//'.bin'
          var(1) = 1.*istep
          var(2) = time
          var(3) = 1.*ichkptnum
          call out0d(trim(datadir)//'log_saves.out',3,myid,var)
        end if
      end if
      call load('w',trim(datadir)//trim(filename),comm_block,ng,[1,1,1],lo_1,hi_1,u,v,w,p,po,time,istep)
      if(.not.is_overwrite_save) then
        !
        ! fld.bin -> last checkpoint file (symbolic link)
        !
        if(myid_block == 0) call execute_command_line('ln -sf '//trim(filename)//' '//trim(datadir)//'fld_b_'//cblock//'.bin')
      end if
      if(myid_block == 0) write(stdout,*) '*** Checkpoint saved at time = ', time, &
                                          'time step = ', istep, 'block = ',my_block,'. ***'
    end if
#ifdef _TIMING
      dt12 = MPI_WTIME()-dt12
      call MPI_ALLREDUCE(dt12,dt12av ,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD)
      call MPI_ALLREDUCE(dt12,dt12min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD)
      call MPI_ALLREDUCE(dt12,dt12max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
      if(myid == 0) write(stdout,*) 'Avrg, min & max elapsed time: '
      if(myid == 0) write(stdout,*) dt12av/(1._rp*nrank),dt12min,dt12max
#endif
  end do
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  call finalize_n_matrices(npsolvers,psolver_fft)
  call fftend(arrplan_p)
#else
  call finalize_matrix(psolver)
#endif
#ifdef _IMPDIFF
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  call finalize_n_matrices(hiu_a(idir)-lo_a(idir)+1,usolver_fft)
  call finalize_n_matrices(hiv_a(idir)-lo_a(idir)+1,vsolver_fft)
  call finalize_n_matrices(hiw_a(idir)-lo_a(idir)+1,wsolver_fft)
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
end program snac
