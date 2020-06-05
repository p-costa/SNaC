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
! Slow CaNS, a.k.a. **SNaC**
! Pedro Costa (p.simoes.costa@gmail.com)
!-----------------------------------------------------------------------------------
program snac
  use mpi
  use mod_bound          , only: bounduvw,boundp,updt_rhs
  use mod_chkdt          , only: chkdt
  use mod_common_mpi     , only: myid,ierr
  use mod_correc         , only: correc
  use mod_initflow       , only: initflow
  use mod_initgrid       , only: initgrid,distribute_grid,save_grid
  use mod_initmpi        , only: initmpi
  use mod_fillps         , only: fillps
  use mod_load           , only: load
  use mod_param          , only: read_input, &
                                 datadir, &
                                 rkcoeff, &
                                 ng,l,gt,gr,cfl,dtmin,uref,lref,rey,visc,            &
                                 inivel,is_wallturb,nstep,time_max,tw_max,stop_type, &
                                 restart,is_overwrite_save,                          &
                                 icheck,iout0d,iout1d,iout2d,iout3d,isave,           &
                                 cbcvel,bcvel,cbcpre,bcpre,                          &
                                 bforce, is_forced,velf,is_outflow,no_outflow,       &
                                 dims,nthreadsmax
  use mod_pressure_update, only: pressure_update
  use mod_rk             , only: rk_mom
  use mod_sanity         , only: test_sanity
  use mod_types
  !$ use omp_lib
  implicit none
  integer, dimension(0:1,3) :: nb
  logical, dimension(0:1,3) :: is_bound
  integer, dimension(3    ) :: halos
  integer , dimension(3) :: lo,hi
  real(rp), allocatable, dimension(:,:,:) :: u,v,w,p,up,vp,wp,pp,po
  real(rp), allocatable, dimension(:,:,:) :: dudtrko,dvdtrko,dwdtrko
  type rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type rhs_bound
  type(rhs_bound) :: rhsbp
  real(rp) :: alpha
#ifdef _IMPDIFF
  type(rhs_bound) :: rhsbu,rhsbv,rhsbw
#endif
  real(rp) :: dt,dtmax,time,dtrk,divtot,divmax
  integer  :: irk,istep
  real(rp), allocatable, dimension(:) :: dxc  ,dxf  ,xc  ,xf  , &
                                         dyc  ,dyf  ,yc  ,yf  , &
                                         dzc  ,dzf  ,zc  ,zf  , &
                                         dxc_g,dxf_g,xc_g,xf_g, &
                                         dyc_g,dyf_g,yc_g,yf_g, &
                                         dzc_g,dzf_g,zc_g,zf_g
  !
  real(rp) :: meanvelu,meanvelv,meanvelw
  real(rp), dimension(3) :: f,dpdl
  !
  real(rp), dimension(100) :: var
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
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,dims,cbcpre,lo,hi,nb,is_bound,halos)
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
           pp(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
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
  allocate(dxc_g(1-1:ng(1)+1), &
           dxf_g(1-1:ng(1)+1), &
            xc_g(1-1:ng(1)+1), &
            xf_g(1-1:ng(1)+1), &
           dyc_g(1-1:ng(2)+1), &
           dyf_g(1-1:ng(2)+1), &
            yc_g(1-1:ng(2)+1), &
            yf_g(1-1:ng(2)+1), &
           dzc_g(1-1:ng(3)+1), &
           dzf_g(1-1:ng(3)+1), &
            zc_g(1-1:ng(3)+1), &
            zf_g(1-1:ng(3)+1))
  allocate(rhsbp%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsbp%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsbp%z(lo(1):hi(1),lo(2):hi(2),0:1))
#ifdef _IMPDIFF
  allocate(rhsbu%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsbu%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsbu%z(lo(1):hi(1),lo(2):hi(2),0:1), &
           rhsbv%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsbv%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsbv%z(lo(1):hi(1),lo(2):hi(2),0:1), &
           rhsbw%x(lo(2):hi(2),lo(3):hi(3),0:1), &
           rhsbw%y(lo(1):hi(1),lo(3):hi(3),0:1), &
           rhsbw%z(lo(1):hi(1),lo(2):hi(2),0:1))
#endif
  !
  if(myid.eq.0) then
    write(stdout,*) '*******************************'
    write(stdout,*) '*** Beginning of simulation ***'
    write(stdout,*) '*******************************'
    write(stdout,*) ''
  endif
  !
  ! generate grid
  !
  call initgrid(ng(1),lo(1),hi(1),gt(1),gr(1),l(1), &
                dxc_g,dxf_g,xc_g,xf_g)
  call initgrid(ng(2),lo(2),hi(2),gt(2),gr(2),l(2), &
                dyc_g,dyf_g,yc_g,yf_g)
  call initgrid(ng(3),lo(3),hi(3),gt(3),gr(3),l(3), &
                dzc_g,dzf_g,zc_g,zf_g)
  call save_grid(trim(datadir)//'grid_x',ng(1),xf_g,xc_g,dxf_g,dxc_g)
  call save_grid(trim(datadir)//'grid_y',ng(2),yf_g,yc_g,dyf_g,dyc_g)
  call save_grid(trim(datadir)//'grid_z',ng(3),zf_g,zc_g,dzf_g,dzc_g)
  open(newunit=iunit,status='replace',file=trim(datadir)//'geometry.out')
    write(iunit,*) ng(1),ng(2),ng(3) 
    write(iunit,*) l(1),l(2),l(3) 
  close(iunit)
  call distribute_grid(lo(1),hi(1),dxc_g,dxc)
  call distribute_grid(lo(1),hi(1),dxf_g,dxf)
  call distribute_grid(lo(1),hi(1), xc_g, xc)
  call distribute_grid(lo(1),hi(1), xf_g, xf)
  call distribute_grid(lo(2),hi(2),dyc_g,dyc)
  call distribute_grid(lo(2),hi(2),dyf_g,dyf)
  call distribute_grid(lo(2),hi(2), yc_g, yc)
  call distribute_grid(lo(2),hi(2), yf_g, yf)
  call distribute_grid(lo(3),hi(3),dzc_g,dzc)
  call distribute_grid(lo(3),hi(3),dzf_g,dzf)
  call distribute_grid(lo(3),hi(3), zc_g, zc)
  call distribute_grid(lo(3),hi(3), zf_g, zf)
  !
  ! initialization of the flow fields
  !
  if(.not.restart) then
    istep = 0
    time = 0._rp
    call initflow(inivel,is_wallturb,lo,hi,ng,l,uref,lref,visc,bforce(1), &
                  xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,p)
    if(myid.eq.0) write(stdout,*) '*** Initial condition succesfully set ***'
  else
    call load('r',trim(datadir)//'fld.bin',ng,[1,1,1],lo,hi,u,v,w,p,time,istep)
    if(myid.eq.0) write(stdout,*) '*** Checkpoint loaded at time = ', time, 'time step = ', istep, '. ***'
  endif
  call bounduvw(cbcvel,lo,hi,bcvel,no_outflow,halos,is_bound,nb, &
                dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
  call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
  up(:,:,:)      = 0._rp
  vp(:,:,:)      = 0._rp
  wp(:,:,:)      = 0._rp
  pp(:,:,:)      = 0._rp
  dudtrko(:,:,:) = 0._rp
  dvdtrko(:,:,:) = 0._rp
  dwdtrko(:,:,:) = 0._rp
  !
  ! post-process and write initial condition
  !
  write(fldnum,'(i7.7)') istep
  !include 'out1d.h90'
  !include 'out3d.h90'
  !
  ! determine time step
  !
  call chkdt(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,visc,u,v,w,dtmax)
  dt = min(cfl*dtmax,dtmin)
  if(myid.eq.0) write(stdout,*) 'dtmax = ', dtmax, 'dt = ',dt
  !
  ! main loop
  !
  if(myid.eq.0) write(stdout,*) '*** Calculation loop starts now ***'
  kill    = .false.
  is_done = .false.
  do while(.not.is_done)
#ifdef _TIMING
    dt12 = MPI_WTIME()
#endif
    istep = istep + 1
    time  = time  + dt
    if(myid.eq.0) write(stdout,*) 'Timestep #', istep, 'Time = ', time
    dpdl(:)  = 0._rp
    do irk=1,3
      dtrk = sum(rkcoeff(:,irk))*dt
      alpha = visc*dtrk/2._rp
      call rk_mom(rkcoeff(:,irk),lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,l,dt,bforce, &
                  is_forced,velf,visc,u,v,w,p,dudtrko,dvdtrko,dwdtrko,up,vp,wp,f)
      dpdl(:) = dpdl(:) - f(:)/dt
      call bounduvw(cbcvel,lo,hi,bcvel,no_outflow,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
      call fillps(lo,hi,dxf,dyf,dzf,dtrk,up,vp,wp,p)
      call updt_rhs(cbcpre,lo,hi,is_bound,rhsbp%x,rhsbp%y,rhsbp%z,pp)
      !solver
      call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
      call correc(lo,hi,dxc,dyc,dzc,dtrk,p,up,vp,wp,u,v,w)
      call bounduvw(cbcvel,lo,hi,bcvel,is_outflow,halos,is_bound,nb, &
                    dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
      call pressure_update(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,pp,p)
      call boundp(  cbcpre,lo,hi,bcpre,halos,is_bound,nb,dxc,dyc,dzc,p)
    enddo
    !
    ! check simulation stopping criteria
    !
    if(stop_type(1)) then ! maximum number of time steps reached
      if(istep.ge.nstep   ) is_done = is_done.or..true.
    endif
    if(stop_type(2)) then ! maximum simulation time reached
      if(time .ge.time_max) is_done = is_done.or..true.
    endif
    if(stop_type(3)) then ! maximum wall-clock time reached
      tw = (MPI_WTIME()-twi)/3600._rp
      call MPI_ALLREDUCE(MPI_IN_PLACE,tw,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
      if(tw   .ge.tw_max  ) is_done = is_done.or..true.
    endif
  enddo
  if(myid.eq.0.and.(.not.kill)) write(stdout,*) '*** Fim ***'
  call MPI_FINALIZE(ierr)
  call exit
end program snac
