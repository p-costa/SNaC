module mod_param
use mod_types
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
real(rp), parameter :: small = epsilon(pi)*10._rp**(precision(pi)/2._rp)
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape( [32._rp/60._rp,  0._rp        , &
                                                           25._rp/60._rp, -17._rp/60._rp, &
                                                           45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
integer , parameter :: hypre_solver_i_default  = 2
real(rp), parameter :: hypre_tol_default     = real(1.e-4,rp)
integer , parameter :: hypre_maxiter_default = 50
integer , parameter :: max_blocks = 999
integer  :: hypre_solver_i
real(rp) :: hypre_tol
integer  :: hypre_maxiter
!
! parameters to be determined from the input file 'dns.nml'
!
real(rp)               :: cfl,dtmax,dt_f
real(rp)               :: uref,lref,rey,visc
integer                :: nstep
real(rp)               :: time_max,tw_max
logical , dimension(3) :: stop_type
logical                :: restart,is_overwrite_save
integer                :: nsaves_max
integer                :: icheck,iout0d,iout1d,iout2d,iout3d,isave
real(rp), dimension(3) :: bforce
!
! parameters specific to each block
!
integer           , dimension(      3) :: dims
integer           , dimension(      3) :: ng
real(rp)          , dimension(      3) :: lmin,lmax
integer           , dimension(      3) :: gt
real(rp)          , dimension(      3) :: gr
character(len=1  ), dimension(0:1,3,3) :: cbcvel
real(rp)          , dimension(0:1,3,3) ::  bcvel
character(len=1  ), dimension(0:1,  3) :: cbcpre
real(rp)          , dimension(0:1,  3) ::  bcpre
integer           , dimension(0:1,  3) ::  inflow_type
character(len=100)                     :: inivel
!
real(rp) :: vol_all
integer  :: my_block,id_first,nblocks,nrank
logical , dimension(3) :: is_periodic
integer , dimension(3) :: periods
real(rp), dimension(3) :: l_periodic
real(rp), dimension(3) :: lmin_min,lmax_max
contains
  subroutine read_input()
  use, intrinsic :: iso_fortran_env, only: iostat_end
  use mpi_f08
  use mod_common_mpi, only:myid
  implicit none
  integer :: iunit,iblock,ierr
  integer, allocatable, dimension(:) :: nranks
  integer :: q
  character(len=1024) :: iomsg
  character(len=2) :: cbc_pair
  integer           , dimension(      3,max_blocks) :: block_dims,block_ng,block_gt
  real(rp)          , dimension(      3,max_blocks) :: block_lmin,block_lmax,block_gr
  character(len=1  ), dimension(0:1,3,3,max_blocks) :: block_cbcvel
  real(rp)          , dimension(0:1,3,3,max_blocks) :: block_bcvel
  character(len=1  ), dimension(0:1,3,  max_blocks) :: block_cbcpre
  real(rp)          , dimension(0:1,3,  max_blocks) :: block_bcpre
  integer           , dimension(0:1,3,  max_blocks) :: block_inflow_type
  character(len=100), dimension(        max_blocks) :: block_inivel
  namelist /dns/ cfl,dtmax,dt_f, &
                 uref,lref,rey, &
                 nstep,time_max,tw_max, &
                 stop_type, &
                 restart,is_overwrite_save,nsaves_max, &
                 icheck,iout0d,iout1d,iout2d,iout3d,isave, &
                 bforce
  namelist /hypre/ hypre_solver_i,hypre_tol,hypre_maxiter
  namelist /blocks/ nblocks, &
                    block_dims,block_ng,block_lmin,block_lmax,block_gt,block_gr, &
                    block_cbcvel,block_cbcpre,block_bcvel,block_bcpre, &
                    block_inflow_type,block_inivel
  !
  ! defaults
  !
  nsaves_max = 0 ! a good default, for backward compatibility
  dt_f = -1._rp
  hypre_solver_i = hypre_solver_i_default
  hypre_tol      = hypre_tol_default
  hypre_maxiter  = hypre_maxiter_default
  block_dims(:,:) = 1
  block_ng(:,:) = 1
  block_lmin(:,:) = 0._rp
  block_lmax(:,:) = 1._rp
  block_gt(:,:) = 0
  block_gr(:,:) = 0._rp
  block_cbcvel(:,:,:,:) = 'D'
  block_cbcpre(:,:,:) = 'N'
  block_bcvel(:,:,:,:) = 0._rp
  block_bcpre(:,:,:) = 0._rp
  block_inflow_type(:,:,:) = 0
  block_inivel(:) = 'zer'
  nblocks = 0
  !
  ! read global parameters
  !
  open(newunit=iunit,file='dns.nml',status='old',action='read',iostat=ierr,iomsg=iomsg)
    if(ierr /= 0) call abort_input('dns.nml',iomsg)
    rewind(iunit)
    read(iunit,nml=dns,iostat=ierr,iomsg=iomsg)
    if(ierr /= 0) call abort_input('dns namelist',iomsg)
    rewind(iunit)
    read(iunit,nml=hypre,iostat=ierr,iomsg=iomsg)
    if(ierr /= 0 .and. ierr /= iostat_end) &
      call abort_input('hypre namelist',iomsg)
  close(iunit)
  visc = uref*lref/rey
  !
  ! read block-specific parameters
  !
  open(newunit=iunit,file='blocks.nml',status='old',action='read',iostat=ierr,iomsg=iomsg)
    if(ierr /= 0) call abort_input('blocks.nml',iomsg)
    read(iunit,nml=blocks,iostat=ierr,iomsg=iomsg)
    if(ierr /= 0) call abort_input('blocks namelist',iomsg)
  close(iunit)
  if(nblocks < 1 .or. nblocks > max_blocks) then
    if(myid == 0) write(stderr,*) '*** Error: invalid number of blocks in blocks.nml ***'
    if(myid == 0) write(stderr,*) 'Aborting...'
    call MPI_FINALIZE(ierr)
    error stop
  end if
  allocate(nranks(nblocks))
  do iblock = 1,nblocks
    nranks(iblock) = product(block_dims(:,iblock))
  end do
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
  if(nrank /= sum(nranks(1:nblocks))) then
    if(myid == 0) write(stderr,*) '*** Error: invalid number of MPI tasks. ***'
    if(myid == 0) write(stderr,*) 'Expected: ',sum(nranks(1:nblocks)), ' Found: ', nrank
    if(myid == 0) write(stderr,*) 'Aborting...'
    call MPI_FINALIZE(ierr)
    error stop
  end if
  is_periodic(:) = .true.
  do iblock = 1,nblocks
    if( myid >= sum(nranks(1:iblock-1)).and. &
        myid <  sum(nranks(1:iblock  )) ) then
      dims(:) = block_dims(:,iblock)
      ng(:) = block_ng(:,iblock)
      lmin(:) = block_lmin(:,iblock)
      lmax(:) = block_lmax(:,iblock)
      gt(:) = block_gt(:,iblock)
      gr(:) = block_gr(:,iblock)
      cbcvel(:,:,:) = block_cbcvel(:,:,:,iblock)
      cbcpre(:,:) = block_cbcpre(:,:,iblock)
      bcvel(:,:,:) = block_bcvel(:,:,:,iblock)
      bcpre(:,:) = block_bcpre(:,:,iblock)
      inflow_type(:,:) = block_inflow_type(:,:,iblock)
      inivel = block_inivel(iblock)
      my_block = iblock
      id_first = sum(nranks(1:iblock-1))
    end if
    do q=1,3
      cbc_pair = block_cbcpre(0,q,iblock)//block_cbcpre(1,q,iblock)
      is_periodic(q) = is_periodic(q).and.(cbc_pair == 'FF' .or. cbc_pair == 'PP')
    end do
  end do
  call mpi_allreduce(MPI_IN_PLACE,is_periodic,3,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
  deallocate(nranks)
  !
  ! compute volume of all blocks (useful to compute bulk averages)
  !
  vol_all = product(lmax(:)-lmin(:))/product(dims(:))
  call mpi_allreduce(MPI_IN_PLACE,vol_all,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  !
  ! determine length of the domain in the periodic directions
  !
  call MPI_ALLREDUCE(lmin(1),lmin_min(1),3,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(lmax(1),lmax_max(1),3,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
  where(is_periodic(:))
    l_periodic(:) = lmax_max(:)-lmin_min(:)
  elsewhere
    l_periodic(:) = 0._rp
  end where
  !
  ! validate iterative solver parameters
  !
  if(hypre_solver_i < 1 .or. hypre_solver_i > 4) then
    if(myid == 0) write(stderr,*) '*** Error: invalid solver choice [1-4] *** '
    if(myid == 0) write(stderr,*) 'Reverting to the default (2 -> PFMG)...'
    hypre_solver_i = hypre_solver_i_default
  end if
  if(hypre_tol > 1._rp .or. hypre_tol < 0._rp) then
    if(myid == 0) write(stderr,*) '*** Error: iterative error tolerance is too high or negative *** '
    if(myid == 0) write(stderr,*) 'Reverting to the default (1.e-4)...'
    hypre_tol = hypre_tol_default
  end if
  if(hypre_maxiter < 0) then
    if(myid == 0) write(stderr,*) '*** Error: maximum number of iterations needs to be > 0 *** '
    if(myid == 0) write(stderr,*) 'Reverting to the default (50)...'
    hypre_maxiter = hypre_maxiter_default
  end if
  contains
    subroutine abort_input(fname,msg)
      character(len=*), intent(in) :: fname,msg
      if(myid == 0) write(stderr,*) '*** Error reading ', trim(fname), ': ', trim(msg)
      if(myid == 0) write(stderr,*) 'Aborting...'
      call MPI_FINALIZE(ierr)
      error stop
    end subroutine abort_input
  end subroutine read_input
end module mod_param
