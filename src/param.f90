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
integer  :: hypre_solver_i
real(rp) :: hypre_tol
integer  :: hypre_maxiter
!
! parameters to be determined from the input file 'dns.in'
!
real(rp)               :: cfl,dtmin
real(rp)               :: uref,lref,rey,visc
integer                :: nstep
real(rp)               :: time_max,tw_max
logical , dimension(3) :: stop_type
logical                :: restart,is_overwrite_save
integer                :: nsaves_max
integer                :: icheck,iout0d,iout1d,iout2d,iout3d,isave
real(rp), dimension(3) :: bforce
integer                :: nthreadsmax
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
  use mpi_f08
  use mod_common_mpi, only:myid
  implicit none
  integer :: iunit,iblock,ierr
  integer, allocatable, dimension(:) :: nranks
  logical :: exists
  character(len=100) :: filename
  character(len=  3) :: cblock
  integer :: q
  nsaves_max = 0 ! a good default, for backward compatibility
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*) cfl,dtmin
        read(iunit,*) uref,lref,rey
        read(iunit,*) nstep, time_max,tw_max
        read(iunit,*) stop_type(1),stop_type(2),stop_type(3)
        read(iunit,*) restart,is_overwrite_save,nsaves_max
        read(iunit,*) icheck,iout0d,iout1d,iout2d,iout3d,isave
        read(iunit,*)  bforce(1),bforce(2),bforce(3)
        read(iunit,*) nthreadsmax
      else
        if(myid == 0) write(stderr,*) '*** Error reading the input file *** '
        if(myid == 0) write(stderr,*) 'Aborting...'
        call MPI_FINALIZE()
        error stop
      end if
    close(iunit)
    visc = uref*lref/rey
    !
    exists = .true.
    nblocks = 1
    filename = 'geo/block.'
    do while(exists)
      write(cblock,'(i3.3)') nblocks
      inquire(file=trim(filename)//cblock, exist = exists)
      if(exists) nblocks = nblocks + 1
    end do
    nblocks = nblocks - 1
    allocate(nranks(nblocks))
    do iblock = 1,nblocks
      write(cblock,'(i3.3)') iblock
      open(newunit=iunit,file=trim(filename)//cblock,status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
          read(iunit,*) dims(1),dims(2),dims(3)
          nranks(iblock) = product(dims(:))
        else
          if(myid == 0) write(stderr,*) '*** Error reading the input file *** '
          if(myid == 0) write(stderr,*) 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      close(iunit)
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
      write(cblock,'(i3.3)') iblock
      open(newunit=iunit,file=trim(filename)//cblock,status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
          if( myid >= sum(nranks(1:iblock-1)).and. &
              myid <  sum(nranks(1:iblock  )) ) then
            read(iunit,*) dims(1),dims(2),dims(3)
            read(iunit,*) ng(1),ng(2),ng(3)
            read(iunit,*) lmin(1),lmin(2),lmin(3)
            read(iunit,*) lmax(1),lmax(2),lmax(3)
            read(iunit,*) gt(1),gt(2),gt(3)
            read(iunit,*) gr(1),gr(2),gr(3)
            read(iunit,*) cbcvel(0,1,1),cbcvel(1,1,1),cbcvel(0,2,1),cbcvel(1,2,1),cbcvel(0,3,1),cbcvel(1,3,1)
            read(iunit,*) cbcvel(0,1,2),cbcvel(1,1,2),cbcvel(0,2,2),cbcvel(1,2,2),cbcvel(0,3,2),cbcvel(1,3,2)
            read(iunit,*) cbcvel(0,1,3),cbcvel(1,1,3),cbcvel(0,2,3),cbcvel(1,2,3),cbcvel(0,3,3),cbcvel(1,3,3)
            read(iunit,*) cbcpre(0,1  ),cbcpre(1,1  ),cbcpre(0,2  ),cbcpre(1,2  ),cbcpre(0,3  ),cbcpre(1,3  )
            read(iunit,*)  bcvel(0,1,1), bcvel(1,1,1), bcvel(0,2,1), bcvel(1,2,1), bcvel(0,3,1), bcvel(1,3,1)
            read(iunit,*)  bcvel(0,1,2), bcvel(1,1,2), bcvel(0,2,2), bcvel(1,2,2), bcvel(0,3,2), bcvel(1,3,2)
            read(iunit,*)  bcvel(0,1,3), bcvel(1,1,3), bcvel(0,2,3), bcvel(1,2,3), bcvel(0,3,3), bcvel(1,3,3)
            read(iunit,*)  bcpre(0,1  ), bcpre(1,1  ), bcpre(0,2  ), bcpre(1,2  ), bcpre(0,3  ), bcpre(1,3  )
            read(iunit,*)  inflow_type(0,1),inflow_type(1,1),inflow_type(0,2),inflow_type(1,2),inflow_type(0,3),inflow_type(1,3)
            read(iunit,*) inivel
            my_block = iblock
            id_first = sum(nranks(1:iblock-1))
            do q=1,3
              is_periodic(q) = is_periodic(q).and.(cbcpre(0,q)//cbcpre(1,q) == 'FF')
            end do
          end if
        else
          if(myid == 0) write(stderr,*) '*** Error reading the input file *** '
          if(myid == 0) write(stderr,*) 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        end if
      close(iunit)
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
    ! read iterative solver parameter file hypre.in, if it exists
    !
    inquire(file='hypre.in', exist = exists)
    if(exists) then
      open(newunit=iunit,file='hypre.in',status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
          read(iunit,*) hypre_solver_i,hypre_tol,hypre_maxiter
          if(hypre_solver_i < 1 .or. hypre_solver_i > 4) then
            if(myid == 0) write(stderr,*) '*** Error: invalid solver choice [1-4] *** '
            if(myid == 0) write(stderr,*) 'Reverting to the default (2 -> PFMG)...'
            hypre_solver_i = hypre_solver_i_default
          end if
          if(hypre_tol > 1._rp .or. hypre_tol < 0._rp) then
            if(myid == 0) write(stderr,*) '*** Error: iterative error tolerance is too high or negative *** '
            if(myid == 0) write(stderr,*) 'Reverting to the default (1.e-4)...'
            hypre_tol     = hypre_tol_default
          end if
          if(hypre_maxiter < 0) then
            if(myid == 0) write(stderr,*) '*** Error: maximum number of iterations needs to be > 0 *** '
            if(myid == 0) write(stderr,*) 'Reverting to the default (50)...'
            hypre_maxiter = hypre_maxiter_default
          end if
        else
          if(myid == 0) write(stderr,*) '*** Error reading the input file *** '
          if(myid == 0) write(stderr,*) 'Aborting...'
          call MPI_FINALIZE()
          error stop
        end if
      close(iunit)
    else
      !
      ! defaults
      !
      hypre_solver_i = hypre_solver_i_default
      hypre_tol      = hypre_tol_default
      hypre_maxiter  = hypre_maxiter_default
    end if
  end subroutine read_input
end module mod_param
