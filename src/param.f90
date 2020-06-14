module mod_param
use mod_types
implicit none
public
!
! parameters
!
real(rp), parameter :: pi = acos(-1._rp)
real(rp), parameter :: small = epsilon(pi)*10._rp**(precision(pi)/2._rp)
logical , parameter, dimension(2,3) :: no_outflow = & 
    reshape([.false.,.false.,   & ! no outflow in x lower,upper bound
             .false.,.false.,   & ! no outflow in y lower,upper bound
             .false.,.false.], & ! no outflow in z lower,upper bound
              shape(no_outflow))
character(len=100), parameter :: datadir = 'data/'
real(rp), parameter, dimension(2,3) :: rkcoeff = reshape( [32._rp/60._rp,  0._rp        , &
                                                           25._rp/60._rp, -17._rp/60._rp, &
                                                           45._rp/60._rp, -25._rp/60._rp], shape(rkcoeff))
real(rp), parameter, dimension(3)   :: rkcoeff12 = rkcoeff(1,:)+rkcoeff(2,:)
!
! parameters to be determined from the input file 'dns.in'
!
real(rp)                             :: cfl,dtmin
real(rp)                             :: uref,lref,rey,visc
integer                              :: nstep
real(rp)                             :: time_max,tw_max
logical , dimension(3)               :: stop_type
logical                              :: restart,is_overwrite_save
integer                              :: icheck,iout0d,iout1d,iout2d,iout3d,isave
real(rp), dimension(3)               :: bforce
integer                              :: nthreadsmax
!
! parameters specific to each block
!
integer           , dimension(      3) :: dims
integer           , dimension(      3) :: lo,hi
real(rp)          , dimension(      3) :: lmin,lmax
integer           , dimension(      3) :: gt
real(rp)          , dimension(      3) :: gr
character(len=1  ), dimension(0:1,3,3) :: cbcvel
real(rp)          , dimension(0:1,3,3) ::  bcvel
character(len=1  ), dimension(0:1,  3) :: cbcpre
real(rp)          , dimension(0:1,  3) ::  bcpre
logical           , dimension(0:1,  3) :: is_outflow
integer           , dimension(      3) :: periods
character(len=100)                     :: inivel
!
real(rp) :: vol_all
integer  :: my_block
contains 
  subroutine read_input()
  use mpi
  use mod_common_mpi, only:myid,ierr
  implicit none
  integer :: iunit,iblock,nblocks,maxrank
  logical :: exists
  character(len=100) :: filename
  character(len=  3) :: cblock
    open(newunit=iunit,file='dns.in',status='old',action='read',iostat=ierr)
      if( ierr == 0 ) then
        read(iunit,*) cfl,dtmin
        read(iunit,*) uref,lref,rey
        read(iunit,*) nstep, time_max,tw_max
        read(iunit,*) stop_type(1),stop_type(2),stop_type(3)
        read(iunit,*) restart,is_overwrite_save
        read(iunit,*) icheck,iout0d,iout1d,iout2d,iout3d,isave
        read(iunit,*)  bforce(1),bforce(2),bforce(3)
        read(iunit,*) nthreadsmax
      else
        if(myid == 0) write(stderr,*) '*** Error reading the input file *** ' 
        if(myid == 0) write(stderr,*) 'Aborting...'
        call MPI_FINALIZE(ierr)
        error stop
      endif
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
    enddo
    nblocks = nblocks - 1
    maxrank = 0
    do iblock = 1,nblocks
      write(cblock,'(i3.3)') iblock
      open(newunit=iunit,file=trim(filename)//cblock,status='old',action='read',iostat=ierr)
        if( ierr == 0 ) then
          read(iunit,*) dims(1),dims(2),dims(3)
          maxrank = maxrank + product(dims(:))-1
          !
          ! WRONG: NEED AN MPI_SCAN HERE TO DETERMINE THE EXTENTS 
          !
          if(myid <= maxrank) then
            read(iunit,*) lo(1),lo(2),lo(3)
            read(iunit,*) hi(1),hi(2),hi(3)
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
            read(iunit,*) periods(1),periods(2),periods(3)
            read(iunit,*) is_outflow(0,1),is_outflow(1,1),is_outflow(0,2),is_outflow(1,2),is_outflow(0,3),is_outflow(1,3)
            read(iunit,*) inivel
            my_block = iblock
          endif
        else
          if(myid == 0) write(stderr,*) '*** Error reading the input file *** ' 
          if(myid == 0) write(stderr,*) 'Aborting...'
          call MPI_FINALIZE(ierr)
          error stop
        endif
      close(iunit)
    enddo
    !
    ! compute volume of all blocks (useful to compute bulk averages)
    !
    vol_all = product(lmax(:)-lmin(:))/product(dims(:))
    call mpi_allreduce(MPI_IN_PLACE,vol_all,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    !
    print*,'delete later; total volume = ', vol_all
  return
  end subroutine read_input
end module mod_param
