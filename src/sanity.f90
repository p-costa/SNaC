module mod_sanity
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_types
  implicit none
  private
  public test_sanity
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  public test_sanity_fft
#endif
  contains
  subroutine test_sanity(gr,stop_type,cbcvel,cbcpre)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    real(rp)        , intent(in), dimension(3      ) :: gr
    logical         , intent(in), dimension(3      ) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
    logical :: passed
    !
    call chk_grid(gr,passed);                    if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);        if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,passed);           if(.not.passed) call abortit
  end subroutine test_sanity
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  subroutine test_sanity_fft(dims,lo,hi,lmin,lmax,gr)
    integer , intent(in ) :: dims,lo,hi
    real(rp), intent(in ) :: lmin,lmax,gr
    logical :: passed
    call chk_grid_fft(dims,lo,hi,lmin,lmax,gr,passed); if(.not.passed) call abortit
  end subroutine test_sanity_fft
#endif
  !
  subroutine chk_stop_type(stop_type,passed)
    implicit none
    logical, intent(in ), dimension(3) :: stop_type
    logical, intent(out) :: passed
    logical :: passed_loc
    passed = .true.
    passed_loc = any(stop_type(:))
    if(.not.passed_loc) &
      call write_error('stopping criterion not chosen.')
    passed = passed.and.passed_loc
  end subroutine chk_stop_type
  !
  subroutine chk_grid(gr,passed)
    implicit none
    real(rp), intent(in ), dimension(3) :: gr
    logical , intent(out) :: passed
    logical :: passed_loc
    passed = .true.
    passed_loc = all(gr(:) >= 0._rp)
    if(.not.passed_loc) & 
      call write_error('grid growth parameter must be positive.')
    passed = passed.and.passed_loc
  end subroutine chk_grid
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  subroutine chk_grid_fft(dims,lo,hi,lmin,lmax,gr,passed)
    implicit none
    integer , intent(in ) :: dims,lo,hi
    real(rp), intent(in ) :: lmin,lmax,gr
    logical , intent(out) :: passed
    logical :: passed_loc
    integer  :: lo_min  ,lo_max  ,hi_min  ,hi_max
    real(rp) :: lmin_min,lmin_max,lmax_min,lmax_max
    passed = .true.
    passed_loc = dims == 1
    call MPI_ALLREDUCE(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,ierr)
    if(.not.passed_loc) &
      call write_error('no domain decomposition allowed in the uniform (FFT) direction')
    passed = passed.and.passed_loc
    !
    call MPI_ALLREDUCE(lo  ,lo_min  ,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(lo  ,lo_max  ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(hi  ,hi_min  ,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(hi  ,hi_max  ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(lmin,lmin_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(lmin,lmin_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(lmax,lmax_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD,ierr) 
    call MPI_ALLREDUCE(lmax,lmax_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr) 
    passed_loc = (lo_min == lo).and.(lo_max == lo).and. &
                 (hi_min == hi).and.(hi_max == hi)
             print*,lo_min,lo,lo_max,lo,hi_min,hi
    passed_loc = passed_loc.and. &
                 (lmin_min == lmin).and.(lmin_max == lmin).and. &
                 (lmax_min == lmax).and.(lmax_max == lmax)
    if(.not.passed_loc) &
      call write_error('domain size and number of points cannot vary among blocks in the FFT direction')
    passed = passed.and.passed_loc
    !
    passed_loc = gr == 0._rp
    if(.not.passed_loc) &
      call write_error('grid growth parameter in the FFT direction must be zero.')
    passed = passed.and.passed_loc
  end subroutine chk_grid_fft
#endif
  !
  subroutine chk_bc(cbcvel,cbcpre,passed)
    implicit none
    character(len=2), parameter, dimension(9) :: bcs = ['ND', 'DN', 'NN', 'DD', &
                                                        'FD', 'DF', 'FF', 'FN', 'NF']
    character(len=1), intent(in ), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in ), dimension(0:1,3  ) :: cbcpre
    logical         , intent(out) :: passed
    character(len=2) :: bc01v,bc01p
    integer :: ivel,idir
    logical :: passed_loc
    passed = .true.
    !
    ! check validity of pressure and velocity BCs
    !
    passed_loc = .true.
    do ivel = 1,3
      do idir=1,3
        bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
        passed_loc = passed_loc.and.any(bc01v == bcs)
      enddo
    enddo
    if(.not.passed_loc) call write_error('velocity BCs not valid.')
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.any(bc01p == bcs)
    enddo
    if(.not.passed_loc) call write_error('pressure BCs not valid.')
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'FF'.and.bc01p == 'FF').or. &
                                    (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                    (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                    (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                    (bc01v == 'FD'.and.bc01p == 'FN').or. &
                                    (bc01v == 'DF'.and.bc01p == 'NF').or. &
                                    (bc01v == 'FN'.and.bc01p == 'FD').or. &
                                    (bc01v == 'NF'.and.bc01p == 'DF').or. &
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    enddo
    if(.not.passed_loc) call write_error('velocity and pressure BCs not compatible.')
    passed = passed.and.passed_loc
    !
  end subroutine chk_bc
  !
  subroutine abortit
    implicit none
    if(myid == 0) write(stderr,*) ''
    if(myid == 0) write(stderr,*) '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) write(stderr,*) '    check dns.in'
    call MPI_FINALIZE(ierr)
    error stop
  end subroutine abortit
  subroutine write_error(message)
    character(len=*), intent(in) :: message
    if(myid == 0) write(stderr,*) 'ERROR: '//message
  end subroutine write_error
end module mod_sanity
