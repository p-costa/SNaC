module mod_sanity
  use mpi_f08
  use mod_common_mpi, only: myid
  use mod_types
  implicit none
  private
  public test_sanity
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  public test_sanity_fft
#endif
  contains
  subroutine test_sanity(lo,hi,dims,gr,stop_type,cbcvel,cbcpre,periods,inflow_type,outflow_type)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer         , intent(in), dimension(3      ) :: lo,hi,dims
    real(rp)        , intent(in), dimension(3      ) :: gr
    logical         , intent(in), dimension(3      ) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
    integer         , intent(in), dimension(    3  ) :: periods
    integer         , intent(in), dimension(0:1,3  ) :: inflow_type,outflow_type
    logical :: passed
    !
    call chk_grid(gr,passed);             if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed); if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,passed);    if(.not.passed) call abortit
    call chk_dims(lo,hi,dims,passed);     if(.not.passed) call abortit
    call chk_inoutflow(periods,inflow_type,outflow_type,passed); if(.not.passed) call abortit
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
    call mpi_allreduce(MPI_IN_PLACE,passed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
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
    call MPI_ALLREDUCE(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
    if(.not.passed_loc) &
      call write_error('no domain decomposition allowed in the uniform (FFT) direction')
    passed = passed.and.passed_loc
    !
    call MPI_ALLREDUCE(lo  ,lo_min  ,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(lo  ,lo_max  ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(hi  ,hi_min  ,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(hi  ,hi_max  ,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(lmin,lmin_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(lmin,lmin_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(lmax,lmax_min,1,MPI_REAL_RP,MPI_MIN,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(lmax,lmax_max,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
    passed_loc = (lo_min == lo).and.(lo_max == lo).and. &
                 (hi_min == hi).and.(hi_max == hi)
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
    call mpi_allreduce(MPI_IN_PLACE,passed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
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
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
    if(.not.passed_loc) call write_error('velocity BCs not valid.')
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.any(bc01p == bcs)
    end do
    call mpi_allreduce(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
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
                                    (bc01v == 'FN'.and.bc01p == 'FN').or. & ! for inflow/outflow
                                    (bc01v == 'NF'.and.bc01p == 'NF').or. & ! for inflow/outflow
                                    (bc01v == 'NN'.and.bc01p == 'NN').or. & ! for inflow/outflow
                                    (bc01v == 'DN'.and.bc01p == 'NN').or. & ! for inflow/outflow
                                    (bc01v == 'ND'.and.bc01p == 'NN').or. & ! for inflow/outflow
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    end do
    call mpi_allreduce(MPI_IN_PLACE,passed_loc,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
    if(.not.passed_loc) call write_error('velocity and pressure BCs not compatible.')
    passed = passed.and.passed_loc
    !
  end subroutine chk_bc
  !
  subroutine chk_dims(lo,hi,dims,passed)
    implicit none
    integer, intent(in), dimension(3) :: lo,hi,dims
    logical         , intent(out) :: passed
    passed = .true.
    passed = passed.and.(all(dims(:) <= hi(:)-lo(:)+1))
    if(.not.passed) call write_error('MPI task partitions cannot exceed the number of grid points.')
    call mpi_allreduce(MPI_IN_PLACE,passed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
  end subroutine chk_dims
  !
  subroutine chk_inoutflow(periods,inflow_type,outflow_type,passed)
    implicit none
    integer, intent(in ), dimension(    3) :: periods
    integer, intent(in ), dimension(0:1,3) :: inflow_type,outflow_type
    logical, intent(out)                   :: passed
    logical,              dimension(0:1,3) :: is_inflow,is_outflow
    integer :: ib,idir
    passed = .true.
    is_inflow( :,:) = inflow_type( :,:) > 0
    is_outflow(:,:) = outflow_type(:,:) > 0 .and. outflow_type(:,:) <= 2
    do idir = 1,3
      do ib = 0,1
        if( is_inflow( ib,idir) ) passed = passed .and. periods(idir) == 0
        if( is_outflow(ib,idir) ) passed = passed .and. periods(idir) == 0
      enddo
    enddo
    if(.not.passed) call write_error('Periodicity and in/outflow BCs cannot be combined.')
    do idir = 1,3
      do ib = 0,1
        if( is_inflow( ib,idir) .and. is_outflow(ib,idir ) ) passed = passed .and. .false.
      enddo
    enddo
    if(.not.passed) call write_error('BC cannot be both inflow and outflow.')
    call mpi_allreduce(MPI_IN_PLACE,passed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD)
  end subroutine chk_inoutflow
  !
  subroutine abortit
    implicit none
    if(myid == 0) write(stderr,*) ''
    if(myid == 0) write(stderr,*) '*** Simulation aborted due to errors in the input file ***'
    if(myid == 0) write(stderr,*) '    check dns.in'
    call MPI_FINALIZE()
    error stop
  end subroutine abortit
  subroutine write_error(message)
    character(len=*), intent(in) :: message
    if(myid == 0) write(stderr,*) 'ERROR: '//message
  end subroutine write_error
end module mod_sanity
