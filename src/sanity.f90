module mod_sanity
  use mpi_f08
  use mod_common_mpi, only: myid
  use mod_types
  implicit none
  private
  public test_sanity
  contains
  subroutine test_sanity(ng,gr,stop_type,cbcvel,cbcpre,is_forced)
    !
    ! performs some a priori checks of the input files before the calculation starts
    !
    implicit none
    integer         , intent(in), dimension(3      ) :: ng
    real(rp)        , intent(in), dimension(3      ) :: gr
    logical         , intent(in), dimension(3      ) :: stop_type
    character(len=1), intent(in), dimension(0:1,3,3) :: cbcvel
    character(len=1), intent(in), dimension(0:1,3  ) :: cbcpre
    logical         , intent(in), dimension(3      ) :: is_forced
    logical :: passed
    !
    call chk_grid(gr,passed);                    if(.not.passed) call abortit
    call chk_stop_type(stop_type,passed);        if(.not.passed) call abortit
    call chk_bc(cbcvel,cbcpre,passed);           if(.not.passed) call abortit
    call chk_forcing(cbcpre,is_forced,passed);   if(.not.passed) call abortit 
  end subroutine test_sanity
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
  !
  subroutine chk_bc(cbcvel,cbcpre,passed)
    implicit none
    character(len=2), parameter, dimension(5) :: bcs = ['ND', 'DN', 'NN', 'DD', 'PP']
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
        passed_loc = passed_loc.and.any(bc01v == bcs(:))
      enddo
    enddo
    if(.not.passed_loc) call write_error('velocity BCs not valid.')
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.any(bc01v == bcs(:))
    enddo
    if(.not.passed_loc) call write_error('pressure BCs not valid.')
    passed = passed.and.passed_loc
    !
    passed_loc = .true.
    do idir=1,3
      ivel = idir
      bc01v = cbcvel(0,idir,ivel)//cbcvel(1,idir,ivel)
      bc01p = cbcpre(0,idir)//cbcpre(1,idir)
      passed_loc = passed_loc.and.( (bc01v == 'PP'.and.bc01p == 'PP').or. &
                                    (bc01v == 'ND'.and.bc01p == 'DN').or. &
                                    (bc01v == 'DN'.and.bc01p == 'ND').or. &
                                    (bc01v == 'DD'.and.bc01p == 'NN').or. &
                                    (bc01v == 'NN'.and.bc01p == 'DD') )
    enddo
    if(.not.passed_loc) call write_error('velocity and pressure BCs not compatible.')
    passed = passed.and.passed_loc
    !
  end subroutine chk_bc
  !
  subroutine chk_forcing(cbcpre,is_forced,passed)
    implicit none
    character(len=1), intent(in ), dimension(0:1,3) :: cbcpre
    logical         , intent(in ), dimension(3    ) :: is_forced
    logical         , intent(out) :: passed
    integer :: idir
    passed = .true.
    !
    ! 1) check for compatibility between pressure BCs and forcing BC
    !
    do idir=1,3
      if(is_forced(idir)) then
        passed = passed.and.(cbcpre(0,idir)//cbcpre(1,idir) == 'PP')
      endif
    enddo
    if(.not.passed) &
      call write_error('Flow cannot be forced in a non-periodic direction; check the BCs and is_forced in dns.in.')
  end subroutine chk_forcing
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
