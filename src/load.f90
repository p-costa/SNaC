module mod_load
  use mpi
  use mod_common_mpi, only:myid,ierr
  use mod_types
  implicit none
  private
  public load,io_field
  contains
  subroutine load(io,filename,ng,nghost,lo,hi,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    character(len=*), intent(in) :: filename
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in)               :: nghost
    real(rp), intent(inout), dimension(lo(1)-nghost:,lo(2)-nghost:,lo(3)-nghost:) :: u,v,w,p
    real(rp), intent(inout) :: time
    integer , intent(inout) :: istep
    real(rp), dimension(2) :: fldinfo
    integer :: fh,nreals_myid
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    !
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      good = (product(ng)*4+2)*sizeof(1._rp)
      if(filesize.ne.good) then
        if(myid.eq.0) write(stderr,*) ''
        if(myid.eq.0) write(stderr,*) '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) write(stderr,*) '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call MPI_FINALIZE(ierr)
        error stop
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call io_field('r',fh,ng,lo,hi,nghost,disp,u)
      call io_field('r',fh,ng,lo,hi,nghost,disp,v)
      call io_field('r',fh,ng,lo,hi,nghost,disp,w)
      call io_field('r',fh,ng,lo,hi,nghost,disp,p)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      nreals_myid = 0
      if(myid.eq.0) nreals_myid = 2
      call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_BCAST(fldinfo,2,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
      time  =      fldinfo(1)
      istep = nint(fldinfo(2))
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)
      disp = 0_MPI_OFFSET_KIND
      call io_field('w',fh,ng,lo,hi,nghost,disp,u)
      call io_field('w',fh,ng,lo,hi,nghost,disp,v)
      call io_field('w',fh,ng,lo,hi,nghost,disp,w)
      call io_field('w',fh,ng,lo,hi,nghost,disp,p)
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      fldinfo = [time,1._rp*istep]
      nreals_myid = 0
      if(myid.eq.0) nreals_myid = 2
      call MPI_FILE_WRITE(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
    end select
    return
  end subroutine load
  subroutine io_field(io,fh,ng,lo,hi,nghost,disp,var)
    implicit none
    character(len=1), intent(in)                 :: io
    integer , intent(in)                         :: fh
    integer , intent(in), dimension(3)           :: ng,lo,hi
    integer , intent(in)                         :: nghost
    integer(kind=MPI_OFFSET_KIND), intent(inout) :: disp
    real(rp), intent(in), dimension(lo(1)-nghost:,lo(2)-nghost:,lo(3)-nghost:) :: var
    integer :: ierr
    integer , dimension(3) :: n
    integer , dimension(3) :: sizes,subsizes,starts
    integer :: type_glob,type_loc
    n(:)        = hi(:)-lo(:)+1
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_glob,ierr)
    sizes(:)    = n(:) + nghost
    subsizes(:) = n(:)
    starts(:)   = 0 + nghost
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_loc ,ierr)
    call MPI_TYPE_COMMIT(type_loc,ierr)
      select case(io)
      case('r')
        call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
        call MPI_FILE_READ_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE,ierr)
      case('w')
        call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
        call MPI_FILE_WRITE_ALL(fh,var,1,type_loc,MPI_STATUS_IGNORE,ierr)
      end select
        disp = disp+product(ng)*sizeof(1._rp)
      return
  end subroutine io_field
end module mod_load
