module mod_load
  use mpi
  use mod_types
  implicit none
  private
  public load
  contains
  subroutine load(io,myid,filename,ng,nghost,lo,hi,u,v,w,p,time,istep)
    !
    ! reads/writes a restart file
    !
    implicit none
    character(len=1), intent(in) :: io
    integer         , intent(in) :: myid
    character(len=*), intent(in) :: filename
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in)               :: nghost
    real(rp), intent(inout), dimension(lo(1)-nghost:,lo(2)-nghost:,lo(3)-nghost:) :: u,v,w,p
    real(rp), intent(inout) :: time,istep
    real(rp), dimension(2) :: fldinfo
    integer , dimension(3) :: n
    integer :: fh,ierr,lenr,nreals_myid
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp,good
    integer, dimension(3) :: sizes,subsizes,starts
    integer :: type_glob,type_loc
    !
    n(:) = hi(:)-lo(:)+1
    lenr = sizeof(time)
    !
    sizes(:)    = ng(:)
    subsizes(:) = n(:)
    starts(:)   = lo(:) - 1 ! starts from 0
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_glob,ierr)
    sizes(:)    = n(:) + nghost
    subsizes(:) = n(:)
    starts(:)   = 0 + nghost
    call MPI_TYPE_CREATE_SUBARRAY(3,sizes,subsizes,starts,MPI_ORDER_FORTRAN,MPI_REAL_RP,type_glob,ierr)
    call MPI_TYPE_COMMIT(type_loc,ierr)
    select case(io)
    case('r')
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
           MPI_MODE_RDONLY, MPI_INFO_NULL,fh, ierr)
      !
      ! check file size first
      !
      call MPI_FILE_GET_SIZE(fh,filesize,ierr)
      good = (product(ng)*4+2)*lenr
      if(filesize.ne.good) then
        if(myid.eq.0) write(error_unit,*) ''
        if(myid.eq.0) write(error_unit,*) '*** Simulation aborted due a checkpoint file with incorrect size ***'
        if(myid.eq.0) write(error_unit,*) '    file: ', filename, ' | expected size: ', good, '| actual size: ', filesize
        call MPI_FINALIZE(ierr)
        error stop
      endif
      !
      ! read
      !
      disp = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_READ_ALL(fh,u,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_READ_ALL(fh,v,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_READ_ALL(fh,w,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_READ_ALL(fh,p,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      !
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      nreals_myid = 0
      if(myid.eq.0) nreals_myid = 2
      time  = fldinfo(1)
      istep = fldinfo(2)
      call MPI_FILE_READ(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
      call MPI_BROADCAST(fldinfo,2,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
    case('w')
      !
      ! write
      !
      call MPI_FILE_OPEN(MPI_COMM_WORLD, filename                 , &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh, ierr)
      filesize = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_SIZE(fh,filesize,ierr)
      disp = 0_MPI_OFFSET_KIND
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_WRITE_ALL(fh,u,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_WRITE_ALL(fh,v,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_WRITE_ALL(fh,w,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,type_glob,'native',MPI_INFO_NULL,ierr)
      call MPI_FILE_WRITE_ALL(fh,p,1,type_loc,MPI_STATUS_IGNORE,ierr)
      disp = disp+product(ng)*lenr
      !
      call MPI_FILE_SET_VIEW(fh,disp,MPI_REAL_RP,MPI_REAL_RP,'native',MPI_INFO_NULL,ierr)
      fldinfo = [time,istep]
      nreals_myid = 0
      if(myid.eq.0) nreals_myid = 2
      time  = fldinfo(1)
      istep = fldinfo(2)
      call MPI_FILE_WRITE(fh,fldinfo,nreals_myid,MPI_REAL_RP,MPI_STATUS_IGNORE,ierr)
      call MPI_FILE_CLOSE(fh,ierr)
    end select
    return
  end subroutine load
end module mod_load
