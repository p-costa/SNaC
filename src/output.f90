module mod_output
  use mpi_f08
  use mod_common_mpi, only: myid,myid_block
  use mod_load      , only: io_field
  use mod_types
  implicit none
  private
  public out0d,out1d,out2d,out3d,write_log_output,write_visu_3d
  contains
  subroutine out0d(fname,n,var)
    !
    ! appends the first n entries of an array
    ! var to a file
    ! fname -> name of the file
    ! n     -> number of entries
    ! var   -> input array of real values
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in) :: n
    real(rp), intent(in), dimension(:) :: var
    integer :: iunit
    character(len=30) :: cfmt
    integer :: i
    !
    write(cfmt,'(A,I3,A)') '(',n,'E15.7)'
    if (myid == 0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) (var(i),i=1,n)
      close(iunit)
    end if
  end subroutine out0d
  !
  subroutine out1d(fname,lo,hi,lo_g,hi_g,idir,l,dx,dy,dz,x_g,y_g,z_g,comm,myrank,p)
    !
    ! writes the profile of a variable averaged over two domain directions
    !
    ! fname       -> name of the file
    ! lo  ,hi     -> local lower and upper bounds of input array
    !                containing all the points of the block
    ! idir        -> direction of the profile
    ! l           -> domain dimensions
    ! dx,dy,dz    -> grid spacings
    ! x_g,y_g,z_g -> coodinates of grid points for the entire block
    ! mpi_comm    -> communicator pertaining to the group of tasks of each block
    ! p           -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: lo,hi,lo_g,hi_g
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(lo(1)-1:  ) :: dx
    real(rp), intent(in), dimension(lo(2)-1:  ) :: dy
    real(rp), intent(in), dimension(lo(3)-1:  ) :: dz
    real(rp), intent(in), dimension(lo_g(1)-1:) :: x_g
    real(rp), intent(in), dimension(lo_g(2)-1:) :: y_g
    real(rp), intent(in), dimension(lo_g(3)-1:) :: z_g
    type(MPI_COMM), intent(in)                  :: comm
    integer , intent(in)                        :: myrank
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), allocatable, dimension(:) :: p1d
    integer, dimension(3) :: ng
    integer :: i,j,k
    integer :: iunit
    !
    ng(:) = hi_g(:)-lo_g(:) + 1
    select case(idir)
    case(3)
      allocate(p1d(lo_g(3):hi_g(3)))
      p1d(:) = 0._rp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            p1d(k) = p1d(k) + p(i,j,k)*dx(i)*dy(j)/(l(1)*l(2))
          end do
        end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,p1d,ng(3),MPI_REAL_RP,MPI_SUM,comm)
      if(myrank == 0) then
        open(newunit=iunit,file=fname)
        do k=lo_g(3),hi_g(3)
          write(iunit,'(2E15.7)') z_g(k),p1d(k)
        end do
        close(iunit)
      end if
    case(2)
      allocate(p1d(lo_g(2):hi_g(2)))
      p1d(:) = 0._rp
      do j=lo(2),hi(2)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            p1d(j) = p1d(j) + p(i,j,k)*dx(i)*dz(k)/(l(1)*l(3))
          end do
        end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,p1d,ng(2),MPI_REAL_RP,MPI_SUM,comm)
      if(myrank == 0) then
        open(newunit=iunit,file=fname)
        do j=lo_g(2),hi_g(2)
          write(iunit,'(2E15.7)') y_g(j),p1d(j)
        end do
        close(iunit)
      end if
    case(1)
      allocate(p1d(lo_g(1):hi_g(1)))
      p1d(:) = 0._rp
      do i=lo(1),hi(1)
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            p1d(i) = p1d(i) + p(i,j,k)*dy(j)*dz(k)/(l(2)*l(3))
          end do
        end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,p1d,ng(1),MPI_REAL_RP,MPI_SUM,comm)
      if(myrank == 0) then
        open(newunit=iunit,file=fname)
        do i=lo_g(1),hi_g(1)
          write(iunit,'(2E15.7)') x_g(i),p1d(i)
        end do
        close(iunit)
      end if
    end select
    deallocate(p1d)
  end subroutine out1d
  !
  subroutine out2d(fname,lo,hi,lo_g,hi_g,idir,l,dx,dy,dz,x_g,y_g,z_g,comm,myrank,p)
    !
    ! writes the profile of a variable averaged over two domain directions
    !
    ! fname       -> name of the file
    ! lo  ,hi     -> local lower and upper bounds of input array
    !                containing all the points of the block
    ! idir        -> direction along which the field is averaged
    ! l           -> domain dimensions
    ! dx,dy,dz    -> grid spacings
    ! x_g,y_g,z_g -> coodinates of grid points for the entire block
    ! mpi_comm    -> communicator pertaining to the group of tasks of each block
    ! p           -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in), dimension(3) :: lo,hi,lo_g,hi_g
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(lo(1)-1:  ) :: dx
    real(rp), intent(in), dimension(lo(2)-1:  ) :: dy
    real(rp), intent(in), dimension(lo(3)-1:  ) :: dz
    real(rp), intent(in), dimension(lo_g(1)-1:) :: x_g
    real(rp), intent(in), dimension(lo_g(2)-1:) :: y_g
    real(rp), intent(in), dimension(lo_g(3)-1:) :: z_g
    type(MPI_COMM), intent(in)                  :: comm
    integer , intent(in)                        :: myrank
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), allocatable, dimension(:,:) :: p2d
    integer, dimension(3) :: ng
    integer :: i,j,k
    integer :: iunit
    !
    ng(:) = hi_g(:)-lo_g(:) + 1
    select case(idir)
    case(1)
      allocate(p2d(lo_g(2):hi_g(2),lo_g(3):hi_g(3)))
      p2d(:,:) = 0._rp
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          p2d(:,:) = 0._rp
          do i=lo(idir),hi(idir)
            p2d(j,k) = p2d(j,k) + p(i,j,k)*dx(i)/l(idir)
          end do
        end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,p2d,ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,comm)
      if(myrank == 0) then
        open(newunit=iunit,file=fname)
        do k=lo_g(3),hi_g(3)
          do j=lo_g(2),hi_g(2)
            write(iunit,'(3E15.7)') y_g(j),z_g(k),p2d(j,k)
          end do
        end do
        close(iunit)
      end if
    case(2)
      allocate(p2d(lo_g(1):hi_g(1),lo_g(3):hi_g(3)))
      p2d(:,:) = 0._rp
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          do j=lo(idir),hi(idir)
            p2d(i,k) = p2d(i,k) + p(i,j,k)*dy(j)/l(idir)
          end do
        end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,p2d,ng(1)*ng(3),MPI_REAL_RP,MPI_SUM,comm)
      if(myrank == 0) then
        open(newunit=iunit,file=fname)
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            write(iunit,'(3E15.7)') x_g(i),z_g(k),p2d(i,k)
          end do
        end do
        close(iunit)
      end if
    case(3)
      allocate(p2d(lo_g(1):hi_g(1),lo_g(2):hi_g(2)))
      p2d(:,:) = 0._rp
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          do k=lo(idir),hi(idir)
            p2d(i,j) = p2d(i,j) + p(i,j,k)*dz(k)/l(idir)
          end do
        end do
      end do
      call mpi_allreduce(MPI_IN_PLACE,p2d,ng(1)*ng(2),MPI_REAL_RP,MPI_SUM,comm)
      if(myrank == 0) then
        open(newunit=iunit,file=fname)
        do j=lo_g(2),hi_g(2)
          do i=lo_g(1),hi_g(1)
            write(iunit,'(3E15.7)') x_g(i),y_g(j),p2d(i,j)
          end do
        end do
        close(iunit)
      end if
    end select
    deallocate(p2d)
  end subroutine out2d
  !
  subroutine out3d(fname,comm,lo,hi,ng,nskip,p)
    !
    ! saves a 3D scalar field into a binary file
    !
    ! fname  -> name of the output file
    ! lo,hi  -> local lower and upper bounds of input array
    !           in global coordinates
    ! ng     -> global sizes of the input array
    ! nskip  -> array with the step size for which the
    !           field is written; i.e.: (/1,1,1/)
    !           writes the full field
    !           n.b.: not implemented for now; it will
    !                 always write the full array
    ! p      -> 3D input scalar field
    !
    implicit none
    character(len=*), intent(in) :: fname
    integer , intent(in   ), dimension(3) :: lo,hi,ng,nskip
    type(MPI_COMM) , intent(in   )               :: comm
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    type(MPI_FILE) :: fh
    integer(kind=MPI_OFFSET_KIND) :: filesize,disp
    !
    call MPI_FILE_OPEN(comm, fname, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL,fh)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize)
    disp = 0_MPI_OFFSET_KIND
    call io_field('w',fh,ng,[1,1,1],lo,hi,disp,p)
    call MPI_FILE_CLOSE(fh)
  end subroutine out3d
  !
  subroutine write_log_output(fname,fname_fld,varname,nmin,nmax,nskip,time,istep)
    !
    ! appends information about a saved binary file to a log file
    ! this file is used to generate a xdmf file for visualization of field data
    !
    ! fname     -> name of the output log file
    ! fname_fld -> name of the saved binary file (excluding the directory)
    ! varname   -> name of the variable that is saved
    ! nmin      -> first element of the field that is saved in each direction, e.g. (/1,1,1/)
    ! nmax      -> last  element of the field that is saved in each direction, e.g. (/ng(1),ng(2),ng(3)/)
    ! nskip     -> step size between nmin and nmax, e.g. (/1,1,1/) if the whole array is saved
    ! time      -> physical time
    ! istep     -> time step number
    !
    implicit none
    character(len=*), intent(in)       :: fname,fname_fld,varname
    integer , intent(in), dimension(3) :: nmin,nmax,nskip
    real(rp), intent(in)               :: time
    integer , intent(in)               :: istep
    character(len=100) :: cfmt
    integer :: iunit
    !
    write(cfmt, '(A)') '(A,A,A,9i5,E15.7,i10)'
    if (myid_block == 0) then
      open(newunit=iunit,file=fname,position='append')
      write(iunit,trim(cfmt)) trim(fname_fld),' ',trim(varname),nmin,nmax,nskip,time,istep
      close(iunit)
    end if
  end subroutine write_log_output
  !
  subroutine write_visu_3d(datadir,fname_bin,comm,fname_log,varname,lo,hi,ng,nmin,nmax,nskip,time,istep,p)
    !
    ! wraps the calls of out3d and write_log_output into the same subroutine
    !
    implicit none
    character(len=*), intent(in)             :: datadir,fname_bin,fname_log,varname
    type(MPI_COMM) , intent(in   )           :: comm
    integer , intent(in   ), dimension(3)    :: lo,hi,ng,nmin,nmax,nskip
    real(rp), intent(in   )                  :: time
    integer , intent(in   )                  :: istep
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    !
    call out3d(trim(datadir)//trim(fname_bin),comm,lo,hi,ng,nskip,p)
    call write_log_output(trim(datadir)//trim(fname_log),trim(fname_bin),trim(varname),nmin,nmax,nskip,time,istep)
  end subroutine write_visu_3d
end module mod_output
