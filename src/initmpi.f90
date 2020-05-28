module mod_initmpi
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(ng,dims,lo,hi)
    implicit none
    integer, intent(in), dimension(3) :: ng
    character(len=1), intent(in), dimension(0:1,3) :: bc
    logical, dimension(3) :: periods
    integer, dimension(3) :: n
    integer :: idir
    !
    periods(:) = .false.
    forall (bc(0,:)//bc(1,:).eq.'PP') periods(:) = .true.
    n(:) = ng(:)/dims(:)
    call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,comm_cart,ierr)
    call MPI_CART_COORDS(comm_cart,myid,3,coords)
    lo(:) = 1    + coord(:)*n(:)
    hi(:) = n(:) + coord(:)*n(:)
    !
    call MPI_CART_SHIFT(comm_cart,0,1,left,right,ierr)
    call MPI_CART_SHIFT(comm_cart,1,1,front,back,ierr)
    !
    do idir=1,3
      call makehalo(idir,1,n,halo(idir))
    enddo
    return
  end subroutine initmpi
  subroutine makehalo(idir,nghost,n,halo)
    integer, intent(in ) :: idir,nghost
    integer, intent(in ), dimension(3) :: n
    integer, intent(out) :: halo
    integer, dimension(3) :: nn
    nn(:) = nn(:) + 2*nghost
    select case(idir)
    case(1)
      call MPI_TYPE_VECTOR(nn(2)*nn(3),nghost            ,nn(1)            ,MPI_REAL_RP,halo,ierr)
    case(2)
      call MPI_TYPE_VECTOR(      nn(3),nghost*nn(1)      ,nn(1)*nn(2)      ,MPI_REAL_RP,halo,ierr)
    case(3)
      call MPI_TYPE_VECTOR(          1,nghost*nn(1)*nn(2),nn(1)*nn(2)*nn(3),MPI_REAL_RP,halo,ierr)
    end select
    return
  end subroutine makehalo
end module mod_initmpi
