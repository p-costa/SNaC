module mod_initmpi
  use mpi_f08
  use mod_types
  use mod_common_mpi, only: myid,comm_cart
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(ng,dims,bc,lo,hi,nb,is_bound,halos)
    implicit none
    integer, intent(in    ), dimension(3) :: ng
    integer, intent(inout ), dimension(3) :: dims
    character(len=1), intent(in), dimension(0:1,3) :: bc
    integer, intent(out), dimension(3) :: lo,hi
    integer, intent(out), dimension(0:1,3) :: nb
    logical, intent(out), dimension(0:1,3) :: is_bound
    type(MPI_DATATYPE), intent(out), dimension(3    ) :: halos
    integer :: nrank
    logical, dimension(3) :: periods
    integer, dimension(3) :: n,coords
    integer :: idir
    !
    ! sanity check
    !
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank)
    if(nrank /= product(dims(:))) then
      dims(:) = 0
      call MPI_DIMS_CREATE(nrank,3,dims)
      if(myid == 0) then
        write(stdout,*) 'WARNING: product(dims(:)) not equal to the number of MPI tasks!'
        write(stdout,*) 'dims changed with MPI_DIMS_CREATE to ', dims(1),dims(2),dims(3)
      endif
    endif
    !
    ! generate Cartesian topology
    !
    periods(:) = .false.
    where (bc(0,:)//bc(1,:) == 'PP') periods(:) = .true.
    call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm_cart)
    call MPI_CART_COORDS(comm_cart,myid,3,coords)
    !
    ! determine array extents for possibly uneven data
    !
    n(:) = ng(:)/dims(:)
    where(coords(:)+1 <= mod(ng(:),dims(:))) n(:) = n(:) + 1
    lo(:) = 1    + coords(:)*n(:)
    hi(:) = n(:) + coords(:)*n(:)
    where(coords(:)+1 >  mod(ng(:),dims(:)))
      lo(:) = lo(:) +    mod(ng(:),dims(:))
      hi(:) = hi(:) +    mod(ng(:),dims(:))
    end where
    !
    ! generate neighbor arrays and derived types for halo regions
    !
    do idir=1,3
      call MPI_CART_SHIFT(comm_cart,idir-1,1,nb(0,idir),nb(1,idir))
      is_bound(:,idir) = .false.
      where(nb(0:1,idir) == MPI_PROC_NULL) is_bound(0:1,idir) = .true.
      call makehalo(idir,1,n,halos(idir))
    enddo
  end subroutine initmpi
  subroutine makehalo(idir,nh,n,halo)
    implicit none
    integer, intent(in ) :: idir,nh
    integer, intent(in ), dimension(3) :: n
    type(MPI_DATATYPE), intent(out) :: halo
    integer, dimension(3) :: nn
    nn(:) = n(:) + 2*nh
    select case(idir)
    case(1)
      call MPI_TYPE_VECTOR(nn(2)*nn(3),nh            ,nn(1)            ,MPI_REAL_RP,halo)
    case(2)
      call MPI_TYPE_VECTOR(      nn(3),nh*nn(1)      ,nn(1)*nn(2)      ,MPI_REAL_RP,halo)
    case(3)
      call MPI_TYPE_VECTOR(          1,nh*nn(1)*nn(2),nn(1)*nn(2)*nn(3),MPI_REAL_RP,halo)
    end select
    call MPI_TYPE_COMMIT(halo)
  end subroutine makehalo
end module mod_initmpi
