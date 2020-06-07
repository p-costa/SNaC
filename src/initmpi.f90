module mod_initmpi
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr,myid,comm_cart
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
    integer, intent(out), dimension(3    ) :: halos
    integer :: nrank
    logical, dimension(3) :: periods
    integer, dimension(3) :: n,coords
    integer :: idir
    !
    ! sanity check
    !
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    if(nrank.ne.product(dims(:))) then
      dims(:) = 0
      CALL MPI_DIMS_CREATE(nrank,3,dims,ierr)
      if(myid.eq.0) then
        write(stdout,*) 'WARNING: product(dims(:)) not equal to the number of MPI tasks!'
        write(stdout,*) 'dims changed with MPI_DIMS_CREATE to ', dims(1),dims(2),dims(3)
      endif
    endif
    !
    ! generate Cartesian topology
    !
    periods(:) = .false.
    where (bc(0,:)//bc(1,:).eq.'PP') periods(:) = .true.
    call MPI_CART_CREATE(MPI_COMM_WORLD,3,dims,periods,.true.,comm_cart,ierr)
    call MPI_CART_COORDS(comm_cart,myid,3,coords,ierr)
    !
    ! determine array extents for possibly uneven data
    !
    n(:) = ng(:)/dims(:)
    where(coords(:)+1.le.mod(ng(:),dims(:))) n(:) = n(:) + 1
    lo(:) = 1    + coords(:)*n(:)
    hi(:) = n(:) + coords(:)*n(:)
    where(coords(:)+1.gt.mod(ng(:),dims(:)))
      lo(:) = lo(:) + mod(ng(:),dims(:))
      hi(:) = hi(:) + mod(ng(:),dims(:))
    end where
    !
    ! generate neighbor arrays and derived types for halo regions
    !
    do idir=1,3
      call MPI_CART_SHIFT(comm_cart,idir-1,1,nb(0,idir),nb(1,idir),ierr)
      is_bound(:,idir) = .false.
      where(nb(0:1,idir).eq.MPI_PROC_NULL) is_bound(0:1,idir) = .true.
      call makehalo(idir,1,n,halos(idir))
    enddo
    return
  end subroutine initmpi
  subroutine makehalo(idir,nh,n,halo)
    integer, intent(in ) :: idir,nh
    integer, intent(in ), dimension(3) :: n
    integer, intent(out) :: halo
    integer, dimension(3) :: nn
    nn(:) = n(:) + 2*nh
    select case(idir)
    case(1)
      call MPI_TYPE_VECTOR(nn(2)*nn(3),nh            ,nn(1)            ,MPI_REAL_RP,halo,ierr)
    case(2)
      call MPI_TYPE_VECTOR(      nn(3),nh*nn(1)      ,nn(1)*nn(2)      ,MPI_REAL_RP,halo,ierr)
    case(3)
      call MPI_TYPE_VECTOR(          1,nh*nn(1)*nn(2),nn(1)*nn(2)*nn(3),MPI_REAL_RP,halo,ierr)
    end select
    call MPI_TYPE_COMMIT(halo,ierr)
    return
  end subroutine makehalo
end module mod_initmpi
