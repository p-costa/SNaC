module mod_initmpi
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr,myid,myid_block,comm_block
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(my_block,dims,cbc,bc,lo,hi,ng,periods,nb,is_bound,halos)
    implicit none
    integer         , intent(in   )                   :: my_block
    integer         , intent(inout), dimension(    3) :: dims
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc
    real(rp)        , intent(in   ), dimension(0:1,3) ::  bc
    integer         , intent(inout), dimension(    3) :: lo,hi
    integer         , intent(out  ), dimension(    3) :: ng,periods
    integer         , intent(out  ), dimension(0:1,3) :: nb
    logical         , intent(out  ), dimension(0:1,3) :: is_bound
    integer         , intent(out  ), dimension(    3) :: halos
    integer                 :: nrank
    integer         , dimension(0:dims(1)-1,0:dims(2)-1,0:dims(3)-1) :: id_grid
    integer, dimension(  3) :: n,start,coords,lo_g,hi_g
    integer, allocatable, dimension(:,:) :: lo_all,hi_all
    integer, allocatable, dimension(  :) :: blocks_all
    integer, dimension(3,3) :: eye
    integer                 :: i,j,k,idir,iidir,inb,irank,id_coords
    logical                 :: is_nb
    !
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank,ierr)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,my_block,myid,comm_block,ierr)
    call MPI_COMM_RANK(comm_block, myid_block, ierr)
    !
    ! generate Cartesian topology
    !
    hi_g(:) = hi(:)
    lo_g(:) = lo(:)
    ng(:) = hi_g(:)-lo_g(:)+1
    !
    ! determine array extents for possibly uneven data
    !
    coords(:) = 0
    id_coords = 0
    do k=0,dims(3)-1
      do j=0,dims(2)-1
        do i=0,dims(1)-1
          id_coords = id_coords+1
          if( myid == id_coords ) then
            coords(1) = i
            coords(2) = j
            coords(3) = k
          endif
          id_grid(i,j,k) = id_coords
        enddo
      enddo
    enddo
    forall (idir=1:3) eye(idir,idir) = 1
    nb(:,:) = MPI_PROC_NULL
    do idir=1,3
      if(coords(idir)-1 >= 0           ) &
        nb(0,idir) = id_grid(coords(1)-eye(1,idir),coords(2)-eye(2,idir),coords(3)-eye(3,idir))
      if(coords(idir)+1 <= dims(idir)-1) &
        nb(1,idir) = id_grid(coords(1)+eye(1,idir),coords(2)+eye(2,idir),coords(3)+eye(3,idir))
    enddo
    n(:) = ng(:)/dims(:)
    where(coords(:)+1 <= mod(ng(:),dims(:))) n(:) = n(:) + 1
    lo(:) = lo_g(:)   + (coords(:)  )*n(:)
    hi(:) = lo_g(:)-1 + (coords(:)+1)*n(:)
    where(coords(:)+1 >  mod(ng(:),dims(:)))
      lo(:) = lo(:) +    mod(ng(:),dims(:))
      hi(:) = hi(:) +    mod(ng(:),dims(:))
    end where
    !
    call MPI_ALLGATHER(lo      ,3,MPI_INTEGER,lo_all    ,3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(hi      ,3,MPI_INTEGER,hi_all    ,3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(my_block,1,MPI_INTEGER,blocks_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    allocate(lo_all(0:nrank-1,3),hi_all(0:nrank-1,3),blocks_all(0:nrank-1))
    do idir=1,3
      do inb=0,1
        if(nb(inb,idir) == MPI_PROC_NULL.and.cbc(inb,idir) == 'F') then
          do irank=0,nrank-1
            if(nint(bc(inb,idir)) == blocks_all(irank)) then
              is_nb = .true. 
              do iidir = 1,3 ! ensure that the sub blocks blocks share the same extents in the two directions normal to the boundary
                if(iidir /= idir) then
                  is_nb = is_nb.and. &
                  lo(iidir).eq.lo_all(irank,iidir).and. &
                  hi(iidir).eq.hi_all(irank,iidir)
                endif
              enddo
              !
              ! these checks are not really needed because we specify the links
              ! in the block file, but they may be useful in future and to trap
              ! inconsistent in inputs
              !
              if(is_nb) then
                if(      inb == 0 ) then
                  if(    lo(idir) == hi_all(irank,idir)+1) then
                    nb(inb,idir) = irank
                  elseif(lo(idir) <  hi_all(irank,idir)+1) then
                    if(periods(idir) == hi_all(irank,idir)-lo(idir)+1) then
                      nb(inb,idir) = irank
                    else
                      write(stderr,*) 'ERROR: Inconsistent periodic boundary condition.'
                    endif
                  endif
                elseif ( inb == 1 ) then
                  if(    hi(idir) == lo_all(irank,idir)-1) then
                    nb(inb,idir) = irank
                  elseif(hi(idir) >  lo_all(irank,idir)-1) then
                    if(periods(idir) == hi(idir)-lo_all(irank,idir)+1) then
                      nb(inb,idir) = irank
                    else
                      write(stderr,*) 'ERROR: Inconsistent periodic boundary condition.'
                    endif
                  endif
                endif
                if(nb(inb,idir) == MPI_PROC_NULL) then
                  write(stderr,*) 'ERROR: Expected connectivity between blocks not found.'
                  error stop
                endif
              endif
            endif
          enddo
        endif
      enddo
    enddo
    deallocate(lo_all,hi_all)
    !
    do idir=1,3
      is_bound(:,idir) = .false.
      where(nb(0:1,idir) == MPI_PROC_NULL) is_bound(0:1,idir) = .true.
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
