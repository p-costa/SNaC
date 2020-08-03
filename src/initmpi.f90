module mod_initmpi
  use mpi
  use mod_types
  use mod_common_mpi, only: ierr,myid,myid_block,comm_block
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(my_block,id_first,dims,cbc,bc,periods,lo,hi,ng,nb,is_bound,halos)
    implicit none
    integer         , intent(in   )                   :: my_block,id_first
    integer         , intent(inout), dimension(    3) :: dims
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc
    real(rp)        , intent(in   ), dimension(0:1,3) ::  bc
    integer         , intent(in   ), dimension(    3) :: periods
    integer         , intent(inout), dimension(    3) :: lo,hi
    integer         , intent(out  ), dimension(    3) :: ng
    integer         , intent(out  ), dimension(0:1,3) :: nb
    logical         , intent(out  ), dimension(0:1,3) :: is_bound
    integer         , intent(out  ), dimension(    3) :: halos
    integer                 :: nrank
    integer         , dimension(0:dims(1)-1,0:dims(2)-1,0:dims(3)-1) :: id_grid
    integer, dimension(  3) :: n,start,coords,lo_g,hi_g
    integer, allocatable, dimension(:,:  ) :: lo_all,hi_all,l,lo_g_all,hi_g_allo_g_all,hi_g_all
    integer, allocatable, dimension(  :  ) :: blocks_all
    real(rp)        , allocatable, dimension(:,:,:) ::  bc_all
    character(len=1), allocatable, dimension(:,:,:) :: cbc_all
    integer, dimension(3,3) :: eye
    integer                 :: i,j,k,idir,iidir,inb,irank,id_coords
    logical                 :: is_nb,found_friend
    integer                 :: ntot,ntot_max,ntot_min,ntot_sum
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
          if( myid == id_first + id_coords ) then
            coords(:) = [i,j,k]
          endif
          id_grid(i,j,k) = id_first + id_coords
          id_coords = id_coords+1
        enddo
      enddo
    enddo
    eye(:,:) = 0
    do idir=1,3
      eye(idir,idir) = 1
    enddo
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
    allocate(lo_all(3,0:nrank-1),hi_all(3,0:nrank-1),lo_g_all(3,0:nrank-1),hi_g_all(3,0:nrank-1), &
             blocks_all(0:nrank-1),cbc_all(0:1,3,0:nrank-1),bc_all(0:1,3,0:nrank-1))
    call MPI_ALLGATHER(lo      ,3,MPI_INTEGER,lo_all    ,3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(hi      ,3,MPI_INTEGER,hi_all    ,3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(lo_g    ,3,MPI_INTEGER,lo_g_all  ,3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(hi_g    ,3,MPI_INTEGER,hi_g_all  ,3,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(my_block,1,MPI_INTEGER,blocks_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER(cbc,6,MPI_CHARACTER,cbc_all,6,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
    call MPI_ALLGATHER( bc,6,MPI_REAL_RP  , bc_all,6,MPI_REAL_RP  ,MPI_COMM_WORLD,ierr)
    do idir=1,3
      do inb=0,1
        if(cbc(inb,idir) == 'F') then
          if(nint(bc(inb,idir)) < 1 .or. nint(bc(inb,idir)) > blocks_all(nrank-1))  then
            write(stderr,*) 'ERROR: invalid connectivity for block ', my_block, '.'
            write(stderr,*) 'Block ', nint(bc(inb,idir)), ' does not exist.'
            write(stderr,*) ''
            error stop
          endif
        endif
        if(nb(inb,idir) == MPI_PROC_NULL.and.cbc(inb,idir) == 'F') then
          is_nb = .false.
          do irank=0,nrank-1
            if(cbc(inb,idir) == cbc_all(1-inb,idir,irank).and. blocks_all(myid) == nint(bc_all(1-inb,idir,irank))) then
              is_nb = .true.
              do iidir = 1,3 ! ensure that the sub blocks blocks share the same extents in the two directions normal to the boundary
                if(iidir /= idir) then
                  is_nb = is_nb.and. &
                  lo(iidir) == lo_all(iidir,irank).and. &
                  hi(iidir) == hi_all(iidir,irank)
                endif
              enddo
              !
              ! these checks are not really needed because we specify the links
              ! in the block file, but they may be useful in future and to trap
              ! inconsistent in inputs
              !
              if(is_nb) then
                if(      inb == 0 .and. lo(idir) == lo_g(idir)) then
                  if(    lo(idir) == hi_all(idir,irank)+1) then
                    nb(inb,idir) = irank
                  elseif(lo(idir) <  hi_all(idir,irank)+1) then
                    if(periods(idir) == hi_all(idir,irank)-lo(idir)+1) then
                      nb(inb,idir) = irank
                    elseif(hi_all(idir,irank) == hi_g_all(idir,irank)) then
                      write(stderr,*) 'ERROR: Inconsistent periodic boundary condition (?).'
                      write(stderr,*) 'Expected: periods(',idir,') = ',hi_all(idir,irank)-lo(idir)+1
                      write(stderr,*) 'Found   : periods(',idir,') = ',periods(idir)
                      write(stderr,*) 'Error when connecting block: ' ,my_block
                      write(stderr,*) ''
                      error stop
                    endif
                  endif
                elseif ( inb == 1 .and. hi(idir) == hi_g(idir)) then
                  if(    hi(idir) == lo_all(idir,irank)-1) then
                    nb(inb,idir) = irank
                  elseif(hi(idir) >  lo_all(idir,irank)-1) then
                    if(periods(idir) == hi(idir)-lo_all(idir,irank)+1) then
                      nb(inb,idir) = irank
                    elseif(lo_all(idir,irank) == lo_g_all(idir,irank)) then
                      write(stderr,*) 'ERROR: Inconsistent periodic boundary condition (?).'
                      write(stderr,*) 'Expected: periods(',idir,') = ',hi(idir)-lo_all(idir,irank)+1
                      write(stderr,*) 'Found   : periods(',idir,') = ',periods(idir)
                      write(stderr,*) 'Error when connecting block: ' ,my_block
                      write(stderr,*) ''
                      error stop
                    endif
                  endif
                else
                  write(stderr,*) 'ERROR: Expected connectivity between blocks',my_block,' and ',blocks_all(irank),'not found.'
                  write(stderr,*) ''
                  error stop
                endif
              endif
            endif
          enddo
        endif
      enddo
    enddo
    deallocate(lo_all,hi_all,blocks_all)
    do idir=1,3
      do inb=0,1
        if(cbc(inb,idir) == 'F'.and.nb(inb,idir) == MPI_PROC_NULL) then
          write(stderr,*) 'ERROR: Expected connectivity between blocks',my_block,' and ',nint(bc(inb,idir)), ' is not possible.'
          write(stderr,*) 'E.g. Blocks must share the same boundaries.'
          write(stderr,*) ''
          error stop
        endif
      enddo
    enddo
    !
    do idir=1,3
      is_bound(:,idir) = .false.
      where(nb(0:1,idir) == MPI_PROC_NULL) is_bound(0:1,idir) = .true.
      call makehalo(idir,1,n,halos(idir))
    enddo
    !
    ! check distribution of grid points over the different tasks
    !
    ntot = product(hi(:)-lo(:)+1)
    call MPI_ALLREDUCE(ntot,ntot_min,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(ntot,ntot_max,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
    call MPI_ALLREDUCE(ntot,ntot_sum,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
    write(stdout,*) 'Maximum, minimum, average, and normalized average average number of grid points for task ',myid, &
                    '(block',my_block,'): ',ntot_max, &
                                            ntot_min, &
                                            (1.*ntot_sum)/(1._rp*nrank), &
                                            (1._rp*ntot*nrank)/(1._rp*ntot_sum)
    return
  end subroutine initmpi
  subroutine makehalo(idir,nh,n,halo)
    implicit none
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
