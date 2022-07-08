module mod_initmpi
  use mpi_f08
  use mod_types
  use mod_common_mpi, only: myid,myid_block,comm_block
  implicit none
  private
  public initmpi
  contains
  subroutine initmpi(my_block,nblocks,id_first,dims,cbc,bc,periods,lmin,lmax,gt,gr,lo,hi,lo_g,hi_g,ng,nb,is_bound,halos)
    implicit none
    integer         , intent(in   )                   :: my_block,nblocks,id_first
    integer         , intent(inout), dimension(    3) :: dims
    character(len=1), intent(in   ), dimension(0:1,3) :: cbc
    real(rp)        , intent(in   ), dimension(0:1,3) ::  bc
    integer         , intent(in   ), dimension(    3) :: periods
    real(rp)        , intent(in   ), dimension(    3) :: lmin,lmax
    integer         , intent(in   ), dimension(    3) :: gt
    real(rp)        , intent(in   ), dimension(    3) :: gr
    integer         , intent(inout), dimension(    3) :: lo,hi,lo_g,hi_g
    integer         , intent(out  ), dimension(    3) :: ng
    integer         , intent(out  ), dimension(0:1,3) :: nb
    logical         , intent(out  ), dimension(0:1,3) :: is_bound
    type(MPI_DATATYPE), intent(out  ), dimension(    3) :: halos
    integer                 :: nrank
    integer , dimension(  3) :: n,start,coords
    integer , allocatable, dimension(:,:) :: ng_all,lo_all,hi_all,gt_all,l,lo_g_all,hi_g_allo_g_all,hi_g_all
    real(rp), allocatable, dimension(:,:) :: lmin_all,lmax_all,gr_all
    integer , allocatable, dimension(:) :: blocks_all
    logical , allocatable, dimension(:) :: is_done,is_unlocked
    real(rp)        , allocatable, dimension(:,:,:) ::  bc_all
    character(len=1), allocatable, dimension(:,:,:) :: cbc_all
    integer                 :: i,j,k,idir,iidir,inb,irank,iblock,ifriend,ioffset,idir_t(2)
    logical                 :: is_nb,found_friend
    integer(i8)             :: ntot,ntot_max,ntot_min,ntot_sum
    integer                 :: isize
    type(MPI_DATATYPE)      :: MPI_INTEGER_I8
    type(MPI_COMM ) :: comm_leaders
    !
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nrank)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,my_block,myid,comm_block)
    call MPI_COMM_RANK(comm_block, myid_block)
    !
    ! generate Cartesian topology
    !
    ! determine lo_g and hi_g from the prescribed decomposition
    !
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,myid_block,myid,comm_leaders)
    if(myid_block == 0) then
      !
      ! gather block size and connectivity information
      !
      allocate(lo_all(3,nblocks),ng_all(3,nblocks),cbc_all(0:1,3,nblocks),bc_all(0:1,3,nblocks),lmin_all(3,nblocks))
      call MPI_ALLGATHER(ng ,3,MPI_INTEGER  ,ng_all ,3,MPI_INTEGER  ,comm_leaders)
      call MPI_ALLGATHER(cbc,6,MPI_CHARACTER,cbc_all,6,MPI_CHARACTER,comm_leaders)
      call MPI_ALLGATHER( bc,6,MPI_REAL_RP  , bc_all,6,MPI_REAL_RP  ,comm_leaders)
      !
      ! determine lower bounds, taking block #1 as reference
      !
      allocate(is_done(nblocks),is_unlocked(nblocks))
      is_done(:)     = .false.
      is_unlocked(:) = .false.
      iblock = 1
      lo_all(:,iblock) = 1
      is_unlocked(iblock) = .true.
      do while(.not.all(is_done))
        if(iblock < 1 .or. iblock > nblocks) then
          iblock = findloc(.not.is_done,.true.,1)
            write(stderr,*) 'ERROR: invalid connectivity for block ', iblock, '.'
            write(stderr,*) ''
            error stop
        endif
        do idir=1,3
          idir_t(:) = pack([1,2,3],[1,2,3] /= idir) ! extract tangential directions
          do inb=0,1
            if(cbc_all(inb,idir,iblock) == 'F') then
              ifriend = nint(bc_all(inb,idir,iblock))
              if(.not.is_unlocked(ifriend)) then
                ! periodicity check
                if(inb == 0 .and. lmin_all(idir,ifriend) > lmin_all(idir,iblock)) cycle
                if(inb == 1 .and. lmin_all(idir,ifriend) < lmin_all(idir,iblock)) cycle
                ! normal direction
                ! inb = 0 -> shift - ng_all(idir,ifriend) | inb = 1 -> shift + ng_all(idir,iblock )
                ioffset = -(1-inb)*ng_all(idir,ifriend) + inb*ng_all(idir,iblock)
                lo_all(idir,ifriend) = lo_all(idir,iblock) + ioffset
                ! tangential directions
                lo_all(idir_t(:),ifriend) = lo_all(idir_t(:),iblock)
                is_unlocked(ifriend) = .true.
              end if
            end if
          end do
        end do
        is_done(iblock) = .true.
        iblock = findloc(.not.is_done.and.is_unlocked,.true.,1)
      end do
      lo_g(:) = lo_all(:,my_block)
      deallocate(lo_all,ng_all,cbc_all,bc_all,lmin_all)
    endif
    call MPI_BCAST(lo_g,3,MPI_INTEGER,0,comm_block)
    hi_g(:) = lo_g(:)+ng(:)-1
    !
    ! determine array extents for possibly uneven data
    !
    coords(:) = get_coords(myid-id_first,dims)
    nb(:,:) = MPI_PROC_NULL
    do idir=1,3
      if(coords(idir)-1 >= 0           ) &
        nb(0,idir) = id_first + get_id([coords(1)-eye(1,idir),coords(2)-eye(2,idir),coords(3)-eye(3,idir)],dims(:))
      if(coords(idir)+1 <= dims(idir)-1) &
        nb(1,idir) = id_first + get_id([coords(1)+eye(1,idir),coords(2)+eye(2,idir),coords(3)+eye(3,idir)],dims(:))
    end do
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
             lmin_all(3,0:nrank-1),lmax_all(3,0:nrank-1),gr_all(3,0:nrank-1),gt_all(3,0:nrank-1), &
             blocks_all(0:nrank-1),cbc_all(0:1,3,0:nrank-1),bc_all(0:1,3,0:nrank-1))
    call MPI_ALLGATHER(lo      ,3,MPI_INTEGER,lo_all    ,3,MPI_INTEGER,MPI_COMM_WORLD)
    call MPI_ALLGATHER(hi      ,3,MPI_INTEGER,hi_all    ,3,MPI_INTEGER,MPI_COMM_WORLD)
    call MPI_ALLGATHER(lo_g    ,3,MPI_INTEGER,lo_g_all  ,3,MPI_INTEGER,MPI_COMM_WORLD)
    call MPI_ALLGATHER(hi_g    ,3,MPI_INTEGER,hi_g_all  ,3,MPI_INTEGER,MPI_COMM_WORLD)
    call MPI_ALLGATHER(lmin    ,3,MPI_REAL_RP,lmin_all  ,3,MPI_REAL_RP,MPI_COMM_WORLD)
    call MPI_ALLGATHER(lmax    ,3,MPI_REAL_RP,lmax_all  ,3,MPI_REAL_RP,MPI_COMM_WORLD)
    call MPI_ALLGATHER(gr      ,3,MPI_REAL_RP,gr_all    ,3,MPI_REAL_RP,MPI_COMM_WORLD)
    call MPI_ALLGATHER(gt      ,3,MPI_INTEGER,gt_all    ,3,MPI_INTEGER,MPI_COMM_WORLD)
    call MPI_ALLGATHER(my_block,1,MPI_INTEGER,blocks_all,1,MPI_INTEGER,MPI_COMM_WORLD)
    call MPI_ALLGATHER(cbc,6,MPI_CHARACTER,cbc_all,6,MPI_CHARACTER,MPI_COMM_WORLD)
    call MPI_ALLGATHER( bc,6,MPI_REAL_RP  , bc_all,6,MPI_REAL_RP  ,MPI_COMM_WORLD)
    do idir=1,3
      do inb=0,1
        if(cbc(inb,idir) == 'F') then
          if(nint(bc(inb,idir)) < 1 .or. nint(bc(inb,idir)) > blocks_all(nrank-1))  then
            write(stderr,*) 'ERROR: invalid connectivity for block ', my_block, '.'
            write(stderr,*) 'Block ', nint(bc(inb,idir)), ' does not exist.'
            write(stderr,*) ''
            error stop
          end if
        end if
        if(nb(inb,idir) == MPI_PROC_NULL.and.cbc(inb,idir) == 'F') then
          is_nb = .false.
          do irank=0,nrank-1
            if(cbc(inb,idir) == cbc_all(1-inb,idir,irank).and.blocks_all(myid) == nint(bc_all(1-inb,idir,irank))) then
              is_nb = .true.
              do iidir = 1,3 ! ensure that the sub blocks blocks share the same extents and grid in the two directions normal to the boundary
                if(iidir /= idir) then
                  is_nb = is_nb.and. &
                  lo(iidir) == lo_all(iidir,irank).and. &
                  hi(iidir) == hi_all(iidir,irank).and. &
                  lmin(iidir) == lmin_all(iidir,irank).and. &
                  lmax(iidir) == lmax_all(iidir,irank).and. &
                  gr(iidir) == gr_all(iidir,irank).and. &
                  gt(iidir) == gt_all(iidir,irank)
                end if
              end do
              !
              ! these checks are not really needed because we specify the links
              ! in the block file, but they may be useful in future and to trap
              ! inconsistent in inputs
              !
              if(is_nb) then
                if(      inb == 0 .and. lo(idir) == lo_g(idir)) then
                  if(     lo(idir) == hi_all(idir,irank)+1) then
                    nb(inb,idir) = irank
                  else if(lo(idir) <  hi_all(idir,irank)+1) then
                    if(periods(idir) == hi_all(idir,irank)-lo(idir)+1) then
                      nb(inb,idir) = irank
                    else if(hi_all(idir,irank) == hi_g_all(idir,irank)) then
                      write(stderr,*) 'ERROR: Inconsistent periodic boundary condition (?).'
                      write(stderr,*) 'Expected: periods(',idir,') = ',hi_all(idir,irank)-lo(idir)+1
                      write(stderr,*) 'Found   : periods(',idir,') = ',periods(idir)
                      write(stderr,*) 'Error when connecting block: ' ,my_block
                      write(stderr,*) ''
                      error stop
                    end if
                  end if
                else if (inb == 1 .and. hi(idir) == hi_g(idir)) then
                  if(     hi(idir) == lo_all(idir,irank)-1) then
                    nb(inb,idir) = irank
                  else if(hi(idir) >  lo_all(idir,irank)-1) then
                    if(periods(idir) == hi(idir)-lo_all(idir,irank)+1) then
                      nb(inb,idir) = irank
                    else if(lo_all(idir,irank) == lo_g_all(idir,irank)) then
                      write(stderr,*) 'ERROR: Inconsistent periodic boundary condition (?).'
                      write(stderr,*) 'Expected: periods(',idir,') = ',hi(idir)-lo_all(idir,irank)+1
                      write(stderr,*) 'Found   : periods(',idir,') = ',periods(idir)
                      write(stderr,*) 'Error when connecting block: ' ,my_block
                      write(stderr,*) ''
                      error stop
                    end if
                  end if
                else
                  write(stderr,*) 'ERROR: Expected connectivity between blocks',my_block,' and ',blocks_all(irank),'not found.'
                  write(stderr,*) ''
                  error stop
                end if
              end if
            end if
          end do
        end if
      end do
    end do
    deallocate(lo_all,hi_all,lmin_all,lmax_all,gr_all,gt_all,blocks_all)
    do idir=1,3
      do inb=0,1
        if(cbc(inb,idir) == 'F'.and.nb(inb,idir) == MPI_PROC_NULL) then
          write(stderr,*) 'ERROR: Expected connectivity between blocks',my_block,' and ',nint(bc(inb,idir)), ' is not possible.'
          write(stderr,*) 'E.g. Blocks must share the same boundaries.'
          write(stderr,*) ''
          error stop
        end if
      end do
    end do
    !
    do idir=1,3
      is_bound(:,idir) = .false.
      where(nb(0:1,idir) == MPI_PROC_NULL) is_bound(0:1,idir) = .true.
      call makehalo(idir,1,n,halos(idir))
    end do
    !
    ! check distribution of grid points over the different tasks
    !
    ntot = product(hi(:)-lo(:)+1_i8)
    call MPI_SIZEOF(1_i8,isize)
    call MPI_TYPE_MATCH_SIZE(MPI_TYPECLASS_INTEGER,isize,MPI_INTEGER_I8)
    call MPI_ALLREDUCE(ntot,ntot_min,1,MPI_INTEGER_I8,MPI_MIN,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(ntot,ntot_max,1,MPI_INTEGER_I8,MPI_MAX,MPI_COMM_WORLD)
    call MPI_ALLREDUCE(ntot,ntot_sum,1,MPI_INTEGER_I8,MPI_SUM,MPI_COMM_WORLD)
    write(stdout,*) 'Maximum, minimum, average, and normalized average number of grid points for task ',myid, &
                    '(block',my_block,'): ',ntot_max, &
                                            ntot_min, &
                                            (1._rp*ntot_sum)/(1._rp*nrank), &
                                            (1._rp*ntot*nrank)/(1._rp*ntot_sum)
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
  function get_id(coords,dims,periods) result(id)
    use mpi_f08, only: MPI_PROC_NULL
    implicit none
    integer :: id
    integer, intent(in), dimension(3) :: coords,dims
    logical, intent(in), dimension(3), optional :: periods
    integer, dimension(3) :: coords_aux,shift
    coords_aux(:) = coords(:)
    if(present(periods)) then
      shift(:) = 0
      where(periods(:))
        where(coords_aux(:)>dims(:)-1) shift(:) =   (0        -coords_aux(:))/dims(:)
        where(coords_aux(:)<0        ) shift(:) =   (dims(:)-1-coords_aux(:))/dims(:)
        coords_aux(:) = coords_aux(:) + shift(:)*dims(:)
      end where
    end if
    if(all(coords_aux(:)<=dims(:)-1).and.all(coords_aux(:)>=0)) then
      id = coords_aux(1)+coords_aux(2)*dims(1)+coords_aux(3)*dims(1)*dims(2)
    else
      id = MPI_PROC_NULL
    end if
  end function get_id
  function get_coords(id,dims) result(coords)
    integer :: coords(3)
    integer, intent(in) :: id, dims(3)
    coords(:) = [mod(id,dims(1)),mod(id/dims(1),dims(2)),mod(id/(dims(1)*dims(2)),dims(3))]
  end function get_coords
  pure integer function eye(i,j)
    integer, intent(in) :: i,j
    eye = (i/j)*(j/i)
  end function eye
end module mod_initmpi
