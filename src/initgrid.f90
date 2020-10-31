module mod_initgrid
  use mpi
  use mod_common_mpi, only:ierr
  use mod_param     , only:pi
  use mod_types
  implicit none
  private
  public initgrid,distribute_grid,bound_grid,save_grid
  contains
  subroutine initgrid(lo,hi,gt,gr,lmin,lmax,drc_g,drf_g,rc_g,rf_g)
    !
    ! initializes a non-uniform grid
    !
    ! NEEDS TO BE PROBABLY SLIGHTLY ADAPTED IN CASE IT IS IN THE MIDDLE OF A BLOCK
    !
    implicit none
    integer, parameter :: CLUSTER_TWO_END              = 0, &
                          CLUSTER_ONE_END              = 1, &
                          CLUSTER_MIDDLE               = 2, &
                          CLUSTER_ONE_END_R            = 3, &
                          CLUSTER_GEOMETRIC_ONE_END    = 4, &
                          CLUSTER_GEOMETRIC_ONE_END_R  = 5, &
                          CLUSTER_GEOMETRIC_TWO_ENDS   = 6, &
                          CLUSTER_GEOMETRIC_MIDDLE     = 7
    integer         , intent(in )                  :: lo,hi,gt
    real(rp)        , intent(in )                  :: gr,lmin,lmax
    real(rp)        , intent(out), dimension(lo-1:) :: drc_g,drf_g,rc_g,rf_g
    real(rp) :: r0
    integer :: q,n
    procedure (), pointer :: gridpoint => null()
    select case(gt)
    case(CLUSTER_TWO_END)
      gridpoint => gridpoint_cluster_two_end
    case(CLUSTER_ONE_END)
      gridpoint => gridpoint_cluster_one_end
    case(CLUSTER_MIDDLE)
      gridpoint => gridpoint_cluster_middle
    case(CLUSTER_ONE_END_R)
      gridpoint => gridpoint_cluster_one_end_r
    case(CLUSTER_GEOMETRIC_ONE_END)
     gridpoint => gridpoint_cluster_geometric_one_end
    CASE(CLUSTER_GEOMETRIC_ONE_END_R)
     gridpoint => gridpoint_cluster_geometric_one_end_r
    CASE(CLUSTER_GEOMETRIC_TWO_ENDS)
     gridpoint => gridpoint_cluster_geometric_two_ends
    CASE(CLUSTER_GEOMETRIC_MIDDLE)
     gridpoint => gridpoint_cluster_geometric_middle
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces rf
    !
    do q=lo,hi
      n  = hi-lo+1
      r0 = (q-lo+1-0._rp)/(1._rp*n)
      call gridpoint(gr,r0,rf_g(q))
      rf_g(q) = lmin + rf_g(q)*(lmax-lmin)
    enddo
    rf_g(lo-1) = lmin
    !
    ! step 2) determine grid spacing between faces drf
    !
    do q=lo,hi
      drf_g(q) = rf_g(q)-rf_g(q-1)
    enddo
    drf_g(lo-1) = drf_g(lo)
    drf_g(hi+1) = drf_g(hi)
    !
    ! step 3) determine grid spacing between centers drc
    !
    do q=lo-1,hi
      drc_g(q) = .5_rp*(drf_g(q)+drf_g(q+1))
    enddo
    drc_g(hi+1) = drc_g(hi)
    !
    ! step 4) compute coordinates of cell centers rc and faces rf
    !
    rc_g(lo-1)    = lmin - drc_g(lo-1)/2._rp
    rf_g(lo-1)    = lmin
    do q=lo,hi+1
      rc_g(q) = rc_g(q-1) + drc_g(q-1)
      rf_g(q) = rf_g(q-1) + drf_g(q  )
    enddo
    !
  end subroutine initgrid
  subroutine distribute_grid(lo_g,lo,hi,grid_g,grid)
    implicit none
    integer , intent(in )                     :: lo_g,lo,hi
    real(rp), intent(in ), dimension(lo_g-1:) :: grid_g
    real(rp), intent(out), dimension(lo  -1:) :: grid
    grid(lo-1:hi+1) = grid_g(lo-1:hi+1)
  end subroutine distribute_grid
  subroutine bound_grid(lo_g,hi_g,lo,hi,nb,is_periodic,lo_min,hi_max,grid_f,grid_c)
    implicit none
    integer , intent(in   )                   :: lo_g,hi_g,lo,hi
    integer , intent(in   ), dimension(0:1)   :: nb
    logical , intent(in   )                   :: is_periodic
    integer , intent(in   )                   :: lo_min,hi_max
    real(rp), intent(inout), dimension(lo-1:) :: grid_f,grid_c
    integer                                   :: shift
    integer , dimension(4)                    :: requests
    integer                                   :: nrequests
    nrequests = 0
    shift = 0
    if(lo == lo_g) then
      if(is_periodic.and.lo == lo_min) shift = hi_max-lo_min + 1
      call MPI_ISEND(grid_f(lo  ),1,MPI_REAL_RP,nb(0),lo_g+shift, &
                     MPI_COMM_WORLD,requests(nrequests+1),ierr)
      call MPI_IRECV(grid_f(lo-1),1,MPI_REAL_RP,nb(0),lo_g-1    , &
                     MPI_COMM_WORLD,requests(nrequests+2),ierr)
      nrequests = nrequests + 2
    endif
    shift = 0
    if(hi == hi_g) then
      if(is_periodic.and.hi == hi_max) shift = hi_max-lo_min + 1
      call MPI_ISEND(grid_f(hi  ),1,MPI_REAL_RP,nb(1),hi_g-shift, &
                     MPI_COMM_WORLD,requests(nrequests+1),ierr)
      call MPI_IRECV(grid_f(hi+1),1,MPI_REAL_RP,nb(1),hi_g+1    , &
                     MPI_COMM_WORLD,requests(nrequests+2),ierr)
      nrequests = nrequests + 2
    endif
    call MPI_WAITALL(nrequests,requests,MPI_STATUSES_IGNORE,ierr)
    if(lo == lo_g) grid_c(lo-1) = (grid_f(lo  )+grid_f(lo-1))/2._rp
    if(hi == hi_g) grid_c(hi  ) = (grid_f(hi+1)+grid_f(hi  ))/2._rp
    if(hi == hi_g) grid_c(hi+1) = grid_c(hi) ! not needed
  end subroutine bound_grid
  subroutine save_grid(fname,lo_g,hi_g,rf_g,rc_g,drf_g,drc_g)
    implicit none
    character(len=*), intent(in) :: fname
    integer         , intent(in) :: lo_g,hi_g
    real(rp)        , intent(in), dimension(lo_g-1:) :: rf_g,rc_g,drf_g,drc_g
    integer :: iunit,q,reclen
    inquire(iolength=reclen) rf_g(lo_g:hi_g),rc_g(lo_g:hi_g),drf_g(lo_g:hi_g),drc_g(lo_g:hi_g)
    open(newunit=iunit,file=trim(fname)//'.bin',status='replace',access='direct',recl=reclen)
    write(iunit,rec=1) rf_g(lo_g:hi_g),rc_g(lo_g:hi_g),drf_g(lo_g:hi_g),drc_g(lo_g:hi_g)
    close(iunit)
    open(newunit=iunit,status='replace',file=trim(fname)//'.out')
    do q=lo_g-1,hi_g+1
      write(iunit,'(5E15.7)') 0._rp,rf_g(q),rc_g(q),drf_g(q),drc_g(q)
    enddo
    close(iunit)
  end subroutine save_grid
  !
  ! grid stretching functions 
  ! see e.g., Fluid Flow Phenomena -- A Numerical Toolkit, by P. Orlandi
  !           Pirozzoli et al. JFM 788, 614–639 (commented)
  !
  subroutine gridpoint_cluster_two_end(alpha,r0,r)
    !
    ! clustered at the two sides
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha /= 0._rp) then
      r = 0.5_rp*(1._rp+tanh((r0-0.5_rp)*alpha)/tanh(alpha/2._rp))
      !r = 0.5_rp*(1._rp+erf( (r0-0.5_rp)*alpha)/erf( alpha/2._rp))
    else
      r = r0
    endif
  end subroutine gridpoint_cluster_two_end
  subroutine gridpoint_cluster_one_end(alpha,r0,r)
    !
    ! clustered at the lower side
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha /= 0._rp) then
      r = 1.0_rp*(1._rp+tanh((r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1.0_rp*(1._rp+erf( (r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    endif
  end subroutine gridpoint_cluster_one_end
  subroutine gridpoint_cluster_one_end_r(alpha,r0,r)
    !
    ! clustered at the lower side
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha /= 0._rp) then
      r = 1._rp-1.0_rp*(1._rp+tanh((1._rp-r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1._rp-1.0_rp*(1._rp+erf( (1._rp-r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    endif
  end subroutine gridpoint_cluster_one_end_r
  subroutine gridpoint_cluster_middle(alpha,r0,r)
    !
    ! clustered in the middle
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha /= 0._rp) then
      if(    r0 <= 0.5_rp) then 
        r = 0.5_rp*(1._rp-1._rp+tanh(2._rp*alpha*(r0-0._rp))/tanh(alpha))
        !r = 0.5_rp*(1._rp-1._rp+erf( 2._rp*alpha*(r0-0._rp))/erf( alpha))
      elseif(r0 >  0.5_rp) then
        r = 0.5_rp*(1._rp+1._rp+tanh(2._rp*alpha*(r0-1._rp))/tanh(alpha))
        !r = 0.5_rp*(1._rp+1._rp+erf( 2._rp*alpha*(r0-1._rp))/erf( alpha))
      endif
    else
      r = r0
    endif
  end subroutine gridpoint_cluster_middle
  subroutine gridpoint_cluster_geometric_one_end(alpha,r0,r)
    !
    ! geometric progression
    !
    implicit none
    !real(rp), parameter   :: power = 1._rp/2._rp
    !real(rp), parameter   :: power = 1._rp/1._rp
    real(rp), parameter   :: power = 1._rp/3._rp
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(    r0 == 1.0_rp) then
      r = r0
    else
      r = r0*((1._rp-r0**(alpha))/(1._rp-r0)/alpha)**power
    endif
  end subroutine gridpoint_cluster_geometric_one_end
  subroutine gridpoint_cluster_geometric_one_end_r(alpha,r0,r)
    !
    ! reversed geometric progression
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    call gridpoint_cluster_geometric_one_end(alpha,1._rp-r0,r)
    r = 1._rp-r
  end subroutine gridpoint_cluster_geometric_one_end_r
  subroutine gridpoint_cluster_geometric_two_ends(alpha,r0,r)
    !
    ! geometric progression towards each end
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(r0 <= 0.5) then
      call gridpoint_cluster_geometric_one_end(alpha/2._rp,2._rp*r0        ,r)
      r = r/2._rp
    else
      call gridpoint_cluster_geometric_one_end(alpha/2._rp,2._rp*(1._rp-r0),r)
      r = 1._rp-r/2._rp
    endif
  end subroutine gridpoint_cluster_geometric_two_ends
  subroutine gridpoint_cluster_geometric_middle(alpha,r0,r)
    !
    ! geometric progression towards each end
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(r0 <= 0.5) then
      call gridpoint_cluster_geometric_one_end(alpha/2._rp,(.5_rp-r0)*2._rp,r)
      r = (1._rp-r)/2._rp
    else
      call gridpoint_cluster_geometric_one_end(alpha/2._rp,(r0-.5_rp)*2._rp,r)
      r = (1._rp+r)/2._rp
    endif
  end subroutine gridpoint_cluster_geometric_middle
end module mod_initgrid
