module mod_initgrid
  use mod_param, only:pi
  use mod_types
  implicit none
  private
  public initgrid,distribute_grid,save_grid
  contains
  subroutine initgrid(n,lo,hi,gt,gr,l,drc_g,drf_g,rc_g,rf_g)
    !
    ! initializes a non-uniform grid
    !
    implicit none
    integer, parameter :: CLUSTER_TWO_END   = 0, &
                          CLUSTER_ONE_END   = 1, &
                          CLUSTER_MIDDLE    = 2, &
                          CLUSTER_ONE_END_R = 3
    integer         , intent(in )                  :: n,lo,hi,gt
    real(rp)        , intent(in )                  :: gr,l
    real(rp)        , intent(out), dimension(1-1:) :: drc_g,drf_g,rc_g,rf_g
    real(rp) :: r0
    integer :: q
    procedure (), pointer :: gridpoint => null()
    select case(gt)
    case(CLUSTER_TWO_END)
      gridpoint => gridpoint_cluster_two_end
    case(CLUSTER_ONE_END)
      gridpoint => gridpoint_cluster_one_end
    case(CLUSTER_ONE_END_R)
      gridpoint => gridpoint_cluster_one_end_r
    case(CLUSTER_MIDDLE)
      gridpoint => gridpoint_cluster_middle
    case default
      gridpoint => gridpoint_cluster_two_end
    end select
    !
    ! step 1) determine coordinates of cell faces rf
    !
    do q=1,n
      r0  = (q-0._rp)/(1._rp*n)
      call gridpoint(gr,r0,rf_g(q))
      rf_g(q) = rf_g(q)*l
    enddo
    rf_g(0) = 0._rp
    !
    ! step 2) determine grid spacing between faces drf
    !
    do q=1,n
      drf_g(q) = rf_g(q)-rf_g(q-1)
    enddo
    drf_g(0  ) = drf_g(1)
    drf_g(n+1) = drf_g(n)
    !
    ! step 3) determine grid spacing between centers drc
    !
    do q=0,n
      drc_g(q) = .5_rp*(drf_g(q)+drf_g(q+1))
    enddo
    drc_g(n+1) = drc_g(n)
    !
    ! step 4) compute coordinates of cell centers rc and faces rf
    !
    rc_g(0)    = -drc_g(0)/2._rp
    rf_g(0)    = 0._rp
    do q=1,n+1
      rc_g(q) = rc_g(q-1) + drc_g(q-1)
      rf_g(q) = rf_g(q-1) + drf_g(q)
    enddo
  end subroutine initgrid
  !
  subroutine distribute_grid(lo,hi,grid_g,grid)
    implicit none
    integer  :: lo,hi
    real(rp), intent(in ), dimension(1-1 :) :: grid_g
    real(rp), intent(out), dimension(lo-1:) :: grid
    grid(lo-1:hi+1) = grid_g(lo-1:hi+1)
  end subroutine distribute_grid
  !
  subroutine save_grid(fname,ng,rf_g,rc_g,drf_g,drc_g)
    implicit none
    character(len=*), intent(in) :: fname
    integer         , intent(in) :: ng
    real(rp)        , intent(in), dimension(1-1:) :: rf_g,rc_g,drf_g,drc_g
    integer :: iunit,q,rlen
    !
    inquire(iolength=rlen) 1._rp
    open(newunit=iunit,file=trim(fname)//'.bin',status='replace',access='direct',recl=ng*4*rlen)
    write(iunit,rec=1) rf_g(1:ng),rc_g(1:ng),drf_g(1:ng),drc_g(1:ng)
    close(iunit)
    open(newunit=iunit,status='replace',file=trim(fname)//'.out')
    do q=0,ng+1
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
end module mod_initgrid
