module mod_initgrid
  use mod_param, only:pi
  use mod_types
  implicit none
  private
  public initgrid
  contains
  subroutine initgrid(n,lo,hi,gt,gr,l,drc,drf,rc,rf,drc_g,drf_g,rc_g,rf_g,fname)
    !
    ! initializes a non-uniform grid
    !
    implicit none
    integer         , intent(in )                   :: n,lo,hi,gt
    real(rp)        , intent(in )                   :: gr,l
    character(len=*), intent(in )                   :: fname
    real(rp)        , intent(out), dimension(lo-1:) :: drc  ,drf  ,rc  ,rf
    real(rp)        , intent(out), dimension(0:   ) :: drc_g,drf_g,rc_g,rf_g
    real(rp) :: r0
    integer :: q
    procedure (), pointer :: gridpoint => null()
    select case(gt)
    case(0)
      gridpoint => gridpoint_cluster_two_end
    case(1)
      gridpoint => gridpoint_cluster_middle
    case(2)
      gridpoint => gridpoint_cluster_one_end
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
    !
    ! step 5) extract portion of the global grid pertaining 
    !         to the specific task
    !
    call distribute_grid(lo,hi,drc_g,drc)
    call distribute_grid(lo,hi,drf_g,drf)
    call distribute_grid(lo,hi, rc_g, rc)
    call distribute_grid(lo,hi, rf_g, rf)
    !
    return
  end subroutine initgrid
  subroutine distribute_grid(lo,hi,grid_g,grid)
    implicit none
    integer  :: lo,hi
    real(rp), intent(in ), dimension(0   :) :: grid_g
    real(rp), intent(out), dimension(lo-1:) :: grid
    grid(lo-1:hi+1) = grid_g(lo-1:hi+1)
    return
  end subroutine distribute_grid
  subroutine save_grid(fname,ng,rf,rc,drf,drc)
    implicit none
    character(len=*), intent(in) :: fname
    integer         , intent(in) :: ng
    real(rp)        , intent(in), dimension(0:) :: rf,rc,drf,drc
    integer :: iunit
    open(newunit=iunit,file=trim(fname),access='direct',recl=4*ng*sizeof(1._rp))
    write(iunit,rec=1) rf(1:ng),rc(1:ng),drf(1:ng),drc(1:ng)
    close(iunit)
    return
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
    if(alpha.ne.0._rp) then
      r = 0.5_rp*(1._rp+tanh((r0-0.5_rp)*alpha)/tanh(alpha/2._rp))
      !r = 0.5_rp*(1._rp+erf( (r0-0.5_rp)*alpha)/erf( alpha/2._rp))
    else
      r = r0
    endif
    return
  end subroutine gridpoint_cluster_two_end
  subroutine gridpoint_cluster_one_end(alpha,r0,r)
    !
    ! clustered at the lower side
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha.ne.0._rp) then
      r = 1.0_rp*(1._rp+tanh((r0-1.0_rp)*alpha)/tanh(alpha/1._rp))
      !r = 1.0_rp*(1._rp+erf( (r0-1.0_rp)*alpha)/erf( alpha/1._rp))
    else
      r = r0
    endif
    return
  end subroutine gridpoint_cluster_one_end
  subroutine gridpoint_cluster_middle(alpha,r0,r)
    !
    ! clustered in the middle
    !
    implicit none
    real(rp), intent(in ) :: alpha,r0
    real(rp), intent(out) :: r
    if(alpha.ne.0._rp) then
      if(    r0.le.0.5_rp) then 
        r = 0.5_rp*(1._rp-1._rp+tanh(2._rp*alpha*(r0-0._rp))/tanh(alpha))
        !r = 0.5_rp*(1._rp-1._rp+erf( 2._rp*alpha*(r0-0._rp))/erf( alpha))
      elseif(r0.gt.0.5_rp) then
        r = 0.5_rp*(1._rp+1._rp+tanh(2._rp*alpha*(r0-1._rp))/tanh(alpha))
        !r = 0.5_rp*(1._rp+1._rp+erf( 2._rp*alpha*(r0-1._rp))/erf( alpha))
      endif
    else
      r = r0
    endif
    return
  end subroutine gridpoint_cluster_middle
end module mod_initgrid
