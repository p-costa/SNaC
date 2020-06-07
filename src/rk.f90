module mod_rk
  use mod_types
  implicit none
  private
  public rk_mom,rk_scal
  contains
  subroutine rk_mom(rkpar,lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,l,dt,bforce,is_forced,velf,visc,u,v,w,p,dudtrko,dvdtrko,dwdtrko,up,vp,wp,f)
    use mod_mom  , only: momx_a,momy_a,momz_a,momx_d,momy_d,momz_d,momx_p,momy_p,momz_p
    use mod_debug, only: chkmean
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: lo,hi
    real(rp), intent(in   ), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in   ), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in   ), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in   ), dimension(3) :: l
    real(rp), intent(in   )               :: visc,dt
    real(rp), intent(in   ), dimension(3) :: bforce
    logical , intent(in   ), dimension(3) :: is_forced
    real(rp), intent(in   ), dimension(3) :: velf
    real(rp), intent(in   ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u ,v ,w ,p
    real(rp), intent(inout), dimension(lo(1):  ,lo(2):  ,lo(3):  ) :: dudtrko,dvdtrko,dwdtrko
    real(rp), intent(out  ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: up,vp,wp
    real(rp), intent(out  ), dimension(3) :: f
    real(rp),                dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) :: dudtrk ,dvdtrk ,dwdtrk
#ifdef _IMPDIFF
    real(rp),                dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) :: dudtrkd,dvdtrkd,dwdtrkd
#endif
    real(rp) :: factor1,factor2,factor12
    integer  :: i,j,k
    real(rp) :: mean
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the momentum equations.
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    factor12 = factor1 + factor2
    !$OMP WORKSHARE
    dudtrk(:,:,:) = 0._rp
    dvdtrk(:,:,:) = 0._rp
    dwdtrk(:,:,:) = 0._rp
    !$OMP END WORKSHARE
    call momx_d(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,visc,u,dudtrk)
    call momy_d(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,visc,v,dvdtrk)
    call momz_d(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,visc,w,dwdtrk)
#ifdef _IMPDIFF
    dudtrkd(:,:,:) = dudtrk(:,:,:)
    dvdtrkd(:,:,:) = dvdtrk(:,:,:)
    dwdtrkd(:,:,:) = dwdtrk(:,:,:)
#endif
    call momx_a(lo,hi,dxf,dyf,dzf,u,v,w,dudtrk)
    call momy_a(lo,hi,dxf,dyf,dzf,u,v,w,dvdtrk)
    call momz_a(lo,hi,dxf,dyf,dzf,u,v,w,dwdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(hi,lo,factor1,factor2,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk,dudtrko,dvdtrko,dwdtrko)
    do k=hi(3),lo(3)
      do j=hi(2),lo(2)
        do i=hi(1),lo(1)
          ! could be split in two loops, because factor2=0 for istep=1, but like this reads nicer
          up(i,j,k) = u(i,j,k) + factor1*dudtrk(i,j,k) + factor2*dudtrko(i,j,k)
          vp(i,j,k) = v(i,j,k) + factor1*dvdtrk(i,j,k) + factor2*dvdtrko(i,j,k)
          wp(i,j,k) = w(i,j,k) + factor1*dwdtrk(i,j,k) + factor2*dwdtrko(i,j,k)
          dudtrko(i,j,k) = dudtrk(i,j,k)
          dvdtrko(i,j,k) = dvdtrk(i,j,k)
          dwdtrko(i,j,k) = dwdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP WORKSHARE
    dudtrk(:,:,:) = 0._rp
    dvdtrk(:,:,:) = 0._rp
    dwdtrk(:,:,:) = 0._rp
    !$OMP END WORKSHARE
    call momx_p(lo,hi,dxc,bforce(1),p,dudtrk)
    call momy_p(lo,hi,dyc,bforce(2),p,dvdtrk)
    call momz_p(lo,hi,dzc,bforce(3),p,dwdtrk) ! we could perform the pressure gradient calculation in the loop below instead, but I
                                              ! decided to have more modular, like this, for simplicity.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,factor12,u,v,w,up,vp,wp,dudtrk,dvdtrk,dwdtrk)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          up(i,j,k) = up(i,j,k) + factor12*dudtrk(i,j,k)
          vp(i,j,k) = vp(i,j,k) + factor12*dvdtrk(i,j,k)
          wp(i,j,k) = wp(i,j,k) + factor12*dwdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !
    ! bulk velocity forcing
    !
    f(:) = 0._rp
    if(is_forced(1)) then
      call chkmean(lo,hi,l,dxc,dyf,dzf,up,mean)
      f(1) = velf(1) - mean
      !$OMP WORKSHARE
      up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = up(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) + f(1)
      !$OMP END WORKSHARE
    endif
    if(is_forced(2)) then
      call chkmean(lo,hi,l,dxf,dyc,dzf,vp,mean)
      f(2) = velf(2) - mean
      !$OMP WORKSHARE
      vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = vp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) + f(2)
      !$OMP END WORKSHARE
    endif
    if(is_forced(3)) then
      call chkmean(lo,hi,l,dxf,dyf,dzc,wp,mean)
      f(3) = velf(3) - mean
      !$OMP WORKSHARE
      wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = wp(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) + f(3)
      !$OMP END WORKSHARE
    endif
#ifdef _IMPDIFF
    !
    ! compute rhs of helmholtz equation
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,factor12,factor2,visc,up,vp,wp,dudtrkd,dvdtrkd,dwdtrkd)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          up(i,j,k) = up(i,j,k) - .5_rp*factor12*dudtrkd(i,j,k)
          vp(i,j,k) = vp(i,j,k) - .5_rp*factor12*dvdtrkd(i,j,k)
          wp(i,j,k) = wp(i,j,k) - .5_rp*factor12*dwdtrkd(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
#endif
    return
  end subroutine rk_mom
  subroutine rk_scal(rkpar,lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,dt,alpha,u,v,w,dsdtrko,s)
    use mod_scal, only: scal_a,scal_d
    !
    ! low-storage 3rd-order Runge-Kutta scheme 
    ! for time integration of the scalar field.
    !
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    integer , intent(in   ), dimension(3) :: lo,hi
    real(rp), intent(in   ), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in   ), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in   ), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in   )                      :: dt,alpha
    real(rp), intent(in   ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(inout), dimension(lo(1):  ,lo(2):  ,lo(3):  ) :: dsdtrko
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: s
    real(rp),                dimension(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) :: dsdtrk
    real(rp) :: factor1,factor2
    integer  :: i,j,k
    !
    factor1 = rkpar(1)*dt
    factor2 = rkpar(2)*dt
    dsdtrk(:,:,:) = 0._rp
    call scal_d(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,s,dsdtrk)
    call scal_a(lo,hi,dxf,dyf,dzf,u,v,w,s,dsdtrk)
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,factor1,factor2,s,dsdtrk,dsdtrko)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          s(i,j,k) = s(i,j,k) + factor1*dsdtrk(i,j,k) + factor2*dsdtrko(i,j,k)
          dsdtrko(i,j,k) = dsdtrk(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine rk_scal
end module mod_rk
