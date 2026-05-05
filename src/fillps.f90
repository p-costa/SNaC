module mod_fillps
  use mod_types
  implicit none
  private
  public fillps
  contains
  subroutine fillps(lo,hi,dxf,dyf,dzf,dt,u,v,w,p)
    !
    !  fill the right-hand side of the Poisson equation for the correction pressure.
    !
    !  the discrete divergence is:
    !
    !  w(i,j,k)-w(i,j,k-1)   v(i,j,k)-v(i,j-1,k)   u(i,j,k)-u(i-1,j,k)
    !  ------------------- + ------------------- + -------------------  = div
    !          dz                    dy                    dx
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in ), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in ), dimension(lo(3)-1:) :: dzf
    real(rp), intent(in ) :: dt
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(out), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,p,u,v,w,dt,dxf,dyf,dzf)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = ( &
#if   defined(_FFT_Z)
                      (w(i,j,k)-w(i,j,k-1))/dzf(k)*dxf(i)*dyf(j) + &
                      (v(i,j,k)-v(i,j-1,k))*dxf(i) + &
                      (u(i,j,k)-u(i-1,j,k))*dyf(j) &
#elif defined(_FFT_Y)
                      (w(i,j,k)-w(i,j,k-1))*dxf(i) + &
                      (v(i,j,k)-v(i,j-1,k))/dyf(j)*dxf(i)*dzf(k) + &
                      (u(i,j,k)-u(i-1,j,k))*dzf(k) &
#elif defined(_FFT_X)
                      (w(i,j,k)-w(i,j,k-1))*dyf(j) + &
                      (v(i,j,k)-v(i,j-1,k))*dzf(k) + &
                      (u(i,j,k)-u(i-1,j,k))/dxf(i)*dyf(j)*dzf(k) &
#else
                      (w(i,j,k)-w(i,j,k-1))*dxf(i)*dyf(j) + &
                      (v(i,j,k)-v(i,j-1,k))*dxf(i)*dzf(k) + &
                      (u(i,j,k)-u(i-1,j,k))*dyf(j)*dzf(k) &
#endif
                     )/dt
        end do
      end do
    end do
  end subroutine fillps
end module mod_fillps
