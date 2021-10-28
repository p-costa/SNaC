module mod_fillps
  use mod_types
  implicit none
  private
  public fillps
  contains
  subroutine fillps(lo,hi,dxf,dyf,dzf,dt,up,vp,wp,p)
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
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: up,vp,wp
    real(rp), intent(out), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,p,up,vp,wp,dt,dxf,dyf,dzf) &
    !$OMP PRIVATE(i,j,k)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = ( &
                      (wp(i,j,k)-wp(i,j,k-1))/dzf(k) + &
                      (vp(i,j,k)-vp(i,j-1,k))/dyf(j) + &
                      (up(i,j,k)-up(i-1,j,k))/dxf(i) &
                     )/dt
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine fillps
end module mod_fillps
