module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(lo,hi,dxc,dyc,dzc,dt,p,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ), dimension(lo(1)-1:) :: dxc
    real(rp), intent(in ), dimension(lo(2)-1:) :: dyc
    real(rp), intent(in ), dimension(lo(3)-1:) :: dzc
    real(rp), intent(in ) :: dt
    real(rp), intent(in   ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dt,dxc,u,p)
    do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)
          u(i,j,k) = u(i,j,k) - dt*(p(i+1,j,k)-p(i,j,k))/dxc(i)
        end do
      end do
    end do
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dt,dyc,v,p)
    do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)
        do i=lo(1)-1,hi(1)+1
          v(i,j,k) = v(i,j,k) - dt*(p(i,j+1,k)-p(i,j,k))/dyc(j)
        end do
      end do
    end do
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dt,dzc,w,p)
    do k=lo(3)-1,hi(3)
      do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
          w(i,j,k) = w(i,j,k) - dt*(p(i,j,k+1)-p(i,j,k))/dzc(k)
        end do
      end do
    end do
  end subroutine correc
end module mod_correc
