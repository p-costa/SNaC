module mod_correc
  use mod_types
  implicit none
  private
  public correc
  contains
  subroutine correc(lo,hi,dxc,dyc,dzc,dt,p,up,vp,wp,u,v,w)
    !
    ! corrects the velocity so that it is divergence free
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ), dimension(lo(1)-1:) :: dxc
    real(rp), intent(in ), dimension(lo(2)-1:) :: dyc
    real(rp), intent(in ), dimension(lo(3)-1:) :: dzc
    real(rp), intent(in ) :: dt
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p,up,vp,wp
    real(rp), intent(out), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dt,dxc,u,up,p) &
    !$OMP PRIVATE(i,j,k)
    do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)
          u(i,j,k) = up(i,j,k) - dt*(p(i+1,j,k)-p(i,j,k))/dxc(i)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dt,dyc,v,vp,p) &
    !$OMP PRIVATE(i,j,k)
    do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)
        do i=lo(1)-1,hi(1)+1
          v(i,j,k) = vp(i,j,k) - dt*(p(i,j+1,k)-p(i,j,k))/dyc(j)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dt,dzc,w,wp,p) &
    !$OMP PRIVATE(i,j,k)
    do k=lo(3)-1,hi(3)
      do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
          w(i,j,k) = wp(i,j,k) - dt*(p(i,j,k+1)-p(i,j,k))/dzc(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine correc
end module mod_correc
