module mod_updt_pressure
  use mod_types
  implicit none
  private
  public updt_pressure
  contains
  subroutine updt_pressure(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,pp,p)
    !
    ! calculates the final pressure field
    !
    implicit none
    integer , intent(in   ), dimension(3       ) :: lo,hi
    real(rp), intent(in   ), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in   ), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in   ), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in   ) :: alpha
    real(rp), intent(in   ) , dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: pp
    real(rp), intent(inout) , dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,p,pp)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = p(i,j,k) + pp(i,j,k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
#ifdef _IMPDIFF
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k)            &
    !$OMP SHARED(lo,hi,p,pp,dxc,dxf,dyc,dyf,dzc,dzf,alpha)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = p(i,j,k) + &
                      alpha*( ((pp(i+1,j,k)-pp(i  ,j,k))/dxc(i  ) - &
                               (pp(i  ,j,k)-pp(i-1,j,k))/dxc(i-1))/dxf(i) + &
                              ((pp(i,j+1,k)-pp(i,j  ,k))/dyc(j  ) - &
                               (pp(i,j  ,k)-pp(i,j-1,k))/dyc(j-1))/dyf(j) + &
                              ((pp(i,j,k+1)-pp(i,j,k  ))/dzc(k  ) - &
                               (pp(i,j,k  )-pp(i,j,k-1))/dzc(k-1))/dzf(k) )
        end do
      end do
    end do
    !$OMP END PARALLEL DO
#endif
  end subroutine updt_pressure
end module mod_updt_pressure
