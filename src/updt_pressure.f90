module mod_updt_pressure
  use mod_types
  implicit none
  private
  public updt_pressure
  contains
  subroutine updt_pressure(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,pp,p,alpha_arr)
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
    real(rp), intent(in   ) , dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), optional :: alpha_arr
    integer :: i,j,k
#if defined(_NON_NEWTONIAN) && defined(_IMPDIFF)
    real(rp) :: alphaxm,alphaxp,alphaym,alphayp,alphazm,alphazp
#endif
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,p,pp)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          p(i,j,k) = p(i,j,k) + pp(i,j,k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
#ifdef _IMPDIFF
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k)            &
    !$OMP SHARED(lo,hi,p,pp,dxc,dxf,dyc,dyf,dzc,dzf,alpha)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
#if defined(_NON_NEWTONIAN) && defined(_IMPDIFF)
                      alphaxp = 0.50_rp*(alpha_arr(i,j,k)+alpha_arr(i+1,j,k))
                      alphaxm = 0.50_rp*(alpha_arr(i,j,k)+alpha_arr(i-1,j,k))
                      alphayp = 0.50_rp*(alpha_arr(i,j,k)+alpha_arr(i,j+1,k))
                      alphaym = 0.50_rp*(alpha_arr(i,j,k)+alpha_arr(i,j-1,k))
                      alphazp = 0.50_rp*(alpha_arr(i,j,k)+alpha_arr(i,j,k+1))
                      alphazm = 0.50_rp*(alpha_arr(i,j,k)+alpha_arr(i,j,k-1))
          p(i,j,k) = p(i,j,k) + & 
                      alpha*( (alphaxp*(pp(i+1,j,k)-pp(i  ,j,k))/dxc(i  ) - &
                               alphaxm*(pp(i  ,j,k)-pp(i-1,j,k))/dxc(i-1))/dxf(i) + &
                              (alphayp*(pp(i,j+1,k)-pp(i,j  ,k))/dyc(j  ) - &
                               alphaym*(pp(i,j  ,k)-pp(i,j-1,k))/dyc(j-1))/dyf(j) + &
                              (alphazp*(pp(i,j,k+1)-pp(i,j,k  ))/dzc(k  ) - &
                               alphazm*(pp(i,j,k  )-pp(i,j,k-1))/dzc(k-1))/dzf(k) )
#else
          p(i,j,k) = p(i,j,k) + & 
                      alpha*( ((pp(i+1,j,k)-pp(i  ,j,k))/dxc(i  ) - &
                               (pp(i  ,j,k)-pp(i-1,j,k))/dxc(i-1))/dxf(i) + &
                              ((pp(i,j+1,k)-pp(i,j  ,k))/dyc(j  ) - &
                               (pp(i,j  ,k)-pp(i,j-1,k))/dyc(j-1))/dyf(j) + &
                              ((pp(i,j,k+1)-pp(i,j,k  ))/dzc(k  ) - &
                               (pp(i,j,k  )-pp(i,j,k-1))/dzc(k-1))/dzf(k) )
#endif
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
#endif
  end subroutine updt_pressure
end module mod_updt_pressure
