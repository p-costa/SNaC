module mod_scal
  use mod_types
  implicit none
  private
  public scal_a,scal_d
  contains
  subroutine scal_d(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,s,dsdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in) :: alpha
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in) :: s
    real(rp), dimension(lo(1):  ,lo(2)  :,lo(3):  ), intent(inout) :: dsdt
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(static) DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,s,dsdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))/dxc(i  )
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))/dxc(i-1)
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))/dyc(j  )
          dsdym = (s(i,j  ,k)-s(i,j-1,k))/dyc(j-1)
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))/dzc(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))/dzc(k-1)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + &
                        alpha*(dsdxp-dsdxm)/dxf(i) + &
                        alpha*(dsdyp-dsdym)/dyf(j) + &
                        alpha*(dsdzp-dsdzm)/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine scal_d
  subroutine scal_a(lo,hi,dxf,dyf,dzf,u,v,w,s,dsdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u,v,w,s
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dsdt
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO COLLAPSE(1) SCHEDULE(static) DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP SHARED(lo,hi,dxf,dyf,dzf,u,v,w,s,dsdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          !
          usim  = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          dsdt(i,j,k) = dsdt(i,j,k) + &
                        ( -usip + usim )/dxf(i) + &
                        ( -vsjp + vsjm )/dyf(j) + &
                        ( -wskp + wskm )/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
  end subroutine scal_a
end module mod_scal
