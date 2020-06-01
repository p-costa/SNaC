module mod_mom
  use mpi
  use mod_types
  implicit none
  private
  public momx_a,momy_a,momz_a,momx_d,momy_d,momz_d,momx_p,momy_p,momz_p
  contains
  subroutine momx_a(lo,hi,dxf,dyf,dzf,u,v,w,dudt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u,v,w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dudt
    real(rp) :: uuip,uuim,uvjp,uvjm,uwkp,uwkm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uuip,uuim,uvjp,uvjm,uwkp,uwkm) &
    !$OMP SHARED(lo,hi,dxf,dyf,dzf,u,v,w,dudt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          uuim  = 0.25_rp*( u(i+1,j,k)+u(i,j,k) )*( u(i+1,j  ,k  )+u(i,j  ,k  ) )
          uuip  = 0.25_rp*( u(i-1,j,k)+u(i,j,k) )*( u(i-1,j  ,k  )+u(i,j  ,k  ) )
          uvjm  = 0.25_rp*( u(i,j+1,k)+u(i,j,k) )*( v(i+1,j  ,k  )+v(i,j  ,k  ) )
          uvjp  = 0.25_rp*( u(i,j-1,k)+u(i,j,k) )*( v(i+1,j-1,k  )+v(i,j-1,k  ) )
          uwkm  = 0.25_rp*( u(i,j,k+1)+u(i,j,k) )*( w(i+1,j  ,k  )+w(i,j  ,k  ) )
          uwkp  = 0.25_rp*( u(i,j,k-1)+u(i,j,k) )*( w(i+1,j  ,k-1)+w(i,j  ,k-1) )
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        ( -uuip + uuim )/dxf(i) + &
                        ( -uvjp + uvjm )/dyf(j) + &
                        ( -uwkp + uwkm )/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momx_a
  !
  subroutine momy_a(lo,hi,dxf,dyf,dzf,u,v,w,dvdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u,v,w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dvdt
    real(rp) :: uvip,uvim,vvjp,vvjm,wvkp,wvkm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uvip,uvim,vvjp,vvjm,wvkp,wvkm) &
    !$OMP SHARED(lo,hi,dxf,dyf,dzf,u,v,w,dvdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          uvip  = 0.25_rp*( u(i  ,j,k  )+u(i  ,j+1,k  ) )*( v(i,j,k  )+v(i+1,j  ,k) )
          uvim  = 0.25_rp*( u(i-1,j,k  )+u(i-1,j+1,k  ) )*( v(i,j,k  )+v(i-1,j  ,k) )
          vvjp  = 0.25_rp*( v(i  ,j,k  )+v(i  ,j+1,k  ) )*( v(i,j,k  )+v(i  ,j+1,k) )
          vvjm  = 0.25_rp*( v(i  ,j,k  )+v(i  ,j-1,k  ) )*( v(i,j,k  )+v(i  ,j-1,k) )
          wvkp  = 0.25_rp*( w(i  ,j,k  )+w(i  ,j+1,k  ) )*( v(i,j,k+1)+v(i  ,j  ,k) )
          wvkm  = 0.25_rp*( w(i  ,j,k-1)+w(i  ,j+1,k-1) )*( v(i,j,k-1)+v(i  ,j  ,k) )
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        ( -uvip + uvim )/dxf(i) + &
                        ( -vvjp + vvjm )/dyf(j) + &
                        ( -wvkp + wvkm )/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momy_a
  !
  subroutine momz_a(lo,hi,dxf,dyf,dzf,u,v,w,dwdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u,v,w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dwdt
    real(rp) :: uwip,uwim,vwjp,vwjm,wwkp,wwkm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(uwip,uwim,vwjp,vwjm,wwkp,wwkm) &
    !$OMP SHARED(lo,hi,dxf,dyf,dzf,u,v,w,dwdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          uwip  = 0.25_rp*( w(i,j,k)+w(i+1,j,k) )*( u(i  ,j  ,k)+u(i  ,j  ,k+1) )
          uwim  = 0.25_rp*( w(i,j,k)+w(i-1,j,k) )*( u(i-1,j  ,k)+u(i-1,j  ,k+1) )
          vwjp  = 0.25_rp*( w(i,j,k)+w(i,j+1,k) )*( v(i  ,j  ,k)+v(i  ,j  ,k+1) )
          vwjm  = 0.25_rp*( w(i,j,k)+w(i,j-1,k) )*( v(i  ,j-1,k)+v(i  ,j-1,k+1) )
          wwkp  = 0.25_rp*( w(i,j,k)+w(i,j,k+1) )*( w(i  ,j  ,k)+w(i  ,j  ,k+1) )
          wwkm  = 0.25_rp*( w(i,j,k)+w(i,j,k-1) )*( w(i  ,j  ,k)+w(i  ,j  ,k-1) )
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        ( -uwip + uwim )/dxf(i) + &
                        ( -vwjp + vwjm )/dyf(j) + &
                        ( -wwkp + wwkm )/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momz_a
  !
  subroutine momx_d(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,u,dudt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in)                      :: visc
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP SHARED(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,u,dudt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dudxp = (u(i+1,j,k)-u(i,j,k))/dxf(i+1)
          dudxm = (u(i,j,k)-u(i-1,j,k))/dxf(i  )
          dudyp = (u(i,j+1,k)-u(i,j,k))/dyc(j  )
          dudym = (u(i,j,k)-u(i,j-1,k))/dyc(j-1)
          dudzp = (u(i,j,k+1)-u(i,j,k))/dzc(k  )
          dudzm = (u(i,j,k)-u(i,j,k-1))/dzc(k-1)
          !
          ! Momentum balance
          !
          dudt(i,j,k) = dudt(i,j,k) + &
                        visc*(dudxp-dudxm)/dxc(i) + &
                        visc*(dudyp-dudym)/dyf(j) + &
                        visc*(dudzp-dudzm)/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momx_d
  !
  subroutine momy_d(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,v,dvdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in)                      :: visc
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: v
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP SHARED(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,v,dvdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dvdxp = (v(i+1,j,k)-v(i,j,k))/dxc(i  )
          dvdxm = (v(i,j,k)-v(i-1,j,k))/dxc(i-1)
          dvdyp = (v(i,j+1,k)-v(i,j,k))/dyf(j+1)
          dvdym = (v(i,j,k)-v(i,j-1,k))/dyf(j  )
          dvdzp = (v(i,j,k+1)-v(i,j,k))/dzc(k  )
          dvdzm = (v(i,j,k)-v(i,j,k-1))/dzc(k-1)
          !
          ! Momentum balance
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        visc*(dvdxp-dvdxm)/dxf(i) + &
                        visc*(dvdyp-dvdym)/dyc(j) + &
                        visc*(dvdzp-dvdzm)/dzf(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momy_d
  !
  subroutine momz_d(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,w,dwdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in)                      :: visc
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dwdt
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP SHARED(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,w,dwdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dwdxp = (w(i+1,j,k)-w(i,j,k))*dxc(i  )
          dwdxm = (w(i,j,k)-w(i-1,j,k))*dxc(i-1)
          dwdyp = (w(i,j+1,k)-w(i,j,k))*dyc(j  )
          dwdym = (w(i,j,k)-w(i,j-1,k))*dyc(j-1)
          dwdzp = (w(i,j,k+1)-w(i,j,k))*dzf(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))*dzf(k  )
          !
          ! Momentum balance
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        visc*(dwdxp-dwdxm)/dxf(i) + &
                        visc*(dwdyp-dwdym)/dyf(j) + &
                        visc*(dwdzp-dwdzm)/dzc(k)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momz_d
  subroutine momx_p(lo,hi,dxc,bforce,p,dudt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc
    real(rp), intent(in) :: bforce
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in ) :: p
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(out) :: dudt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,dxc,bforce,p,dudt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dudt(i,j,k) = - ( p(i-1,j,k)-p(i,j,k) )/dxc(i) + bforce
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momx_p
  !
  subroutine momy_p(lo,hi,dyc,bforce,p,dvdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc
    real(rp), intent(in) :: bforce
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in ) :: p
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(out) :: dvdt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,dyc,bforce,p,dvdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dvdt(i,j,k) = - ( p(i,j-1,k)-p(i,j,k) )/dyc(j) + bforce
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momy_p
  !
  subroutine momz_p(lo,hi,dzc,bforce,p,dwdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc
    real(rp), intent(in) :: bforce
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in ) :: p
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(out) :: dwdt
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,dzc,bforce,p,dwdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dwdt(i,j,k) = - ( p(i,j,k-1)-p(i,j,k) )/dzc(k) + bforce
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    return
  end subroutine momz_p
end module mod_mom
