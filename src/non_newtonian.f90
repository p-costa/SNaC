module mod_non_newtonian
  use mpi_f08
  use mod_types
  implicit none
  private
  public strain_rate_norm,compute_viscosity,momx_d_nn,momy_d_nn,momz_d_nn
  contains
  subroutine momx_d_nn(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,mu,u,v,w,dudt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: mu,u,v,w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dudt
    real(rp) :: dudxp,dudxm,dudyp,dudym,dudzp,dudzm, &
                dvdxp,dvdxm,dwdxp,dwdxm, &
                muxp,muxm,muyp,muym,muzp,muzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dudxp,dudxm,dudyp,dudym,dudzp,dudzm) &
    !$OMP PRIVATE(dvdxp,dvdxm,dwdxp,dwdxm) &
    !$OMP PRIVATE(muxp,muxm,muyp,muym,muzp,muzm) &
    !$OMP SHARED(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,mu,u,v,w,dudt)
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
          dvdxp = (v(i+1,j  ,k  )-v(i  ,j  ,k  ))/dxc(i)
          dvdxm = (v(i+1,j-1,k  )-v(i  ,j-1,k  ))/dxc(i)
          dwdxp = (w(i+1,j  ,k  )-w(i  ,j  ,k  ))/dxc(i)
          dwdxm = (w(i+1,j  ,k-1)-w(i  ,j  ,k-1))/dxc(i)
          !
          muxp = mu(i+1,j,k)
          muxm = mu(i  ,j,k)
          muyp = 0.25_rp*(mu(i,j,k)+mu(i,j+1,k)+mu(i+1,j+1,k)+mu(i+1,j,k))
          muym = 0.25_rp*(mu(i,j,k)+mu(i,j-1,k)+mu(i+1,j-1,k)+mu(i+1,j,k))
          muzp = 0.25_rp*(mu(i,j,k)+mu(i,j,k+1)+mu(i+1,j,k+1)+mu(i+1,j,k))
          muzm = 0.25_rp*(mu(i,j,k)+mu(i,j,k-1)+mu(i+1,j,k-1)+mu(i+1,j,k))
          dudt(i,j,k) = dudt(i,j,k) + &
                        ((dudxp+dudxp)*muxp-(dudxm+dudxm)*muxm)/dxc(i) + &
                        ((dudyp+dvdxp)*muyp-(dudym+dvdxm)*muym)/dyf(j) + &
                        ((dudzp+dwdxp)*muzp-(dudzm+dwdxm)*muzm)/dzf(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momx_d_nn
  !
  subroutine momy_d_nn(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,mu,u,v,w,dvdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: mu,u,v,w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dvdt
    real(rp) :: dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm, &
                dudyp,dudym,dwdyp,dwdym, &
                muxp,muxm,muyp,muym,muzp,muzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dvdxp,dvdxm,dvdyp,dvdym,dvdzp,dvdzm) &
    !$OMP PRIVATE(dudyp,dudym,dwdyp,dwdym) &
    !$OMP PRIVATE(muxp,muxm,muyp,muym,muzp,muzm) &
    !$OMP SHARED(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,mu,u,v,w,dvdt)
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
          dudyp = (u(i  ,j+1,k  )-u(i  ,j  ,k  ))/dyc(j)
          dudym = (u(i-1,j+1,k  )-u(i-1,j  ,k  ))/dyc(j)
          dwdyp = (w(i  ,j+1,k  )-w(i  ,j  ,k  ))/dyc(j)
          dwdym = (w(i  ,j+1,k-1)-w(i  ,j  ,k-1))/dyc(j)
          !
          muxp = 0.25_rp*(mu(i,j,k)+mu(i+1,j,k)+mu(i+1,j+1,k)+mu(i,j+1,k))
          muxm = 0.25_rp*(mu(i,j,k)+mu(i-1,j,k)+mu(i-1,j+1,k)+mu(i,j+1,k))
          muyp = mu(i,j+1,k)
          muym = mu(i,j  ,k)
          muzp = 0.25_rp*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k+1)+mu(i,j,k+1))
          muzm = 0.25_rp*(mu(i,j,k)+mu(i,j+1,k)+mu(i,j+1,k-1)+mu(i,j,k-1))
          !
          dvdt(i,j,k) = dvdt(i,j,k) + &
                        ((dvdxp+dudyp)*muxp-(dvdxm+dudym)*muxm)/dxf(i) + &
                        ((dvdyp+dvdyp)*muyp-(dvdym+dvdym)*muym)/dyc(j) + &
                        ((dvdzp+dwdyp)*muzp-(dvdzm+dwdym)*muzm)/dzf(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momy_d_nn
  !
  subroutine momz_d_nn(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,mu,u,v,w,dwdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: mu,u,v,w
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dwdt
    real(rp) :: dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm, &
                dudzp,dudzm,dvdzp,dvdzm, &
                muxp,muxm,muyp,muym,muzp,muzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dwdxp,dwdxm,dwdyp,dwdym,dwdzp,dwdzm) &
    !$OMP PRIVATE(dudzp,dudzm,dvdzp,dvdzm) &
    !$OMP PRIVATE(muxp,muxm,muyp,muym,muzp,muzm) &
    !$OMP SHARED(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,mu,u,v,w,dwdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dwdxp = (w(i+1,j,k)-w(i,j,k))/dxc(i  )
          dwdxm = (w(i,j,k)-w(i-1,j,k))/dxc(i-1)
          dwdyp = (w(i,j+1,k)-w(i,j,k))/dyc(j  )
          dwdym = (w(i,j,k)-w(i,j-1,k))/dyc(j-1)
          dwdzp = (w(i,j,k+1)-w(i,j,k))/dzf(k+1)
          dwdzm = (w(i,j,k)-w(i,j,k-1))/dzf(k  )
          !
          dudzp = (u(i  ,j  ,k+1)-u(i  ,j  ,k  ))/dzc(k)
          dudzm = (u(i-1,j  ,k+1)-u(i-1,j  ,k  ))/dzc(k)
          dvdzp = (v(i  ,j  ,k+1)-v(i  ,j  ,k  ))/dzc(k)
          dvdzm = (v(i  ,j-1,k+1)-v(i  ,j-1,k  ))/dzc(k)
          !
          muxp = 0.25_rp*(mu(i,j,k)+mu(i,j,k+1)+mu(i+1,j  ,k+1)+mu(i+1,j  ,k) )
          muxm = 0.25_rp*(mu(i,j,k)+mu(i,j,k+1)+mu(i-1,j  ,k+1)+mu(i-1,j  ,k) )
          muyp = 0.25_rp*(mu(i,j,k)+mu(i,j,k+1)+mu(i  ,j+1,k+1)+mu(i  ,j+1,k) )
          muym = 0.25_rp*(mu(i,j,k)+mu(i,j,k+1)+mu(i  ,j-1,k+1)+mu(i  ,j-1,k) )
          muzp = mu(i,j,k+1)
          muzm = mu(i,j,k  )
          !
          dwdt(i,j,k) = dwdt(i,j,k) + &
                        ((dwdxp+dudzp)*muxp-(dwdxm+dudzm)*muxm)/dxf(i) + &
                        ((dwdyp+dvdzp)*muyp-(dwdym+dvdzm)*muym)/dyf(j) + &
                        ((dwdzp+dwdzp)*muzp-(dwdzm+dwdzm)*muzm)/dzc(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine momz_d_nn
  subroutine strain_rate_norm(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,gamma_dot)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u,v,w
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(inout) :: gamma_dot
    real(rp) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz) &
    !$OMP SHARED(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,gamma_dot)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          !
          dudx = (u(i+1,j,k)-u(i,j,k))/dxf(i)
          dvdy = (v(i,j+1,k)-v(i,j,k))/dyf(j)
          dwdz = (w(i,j,k+1)-w(i,j,k))/dzf(k)
          !
          dudy = 0.25_rp*((u(i  ,j+1,k)-u(i  ,j  ,k))/dyc(j  ) + &
                          (u(i  ,j  ,k)-u(i  ,j-1,k))/dyc(j-1) + &
                          (u(i-1,j+1,k)-u(i-1,j  ,k))/dyc(j  ) + &
                          (u(i-1,j  ,k)-u(i-1,j-1,k))/dyc(j-1) )
          dudz = 0.25_rp*((u(i  ,j,k+1)-u(i  ,j,k  ))/dzc(k  ) + &
                          (u(i  ,j,k  )-u(i  ,j,k-1))/dzc(k-1) + &
                          (u(i-1,j,k+1)-u(i-1,j,k  ))/dzc(k  ) + &
                          (u(i-1,j,k  )-u(i-1,j,k-1))/dzc(k-1) )
          dvdx = 0.25_rp*((v(i+1,j  ,k)-v(i  ,j  ,k))/dxc(i  ) + &
                          (v(i  ,j  ,k)-v(i-1,j  ,k))/dxc(i-1) + &
                          (v(i+1,j-1,k)-v(i  ,j-1,k))/dxc(i  ) + &
                          (v(i  ,j-1,k)-v(i-1,j-1,k))/dxc(i-1) )
          dvdz = 0.25_rp*((v(i,j  ,k+1)-v(i,j  ,k  ))/dzc(k  ) + &
                          (v(i,j  ,k  )-v(i,j  ,k-1))/dzc(k-1) + &
                          (v(i,j-1,k+1)-v(i,j-1,k  ))/dzc(k  ) + &
                          (v(i,j-1,k  )-v(i,j-1,k-1))/dzc(k-1) )
          dwdx = 0.25_rp*((w(i+1,j,k  )-w(i  ,j,k  ))/dxc(i  ) + &
                          (w(i  ,j,k  )-w(i-1,j,k  ))/dxc(i-1) + &
                          (w(i+1,j,k-1)-w(i  ,j,k-1))/dxc(i  ) + &
                          (w(i  ,j,k-1)-w(i-1,j,k-1))/dxc(i-1) )
          dwdy = 0.25_rp*((w(i,j+1,k  )-w(i,j  ,k  ))/dyc(j  ) + &
                          (w(i,j  ,k  )-w(i,j-1,k  ))/dyc(j-1) + &
                          (w(i,j+1,k-1)-w(i,j  ,k-1))/dyc(j  ) + &
                          (w(i,j  ,k-1)-w(i,j-1,k-1))/dyc(j-1) )
          gamma_dot(i,j,k) = sqrt( 2._rp*(dudx**2+dvdy**2+dwdz**2) + &
                                   (dudy+dvdx)**2+(dudz+dwdx)**2+(dvdz+dwdy)**2 )
        end do
      end do
    end do
  end subroutine strain_rate_norm
  subroutine compute_viscosity(lo,hi,kappa,rn,tau0,eps,gamma_dot)
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in) :: kappa,rn,tau0,eps
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(inout) :: gamma_dot
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP SHARED(lo,hi,kappa,rn,tau0,eps,gamma_dot)
    do k=lo(3)-1,hi(3)+1
      do j=lo(2)-1,hi(2)+1
        do i=lo(1)-1,hi(1)+1
          gamma_dot(i,j,k) = herschel_bulkley_visc(kappa,rn,tau0,eps,gamma_dot(i,j,k))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
  end subroutine
  pure real(rp) function herschel_bulkley_visc(kappa,rn,tau0,eps,gamma_dot) result(visc)
    real(rp), parameter :: small = epsilon(0._rp)*10._rp**(precision(0._rp)/2._rp)
    real(rp), intent(in) :: kappa,rn,tau0,eps,gamma_dot
    real(rp) :: gamma_dot_aux
    gamma_dot_aux = gamma_dot + small
    visc = kappa*gamma_dot_aux**(rn-1._rp) + tau0*(1._rp-exp(-gamma_dot_aux/eps))/gamma_dot_aux
  end function herschel_bulkley_visc
end module mod_non_newtonian
