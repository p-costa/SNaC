module mod_debug
  use mpi_f08
  use mod_types
  implicit none
  private
  public chkmean,chk_helmholtz
  contains
  subroutine chkmean(lo,hi,l,dx,dy,dz,p,vol,comm,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ), dimension(3) :: l
    real(rp), intent(in ), dimension(lo(1)-1:) :: dx
    real(rp), intent(in ), dimension(lo(2)-1:) :: dy
    real(rp), intent(in ), dimension(lo(3)-1:) :: dz
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), intent(in ) :: vol
    type(MPI_COMM), intent(in ) :: comm
    real(rp), intent(out) :: mean
    integer :: i,j,k
    mean = 0._rp
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,p,dx,dy,dz,vol) &
    !$OMP REDUCTION(+:mean)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          mean = mean + p(i,j,k)*dx(i)*dy(j)*dz(k)/vol
        end do
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,comm)
  end subroutine chkmean
  subroutine chk_helmholtz(lo,hi,is_centered,dx1,dx2,dy1,dy2,dz1,dz2,alpha,fpp,fp,diffmax)
    !
    ! this subroutine checks if the implementation of implicit diffusion is
    ! correct
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    logical , intent(in ), dimension(3) :: is_centered
    real(rp), intent(in ), dimension(lo(1)-1:) :: dx1,dx2
    real(rp), intent(in ), dimension(lo(2)-1:) :: dy1,dy2
    real(rp), intent(in ), dimension(lo(3)-1:) :: dz1,dz2
    real(rp), intent(in )                      :: alpha
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: fpp,fp
    real(rp), intent(out) :: diffmax
    real(rp) :: val
    integer :: i,j,k
    integer, dimension(3) :: q
    q(:) = 0
    where(.not.is_centered(:)) q(:) = 1
    !
    diffmax = 0._rp
    do k=lo(3)+1,hi(3)-1
      do j=lo(2)+1,hi(2)-1
        do i=lo(1)+1,hi(1)-1
          val = alpha*fpp(i,j,k) + &
                 ((fpp(i+1,j,k)-fpp(i  ,j,k))/dx1(i  +q(1)) - &
                  (fpp(i  ,j,k)-fpp(i-1,j,k))/dx1(i-1+q(1)))/dx2(i) + &
                 ((fpp(i,j+1,k)-fpp(i,j  ,k))/dy1(j  +q(2)) - &
                  (fpp(i,j  ,k)-fpp(i,j-1,k))/dy1(j-1+q(2)))/dy2(j) + &
                 ((fpp(i,j,k+1)-fpp(i,j,k  ))/dz1(k  +q(3)) - &
                  (fpp(i,j,k  )-fpp(i,j,k-1))/dz1(k-1+q(3)))/dz2(k)
          diffmax = max(diffmax,abs(val-fp(i,j,k)))
        end do
      end do
    end do
    call mpi_allreduce(MPI_IN_PLACE,diffmax,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD)
  end subroutine chk_helmholtz
  subroutine mean_boundary_force(dt,factor,l,tau,tauo,f)
    implicit none
    type wallshear
      real(rp), dimension(3) :: x,y,z
    end type wallshear
    real(rp)       , intent(in   )               :: dt,factor(2),l(3)
    type(wallshear), intent(in   )               :: tau
    type(wallshear), intent(inout)               :: tauo
    real(rp)       , intent(inout), dimension(3) :: f
    f(1) = f(1) + dt*(factor(1)*sum(tau%x(:)/l(:)) + factor(2)*sum(tauo%x(:)/l(:)))
    f(2) = f(2) + dt*(factor(1)*sum(tau%y(:)/l(:)) + factor(2)*sum(tauo%y(:)/l(:)))
    f(3) = f(3) + dt*(factor(1)*sum(tau%z(:)/l(:)) + factor(2)*sum(tauo%z(:)/l(:)))
    tauo = tau
  end subroutine mean_boundary_force
  subroutine compute_mean_wall_shear(is_bound,lo,hi,dxc,dyc,dzc,u,v,w,l,visc,tau)
    type wallshear
      real(rp), dimension(3) :: x,y,z
    end type wallshear
    logical        , intent(in ), dimension(0:1,3) :: is_bound
    integer        , intent(in ), dimension(3    ) :: lo,hi
    real(rp)       , intent(in ), dimension(lo(1)-1:) :: dxc
    real(rp)       , intent(in ), dimension(lo(2)-1:) :: dyc
    real(rp)       , intent(in ), dimension(lo(3)-1:) :: dzc
    real(rp)       , intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp)       , intent(in ) :: l(3),visc
    type(wallshear), intent(out) :: tau
    integer :: i,j,k,idir
    tau%x(:) = 0._rp
    idir = 2
    if(is_bound(0,idir)) then
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          tau%x(idir) = tau%x(idir) + &
            (u(i,lo(idir),k)-u(i,lo(idir)-1,k))/dyc(lo(idir)-1)*visc*dxc(i)*dzc(k)/(l(1)*l(3))
        end do
      end do
    end if
    if(is_bound(1,idir)) then
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          tau%x(idir) = tau%x(idir) + &
            (u(i,hi(idir),k)-u(i,hi(idir)+1,k))/dyc(hi(idir)+0)*visc*dxc(i)*dzc(k)/(l(1)*l(3))
        end do
      end do
    end if
    idir = 3
    if(is_bound(0,idir)) then
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          tau%x(idir) = tau%x(idir) + &
            (u(i,j,lo(idir))-u(i,j,lo(idir)-1))/dzc(lo(idir)-1)*visc*dxc(i)*dyc(j)/(l(1)*l(2))
        end do
      end do
    end if
    if(is_bound(1,idir)) then
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          tau%x(idir) = tau%x(idir) + &
            (u(i,j,hi(idir))-u(i,j,hi(idir)+1))/dzc(hi(idir)+0)*visc*dxc(i)*dyc(j)/(l(1)*l(2))
        end do
      end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,tau%x(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD)
    tau%y(:) = 0._rp
    idir = 1
    if(is_bound(0,idir)) then
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          tau%y(idir) = tau%y(idir) + &
            (v(lo(idir),j,k)-v(lo(idir)-1,j,k))/dxc(lo(idir)-1)*visc*dyc(j)*dzc(k)/(l(2)*l(3))
        end do
      end do
    end if
    if(is_bound(1,idir)) then
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          tau%y(idir) = tau%y(idir) + &
            (v(hi(idir),j,k)-v(hi(idir)+1,j,k))/dxc(hi(idir)+0)*visc*dyc(j)*dzc(k)/(l(2)*l(3))
        end do
      end do
    end if
    idir = 3
    if(is_bound(0,idir)) then
      do i=lo(1),hi(1)
        do j=lo(2),hi(2)
          tau%y(idir) = tau%y(idir) + &
            (v(i,j,lo(idir))-v(i,j,lo(idir)-1))/dzc(lo(idir)-1)*visc*dyc(j)*dxc(i)/(l(2)*l(1))
        end do
      end do
    end if
    if(is_bound(1,idir)) then
      do i=lo(1),hi(1)
        do j=lo(2),hi(2)
          tau%y(idir) = tau%y(idir) + &
            (v(i,j,hi(idir))-v(i,j,hi(idir)+1))/dzc(hi(idir)+0)*visc*dyc(j)*dxc(i)/(l(2)*l(1))
        end do
      end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,tau%y(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD)
    tau%z(:) = 0._rp
    idir = 1
    if(is_bound(0,idir)) then
      do j=lo(2),hi(2)
        do k=lo(3),hi(3)
          tau%z(idir) = tau%z(idir) + &
            (w(lo(idir),j,k)-w(lo(idir)-1,j,k))/dxc(lo(idir)-1)*visc*dzc(k)*dyc(j)/(l(3)*l(2))
        end do
      end do
    end if
    if(is_bound(1,idir)) then
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          tau%z(idir) = tau%z(idir) + &
            (w(hi(idir),j,k)-w(hi(idir)+1,j,k))/dxc(hi(idir)+0)*visc*dzc(k)*dyc(j)/(l(3)*l(2))
        end do
      end do
    end if
    idir = 2
    if(is_bound(0,idir)) then
      do i=lo(1),hi(1)
        do k=lo(3),hi(3)
          tau%z(idir) = tau%z(idir) + &
            (w(i,lo(idir),k)-w(i,lo(idir)-1,k))/dyc(lo(idir)-1)*visc*dzc(k)*dxc(i)/(l(3)*l(1))
        end do
      end do
    end if
    if(is_bound(1,idir)) then
      do i=lo(1),hi(1)
        do k=lo(3),hi(3)
          tau%z(idir) = tau%z(idir) + &
            (w(i,hi(idir),k)-w(i,hi(idir)+1,k))/dyc(hi(idir)+0)*visc*dzc(k)*dxc(i)/(l(3)*l(1))
        end do
      end do
    end if
    call mpi_allreduce(MPI_IN_PLACE,tau%z(1),3,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD)
  end subroutine compute_mean_wall_shear
end module mod_debug
