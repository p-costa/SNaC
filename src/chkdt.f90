module mod_chkdt
  use mpi
  use mod_common_mpi, only:ierr
  use mod_types
  implicit none
  private
  public chkdt
  contains
  subroutine chkdt(lo,hi,dxc,dyc,dzc,dxf,dyf,dzf,visc,u,v,w,dtmax)
    !
    ! computes maximum allowed timestep
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in ), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in ), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in ) :: visc
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(out) :: dtmax
    real(rp) :: ux,uy,uz,vx,vy,vz,wx,wy,wz
    real(rp) :: dtix,dtiy,dtiz,dti,dlmin
    integer :: i,j,k
    !
    dti = 0.
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,u,v,w,dxc,dyc,dzc,dxf,dyf,dzf) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP PRIVATE(ux,uy,uz,vx,vy,vz,wx,wy,wz,dtix,dtiy,dtiz) &
    !$OMP REDUCTION(max:dti)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          ux = abs(u(i,j,k))
          vx = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i+1,j,k)+v(i+1,j-1,k) )
          wx = 0.25_rp*abs( w(i,j,k)+w(i,j,k-1)+w(i+1,j,k)+w(i+1,j,k-1) )
          dtix = ux/dxc(i)+vx/dyf(j)+wx/dzf(k)
          uy = 0.25_rp*abs( u(i,j,k)+u(i,j+1,k)+u(i-1,j+1,k)+u(i-1,j,k) )
          vy = abs(v(i,j,k))
          wy = 0.25_rp*abs( w(i,j,k)+w(i,j+1,k)+w(i,j+1,k-1)+w(i,j,k-1) )
          dtiy = uy/dxf(i)+vy/dyc(j)+wy/dzf(k)
          uz = 0.25_rp*abs( u(i,j,k)+u(i-1,j,k)+u(i-1,j,k+1)+u(i,j,k+1) )
          vz = 0.25_rp*abs( v(i,j,k)+v(i,j-1,k)+v(i,j-1,k+1)+v(i,j,k+1) )
          wz = abs(w(i,j,k))
          dtiz = uz/dxf(i)+vz/dyf(j)+wz/dzc(k)
          dti = max(dti,dtix,dtiy,dtiz)
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,dti,1,MPI_REAL_RP,MPI_MAX,MPI_COMM_WORLD,ierr)
    if(dti.eq.0._rp) dti = 1._rp
    dlmin     = min(minval(dxf),minval(dyf),minval(dzf))
#ifdef IMPDIFF
    dtmax = sqrt(3.)/dti
#else
    dtmax = min(1.65_rp/12._rp/visc*dlmin**2,sqrt(3._rp)/dti)
#endif
    return
    !
  end subroutine chkdt
end module mod_chkdt
