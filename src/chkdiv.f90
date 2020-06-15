module mod_chkdiv
  use mpi
  use mod_common_mpi, only: ierr
  use mod_types
  implicit none
  private
  public chkdiv
  contains
    subroutine chkdiv(lo,hi,dxf,dyf,dzf,u,v,w,vol,mpi_comm,divtot,divmax)
    !
    ! checks the divergence of the velocity field
    !
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in ), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in ), dimension(lo(3)-1:) :: dzf
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    real(rp), intent(in ) :: vol 
    integer , intent(in ) :: mpi_comm
    real(rp), intent(out) :: divtot,divmax
    real(rp) :: div
    integer  :: i,j,k
    !
    divtot = 0._rp
    divmax = 0._rp
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,u,v,w,dxf,dyf,dzf,vol) &
    !$OMP PRIVATE(i,j,k,div) &
    !$OMP REDUCTION(+:divtot) &
    !$OMP REDUCTION(max:divmax)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          div = (w(i,j,k)-w(i,j,k-1))/dzf(k) + &
                (v(i,j,k)-v(i,j-1,k))/dyf(j) + &
                (u(i,j,k)-u(i-1,j,k))/dxf(i)
          divmax = max(divmax,abs(div))
          divtot = divtot + div*dxf(i)*dyf(j)*dzf(k)/vol
          if(abs(div).gt.0.001) print*,div,i,j,k
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,divtot,1,MPI_REAL_RP,MPI_SUM,mpi_comm,ierr)
    call mpi_allreduce(MPI_IN_PLACE,divmax,1,MPI_REAL_RP,MPI_MAX,mpi_comm,ierr)
    return
  end subroutine chkdiv
end module mod_chkdiv
