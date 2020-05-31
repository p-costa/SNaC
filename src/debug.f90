module mod_debug
  use mpi
  use mod_types
  implicit none
  private
  public chkmean
  contains
  subroutine chkmean(lo,hi,l,dx,dy,dz,p,mean)
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
    real(rp), intent(out) :: mean
    integer :: i,j,k
    integer :: ierr
    mean = 0._rp
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,p,dx,dy,dz,l) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:mean)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          mean = mean + p(i,j,k)*dx(i)*dy(j)*dz(k)/(l(1)*l(2)*l(3))
        enddo
      enddo
    enddo
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    return
  end subroutine chkmean
end module mod_debug
