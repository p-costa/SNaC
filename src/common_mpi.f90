module mod_common_mpi
  use mpi_f08, only: MPI_COMM
  implicit none
  integer :: myid,ierr
  type(MPI_COMM) :: comm_cart
end module mod_common_mpi
