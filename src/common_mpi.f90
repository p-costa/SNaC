module mod_common_mpi
  use mpi_f08, only: MPI_COMM
  implicit none
  integer :: myid,myid_block,ierr
  type(MPI_COMM) :: comm_block
end module mod_common_mpi
