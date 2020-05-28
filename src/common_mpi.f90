module mod_common_mpi
  use mpi
  implicit none
  integer :: myid,ierr
  integer, dimension(0:1,3) :: nb
  logical, dimension(0:1,3) :: is_bound
  integer, dimension(3) :: halos
end module mod_common_mpi
