module mod_types
  use, intrinsic :: iso_fortran_env
  use, intrinsic :: iso_fortran_env, sp => real32
  use, intrinsic :: iso_fortran_env, dp => real64
#ifdef _SINGLE_PRECISION
  use, intrinsic :: iso_fortran_env, rp          => real32
  use            :: mpi            , MPI_REAL_RP => MPI_REAL
#else
  use, intrinsic :: iso_fortran_env, rp          => real64
  use            :: mpi            , MPI_REAL_RP => MPI_DOUBLE_PRECISION
#endif
end module mod_types
