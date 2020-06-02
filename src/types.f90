module mod_types
  use, intrinsic :: iso_fortran_env, stdin  => input_unit , &
                                     stdout => output_unit, &
                                     stderr => error_unit
  use, intrinsic :: iso_fortran_env, sp => real32, &
                                     dp => real64
  use            :: mpi            , MPI_REAL_DP => MPI_REAL, &
                                     MPI_REAL_SP => MPI_DOUBLE_PRECISION
#ifdef _SINGLE_PRECISION
  use, intrinsic :: iso_fortran_env, rp          => real32
  use            :: mpi            , MPI_REAL_RP => MPI_REAL
#else
  use, intrinsic :: iso_fortran_env, rp          => real64
  use            :: mpi            , MPI_REAL_RP => MPI_DOUBLE_PRECISION
#endif
end module mod_types
