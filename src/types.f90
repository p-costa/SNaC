module mod_types
  use, intrinsic :: iso_fortran_env, only: stdin  => input_unit , &
                                           stdout => output_unit, &
                                           stderr => error_unit
  use, intrinsic :: iso_fortran_env, only: sp => real32, &
                                           dp => real64
  use            :: mpi_f08        , only: MPI_REAL_SP => MPI_REAL, &
                                           MPI_REAL_DP => MPI_DOUBLE_PRECISION
#ifdef _SINGLE_PRECISION
  use, intrinsic :: iso_fortran_env, only: rp          => real32
  use               mpi_f08        , only: MPI_REAL_RP => MPI_REAL
#else
  use, intrinsic :: iso_fortran_env, only: rp          => real64
  use               mpi_f08        , only: MPI_REAL_RP => MPI_DOUBLE_PRECISION
#endif
end module mod_types
