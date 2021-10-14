module mod_types
  use, intrinsic :: iso_fortran_env, only: stdin  => input_unit , &
                                           stdout => output_unit, &
                                           stderr => error_unit
  use            :: mpi_f08        , only: MPI_REAL_SP => MPI_REAL, &
                                           MPI_REAL_DP => MPI_DOUBLE_PRECISION
#ifdef _SINGLE_PRECISION
  use               mpi_f08        , only: MPI_REAL_RP => MPI_REAL
#else
  use               mpi_f08        , only: MPI_REAL_RP => MPI_DOUBLE_PRECISION
#endif
  integer, parameter :: sp = selected_real_kind(6 , 37), &
                        dp = selected_real_kind(15,307), &
                        i8 = selected_int_kind(18)
#ifdef _SINGLE_PRECISION
  integer, parameter :: rp = sp
#else
  integer, parameter :: rp = dp
#endif
end module mod_types
