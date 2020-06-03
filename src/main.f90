!    SSSSSSSSSSSSSSS NNNNNNNN        NNNNNNNN                         CCCCCCCCCCCCC
!  SS:::::::::::::::SN:::::::N       N::::::N                      CCC::::::::::::C
! S:::::SSSSSS::::::SN::::::::N      N::::::N                    CC:::::::::::::::C
! S:::::S     SSSSSSSN:::::::::N     N::::::N                   C:::::CCCCCCCC::::C
! S:::::S            N::::::::::N    N::::::N  aaaaaaaaaaaaa   C:::::C       CCCCCC
! S:::::S            N:::::::::::N   N::::::N  a::::::::::::a C:::::C              
!  S::::SSSS         N:::::::N::::N  N::::::N  aaaaaaaaa:::::aC:::::C              
!   SS::::::SSSSS    N::::::N N::::N N::::::N           a::::aC:::::C              
!     SSS::::::::SS  N::::::N  N::::N:::::::N    aaaaaaa:::::aC:::::C              
!        SSSSSS::::S N::::::N   N:::::::::::N  aa::::::::::::aC:::::C              
!             S:::::SN::::::N    N::::::::::N a::::aaaa::::::aC:::::C              
!             S:::::SN::::::N     N:::::::::Na::::a    a:::::a C:::::C       CCCCCC
! SSSSSSS     S:::::SN::::::N      N::::::::Na::::a    a:::::a  C:::::CCCCCCCC::::C
! S::::::SSSSSS:::::SN::::::N       N:::::::Na:::::aaaa::::::a   CC:::::::::::::::C
! S:::::::::::::::SS N::::::N        N::::::N a::::::::::aa:::a    CCC::::::::::::C
!  SSSSSSSSSSSSSSS   NNNNNNNN         NNNNNNN  aaaaaaaaaa  aaaa       CCCCCCCCCCCCC
!-----------------------------------------------------------------------------------
! Slow CaNS, a.k.a. **SNaC**
! Pedro Costa (p.simoes.costa@gmail.com)
!-----------------------------------------------------------------------------------
program snac
  use mpi
  use mod_common_mpi, only: myid,ierr
  use mod_initmpi   , only: initmpi
  use mod_fillps    , only: fillps
  use mod_load      , only: load
  use mod_param     , only: read_input, &
                            ng,l,gt,gr,cfl,dtmin,uref,lref,rey,                  &
                            inivel,is_wallturb,nstep,time_max,tw_max,stop_type, &
                            restart,is_overwrite_save,                           &
                            icheck,iout0d,iout1d,iout2d,iout3d,isave,            &
                            cbcvel,bcvel,cbcpre,bcpre,                           &
                            bforce, is_forced,velf,is_outflow,                   &
                            dims,nthreadsmax
  use mod_sanity , only: test_sanity
  use mod_types
  implicit none
  integer , dimension(3) :: lo,hi
  logical :: kill
  !
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
  call read_input()
  !
  ! read parameter file
  !
  call read_input
  !
  ! initialize MPI/OpenMP
  !
  !$call omp_set_num_threads(nthreadsmax)
  call initmpi(ng,dims,cbcpre,lo,hi)
  !
  call MPI_FINALIZE(ierr)
  call exit
end program snac
