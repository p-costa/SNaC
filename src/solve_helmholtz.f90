module mod_solve_helmholtz
#ifdef _IMPDIFF
  use, intrinsic :: iso_c_binding, only: C_PTR
#ifdef _FFT_USE_SLABS
  use mpi_f08
#endif
  use mod_bound     , only: updt_rhs
  use mod_solver    , only: hypre_solver, &
                            add_weighted_constant_to_diagonal,create_solver,setup_solver, &
                            solve_helmholtz,finalize_solver
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  use mod_fft       , only: fft
  use mod_solver    , only: add_weighted_constant_to_n_diagonals,create_n_solvers,setup_n_solvers, &
                            solve_n_helmholtz_2d,finalize_n_solvers
#ifdef _FFT_USE_SLABS
  use mod_solver    , only: alltoallw,transpose_slab
#endif
#endif
  use mod_types
  implicit none
  private
  public solve_impdiff_field
  interface solve_impdiff_field
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
    procedure solve_impdiff_field_fft
#else
    procedure solve_impdiff_field_3d
#endif
  end interface
  contains
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
  subroutine solve_impdiff_field_fft(alphai_new,alphai_old,lo,hi,lo_a,hif,hif_a, &
                                     idir,il,iu,iskip,is_bound,rhsx,rhsy,rhsz,dl1_2,dl2_2, &
                                     hypre_maxiter,hypre_tol,hypre_solver_i, &
#ifndef _FFT_USE_SLABS
                                     phi,phio,arrplan,normfft,asolver)
#else
                                     phi,phio,arrplan,normfft,asolver,t_params,comm_block,phi_s)
#endif
    !
    ! solves an implicit-diffusion Helmholtz problem for one field
    !
    implicit none
    real(rp), intent(in) :: alphai_new,alphai_old
    integer , intent(in), dimension(3) :: lo,hi,lo_a,hif,hif_a
    integer , intent(in) :: idir,il,iu,iskip
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(lo(2):,lo(3):,0:) :: rhsx
    real(rp), intent(in), dimension(lo(1):,lo(3):,0:) :: rhsy
    real(rp), intent(in), dimension(lo(1):,lo(2):,0:) :: rhsz
    real(rp), intent(in), dimension(lo(il)-1:) :: dl1_2
    real(rp), intent(in), dimension(lo(iu)-1:) :: dl2_2
    integer , intent(in) :: hypre_maxiter,hypre_solver_i
    real(rp), intent(in) :: hypre_tol
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: phi
    real(rp), intent(inout), dimension(:,:,:) :: phio
    type(C_PTR), intent(in), dimension(2) :: arrplan
    real(rp), intent(in) :: normfft
    type(hypre_solver), intent(inout), dimension(:) :: asolver
#ifdef _FFT_USE_SLABS
    type(alltoallw), intent(in), dimension(:,:) :: t_params
    type(MPI_COMM), intent(in) :: comm_block
    real(rp), intent(inout), dimension(:,:,:) :: phi_s
#endif
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k) COLLAPSE(3)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
#ifdef _FFT_X
          phi(i,j,k) = phi(i,j,k)*alphai_new*dl1_2(j)*dl2_2(k)
#elif  _FFT_Y
          phi(i,j,k) = phi(i,j,k)*alphai_new*dl1_2(i)*dl2_2(k)
#elif  _FFT_Z
          phi(i,j,k) = phi(i,j,k)*alphai_new*dl1_2(i)*dl2_2(j)
#endif
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call updt_rhs(lo,hif,is_bound,rhsx,rhsy,rhsz,phi)
    call add_weighted_constant_to_n_diagonals(hif_a(idir)-lo_a(idir)+1,lo_a(il:iu:iskip),hif_a(il:iu:iskip), &
                                              dl1_2,dl2_2,alphai_new-alphai_old,asolver(:)%mat)
    call create_n_solvers(hif_a(idir)-lo_a(idir)+1,hypre_maxiter,hypre_tol,.true.,hypre_solver_i,asolver)
    call setup_n_solvers(hif_a(idir)-lo_a(idir)+1,asolver)
    call fft(arrplan(1),phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
#ifndef _FFT_USE_SLABS
    call solve_n_helmholtz_2d(asolver,lo(idir),hif(idir),1,lo(il:iu:iskip),hif(il:iu:iskip),phi,phio)
#else
    call transpose_slab(1,0,t_params(:,1:2:1 ),comm_block,phi,phi_s)
    call solve_n_helmholtz_2d(asolver,lo_a(idir),hif_a(idir),0,lo_a(il:iu:iskip),hif_a(il:iu:iskip),phi_s,phio)
    call transpose_slab(0,1,t_params(:,2:1:-1),comm_block,phi_s,phi)
#endif
    call fft(arrplan(2),phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
    !$OMP PARALLEL WORKSHARE
    phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)) = phi(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))*normfft
    !$OMP END PARALLEL WORKSHARE
    call finalize_n_solvers(hif_a(idir)-lo_a(idir)+1,asolver)
  end subroutine solve_impdiff_field_fft
#else
  subroutine solve_impdiff_field_3d(alphai_new,alphai_old,lo,hi,hif,is_bound,rhsx,rhsy,rhsz,dx2,dy2,dz2, &
                                    hypre_maxiter,hypre_tol,hypre_solver_i,phi,phio,asolver)
    !
    ! solves an implicit-diffusion Helmholtz problem for one field
    !
    implicit none
    real(rp), intent(in) :: alphai_new,alphai_old
    integer , intent(in), dimension(3) :: lo,hi,hif
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(lo(2):,lo(3):,0:) :: rhsx
    real(rp), intent(in), dimension(lo(1):,lo(3):,0:) :: rhsy
    real(rp), intent(in), dimension(lo(1):,lo(2):,0:) :: rhsz
    real(rp), intent(in), dimension(lo(1)-1:) :: dx2
    real(rp), intent(in), dimension(lo(2)-1:) :: dy2
    real(rp), intent(in), dimension(lo(3)-1:) :: dz2
    integer , intent(in) :: hypre_maxiter,hypre_solver_i
    real(rp), intent(in) :: hypre_tol
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: phi,phio
    type(hypre_solver), intent(inout) :: asolver
    integer :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(shared) PRIVATE(i,j,k) COLLAPSE(3)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          phi(i,j,k) = phi(i,j,k)*alphai_new*dx2(i)*dy2(j)*dz2(k)
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call updt_rhs(lo,hif,is_bound,rhsx,rhsy,rhsz,phi)
    call add_weighted_constant_to_diagonal(lo,hif,dx2,dy2,dz2,alphai_new-alphai_old,asolver%mat)
    call create_solver(hypre_maxiter,hypre_tol,.true.,hypre_solver_i,asolver)
    call setup_solver(asolver)
    call solve_helmholtz(asolver,lo,hif,phi,phio)
    call finalize_solver(asolver)
  end subroutine solve_impdiff_field_3d
#endif
#endif
end module mod_solve_helmholtz
