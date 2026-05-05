module mod_scal
  use, intrinsic :: iso_c_binding, only: C_PTR
  use mpi_f08
  use mod_common_mpi, only: ierr
#ifdef _IMPDIFF
  use mod_solver, only: hypre_solver
#endif
  use mod_types
  implicit none
  private
  public scalar,initialize_scalars,scal_a,scal_d,bulk_mean_s,bulk_forcing_s
  !
  ! scalar field data
  !
#ifdef _IMPDIFF
  type scalar_rhs_bound
    real(rp), allocatable, dimension(:,:,:) :: x
    real(rp), allocatable, dimension(:,:,:) :: y
    real(rp), allocatable, dimension(:,:,:) :: z
  end type scalar_rhs_bound
#endif
  type scalar
    real(rp), allocatable, dimension(:,:,:) :: val
    real(rp), allocatable, dimension(:,:,:) :: dsdtrko
    real(rp) :: alpha
    character(len=100) :: ini
    character(len=1), dimension(0:1,3) :: cbc
    real(rp),         dimension(0:1,3) :: bc
    real(rp) :: source
    logical  :: is_forced
    real(rp) :: scalf
    real(rp) :: f
    real(rp), dimension(0:1,3) :: fluxo
#ifdef _IMPDIFF
    real(rp), allocatable, dimension(:,:,:) :: val_o
    type(scalar_rhs_bound) :: rhs
    type(hypre_solver) :: solver
    real(rp) :: alphai_o
#if defined(_FFT_X) || defined(_FFT_Y) || defined(_FFT_Z)
    type(C_PTR), dimension(2) :: arrplan
    real(rp) :: normfft
    real(rp), allocatable, dimension(:) :: lambda,lambda_a
    type(hypre_solver), allocatable, dimension(:) :: solver_fft
#ifdef _FFT_USE_SLABS
    real(rp), allocatable, dimension(:,:,:) :: val_s
#endif
#endif
#endif
  end type scalar
  contains
  subroutine initialize_scalars(scalars,nscal,lo,hi)
    use mod_param, only: alphai,iniscal,cbcscal,bcscal,ssource,is_sforced,scalf
    !
    ! initializes/allocates members of an array of scalar derived types
    !
    implicit none
    type(scalar), intent(inout), dimension(:) :: scalars
    integer     , intent(in   ) :: nscal
    integer     , intent(in   ), dimension(3) :: lo,hi
    integer :: iscal
    do iscal=1,nscal
      allocate(scalars(iscal)%val(    lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1), &
               scalars(iscal)%dsdtrko(lo(1)  :hi(1)  ,lo(2)  :hi(2)  ,lo(3)  :hi(3)  ))
#ifdef _IMPDIFF
      allocate(scalars(iscal)%val_o(lo(1)-1:hi(1)+1,lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1))
      allocate(scalars(iscal)%rhs%x(lo(2):hi(2),lo(3):hi(3),0:1), &
               scalars(iscal)%rhs%y(lo(1):hi(1),lo(3):hi(3),0:1), &
               scalars(iscal)%rhs%z(lo(1):hi(1),lo(2):hi(2),0:1))
      scalars(iscal)%val_o(:,:,:) = 0._rp
      scalars(iscal)%rhs%x(:,:,:) = 0._rp
      scalars(iscal)%rhs%y(:,:,:) = 0._rp
      scalars(iscal)%rhs%z(:,:,:) = 0._rp
      scalars(iscal)%alphai_o = 0._rp
#endif
      scalars(iscal)%val(:,:,:) = 0._rp
      scalars(iscal)%dsdtrko(:,:,:) = 0._rp
      scalars(iscal)%alpha     = alphai(iscal)**(-1)
      scalars(iscal)%ini       = iniscal(iscal)
      scalars(iscal)%cbc       = cbcscal(:,:,iscal)
      scalars(iscal)%bc        = bcscal(:,:,iscal)
      scalars(iscal)%source    = ssource(iscal)
      scalars(iscal)%is_forced = is_sforced(iscal)
      scalars(iscal)%scalf     = scalf(iscal)
      scalars(iscal)%f         = 0._rp
      scalars(iscal)%fluxo(:,:) = 0._rp
    end do
  end subroutine initialize_scalars
  !
  subroutine scal_d(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,s,dsdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(in) :: alpha
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in) :: s
    real(rp), dimension(lo(1):  ,lo(2)  :,lo(3):  ), intent(inout) :: dsdt
    real(rp) :: dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(dsdxp,dsdxm,dsdyp,dsdym,dsdzp,dsdzm) &
    !$OMP SHARED(lo,hi,dxc,dxf,dyc,dyf,dzc,dzf,alpha,s,dsdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          dsdxp = (s(i+1,j,k)-s(i  ,j,k))/dxc(i  )
          dsdxm = (s(i  ,j,k)-s(i-1,j,k))/dxc(i-1)
          dsdyp = (s(i,j+1,k)-s(i,j  ,k))/dyc(j  )
          dsdym = (s(i,j  ,k)-s(i,j-1,k))/dyc(j-1)
          dsdzp = (s(i,j,k+1)-s(i,j,k  ))/dzc(k  )
          dsdzm = (s(i,j,k  )-s(i,j,k-1))/dzc(k-1)
          !
          dsdt(i,j,k) = dsdt(i,j,k) + &
                        alpha*(dsdxp-dsdxm)/dxf(i) + &
                        alpha*(dsdyp-dsdym)/dyf(j) + &
                        alpha*(dsdzp-dsdzm)/dzf(k)
        end do
      end do
    end do
  end subroutine scal_d
  !
  subroutine scal_a(lo,hi,dxf,dyf,dzf,u,v,w,s,dsdt)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(in   ) :: u,v,w,s
    real(rp), dimension(lo(1):  ,lo(2):  ,lo(3):  ), intent(inout) :: dsdt
    real(rp) :: usip,usim,vsjp,vsjm,wskp,wskm
    integer  :: i,j,k
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(usip,usim,vsjp,vsjm,wskp,wskm) &
    !$OMP SHARED(lo,hi,dxf,dyf,dzf,u,v,w,s,dsdt)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          !
          usim  = 0.5_rp*( s(i-1,j,k)+s(i,j,k) )*u(i-1,j,k)
          usip  = 0.5_rp*( s(i+1,j,k)+s(i,j,k) )*u(i  ,j,k)
          vsjm  = 0.5_rp*( s(i,j-1,k)+s(i,j,k) )*v(i,j-1,k)
          vsjp  = 0.5_rp*( s(i,j+1,k)+s(i,j,k) )*v(i,j  ,k)
          wskm  = 0.5_rp*( s(i,j,k-1)+s(i,j,k) )*w(i,j,k-1)
          wskp  = 0.5_rp*( s(i,j,k+1)+s(i,j,k) )*w(i,j,k  )
          dsdt(i,j,k) = dsdt(i,j,k) + &
                        ( -usip + usim )/dxf(i) + &
                        ( -vsjp + vsjm )/dyf(j) + &
                        ( -wskp + wskm )/dzf(k)
        end do
      end do
    end do
  end subroutine scal_a
  !
  subroutine bulk_mean_s(lo,hi,vol_all,dxf,dyf,dzf,s,mean)
    implicit none
    integer , intent(in ), dimension(3) :: lo,hi
    real(rp), intent(in ) :: vol_all
    real(rp), intent(in ), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in ), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in ), dimension(lo(3)-1:) :: dzf
    real(rp), intent(in ), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: s
    real(rp), intent(out) :: mean
    integer :: i,j,k
    mean = 0._rp
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,dxf,dyf,dzf,s) &
    !$OMP REDUCTION(+:mean)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          mean = mean + s(i,j,k)*dxf(i)*dyf(j)*dzf(k)
        end do
      end do
    end do
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
    mean = mean/vol_all
  end subroutine bulk_mean_s
  !
  subroutine bulk_forcing_s(lo,hi,is_forced,ff,s)
    implicit none
    integer , intent(in   ), dimension(3) :: lo,hi
    logical , intent(in   ) :: is_forced
    real(rp), intent(in   ) :: ff
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: s
    integer :: i,j,k
    if(is_forced) then
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP SHARED(lo,hi,ff,s)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            s(i,j,k) = s(i,j,k) + ff
          end do
        end do
      end do
    end if
  end subroutine bulk_forcing_s
end module mod_scal
