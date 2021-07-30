module mod_fft
  use iso_c_binding , only: C_INT
  use mod_common_mpi, only: ierr
  use mod_fftw_param
  use mod_types
  !$ use omp_lib
  private
  public eigenvalues,fftini,fftend,fft
  contains
  subroutine fftini(idir,n,bc,is_centered,arrplan,normfft)
    implicit none
    integer         , intent(in )                 :: idir
    integer         , intent(in ), dimension(3  ) :: n
    character(len=1), intent(in ), dimension(0:1) :: bc
    logical         , intent(in )                 :: is_centered
    type(C_PTR)     , intent(out), dimension(2  ) :: arrplan
    real(rp)        , intent(out)                 :: normfft
    real(rp), dimension(n(1),n(2),n(3)) :: arr
    type(fftw_iodim), dimension(1) :: iodim
    type(fftw_iodim), dimension(2) :: iodim_howmany
    integer  :: kind_fwd,kind_bwd
    real(rp), dimension(2) :: norm
    integer  :: q
#ifdef _SINGLE_PRECISION
    !$ call sfftw_init_threads(ierr)
    !$ call sfftw_plan_with_nthreads(omp_get_max_threads())
#else
    !$ call dfftw_init_threads(ierr)
    !$ call dfftw_plan_with_nthreads(omp_get_max_threads())
#endif
    q = 0
    if(bc(0)//bc(1) /= 'FF'.and.(.not.is_centered)) q = 1
    select case(idir)
    case(1)
      iodim(1)%n  = n(idir)-q
      iodim(1)%is = 1
      iodim(1)%os = 1
      iodim_howmany(1)%n  = n(2)
      iodim_howmany(1)%is = n(1)
      iodim_howmany(1)%os = n(1)
      iodim_howmany(2)%n  = n(3)
      iodim_howmany(2)%is = n(1)*n(2)
      iodim_howmany(2)%os = n(1)*n(2)
    case(2)
      iodim(1)%n  = n(idir)-q
      iodim(1)%is = n(1)
      iodim(1)%os = n(1)
      iodim_howmany(1)%n  = n(1)
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1
      iodim_howmany(2)%n  = n(3)
      iodim_howmany(2)%is = n(1)*n(2)
      iodim_howmany(2)%os = n(1)*n(2)
    case(3)
      iodim(1)%n  = n(idir)-q
      iodim(1)%is = n(1)*n(2)
      iodim(1)%os = n(1)*n(2)
      iodim_howmany(1)%n  = n(1)
      iodim_howmany(1)%is = 1
      iodim_howmany(1)%os = 1
      iodim_howmany(2)%n  = n(2)
      iodim_howmany(2)%is = n(1)
      iodim_howmany(2)%os = n(1)
    end select
    call find_fft(bc,is_centered,kind_fwd,kind_bwd,norm)
    arrplan(1)=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arr,arr,kind_fwd,FFTW_ESTIMATE)
    arrplan(2)=fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,arr,arr,kind_bwd,FFTW_ESTIMATE)
    normfft = (norm(1)*(n(idir)+norm(2)-q))**(-1)
  end subroutine fftini
  !
  subroutine fftend(arrplan)
    implicit none
    type(C_PTR), intent(in), dimension(2) :: arrplan
#ifdef _SINGLE_PRECISION
    call sfftw_destroy_plan(arrplan(1))
    call sfftw_destroy_plan(arrplan(2))
    !$call sfftw_cleanup_threads(ierr)
#else
    call dfftw_destroy_plan(arrplan(1))
    call dfftw_destroy_plan(arrplan(2))
    !$call dfftw_cleanup_threads(ierr)
#endif
  end subroutine fftend
  !
  subroutine fft(plan,arr)
    implicit none
    type(C_PTR), intent(in   )                   :: plan
    real(rp)   , intent(inout), dimension(:,:,:) :: arr
#ifdef _SINGLE_PRECISION
    call sfftw_execute_r2r(plan,arr,arr)
#else
    call dfftw_execute_r2r(plan,arr,arr)
#endif
  end subroutine fft
  !
  subroutine find_fft(bc,is_centered,kind_fwd,kind_bwd,norm)
  implicit none
  character(len=1), intent(in ), dimension(0:1) :: bc
  logical         , intent(in )                 :: is_centered
  integer         , intent(out)                 :: kind_fwd,kind_bwd
  real(rp)        , intent(out), dimension(2  ) :: norm
  if(is_centered) then
    select case(bc(0)//bc(1))
    case('FF')
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = [1.,0.]*1._rp
    case('NN')
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
      norm = [2.,0.]*1._rp
    case('DD')
      kind_fwd = FFTW_RODFT10
      kind_bwd = FFTW_RODFT01
      norm = [2.,0.]*1._rp
    case('ND')
      kind_fwd = FFTW_REDFT11
      kind_bwd = FFTW_REDFT11
      norm = [2.,0.]*1._rp
    case('DN')
      kind_fwd = FFTW_RODFT11
      kind_bwd = FFTW_RODFT11
      norm = [2.,0.]*1._rp
    end select
  else
    select case(bc(0)//bc(1))
    case('FF')
      kind_fwd = FFTW_R2HC
      kind_bwd = FFTW_HC2R
      norm = [1.,0.]*1._rp
    case('NN')
      kind_fwd = FFTW_REDFT00
      kind_bwd = FFTW_REDFT00
      norm = [2.,-1.]*1._rp
    case('DD')
      kind_fwd = FFTW_RODFT00
      kind_bwd = FFTW_RODFT00
      norm = [2.,1.]*1._rp
    case('ND')
      kind_fwd = FFTW_REDFT10
      kind_bwd = FFTW_REDFT01
      norm = [2.,0.]*1._rp
    case('DN')
      kind_fwd = FFTW_RODFT01
      kind_bwd = FFTW_RODFT10
      norm = [2.,0.]*1._rp
    end select
  endif
  end subroutine find_fft
  !
  subroutine eigenvalues(n,bc,is_centered,lambda)
  use mod_param, only: pi
    implicit none
    integer , intent(in ) :: n
    character(len=1), intent(in ), dimension(0:1) :: bc
    logical         , intent(in )                 :: is_centered
    real(rp)        , intent(out), dimension(n  ) :: lambda
    integer :: l 
    select case(bc(0)//bc(1))
    case('FF')
      do l=1,n
        lambda(l) = -4._rp*sin((1._rp*(l-1))*pi/(1._rp*n))**2
      enddo
    case('NN')
      if(is_centered) then
        do l=1,n
          lambda(l) = -4._rp*sin((1._rp*(l-1))*pi/(2._rp*n))**2
        enddo
      else
        do l=1,n
          lambda(l) = -4._rp*sin((1._rp*(l-1))*pi/(2._rp*(n-1+1)))**2
        enddo
      endif
    case('DD')
      if(is_centered) then
        do l=1,n
          lambda(l) = -4._rp*sin((1._rp*(l-0))*pi/(2._rp*n))**2
        enddo
      else
        do l=1,n-1
          lambda(l) = -4._rp*sin((1._rp*(l-0))*pi/(2._rp*(n+1-1)))**2
        enddo
      endif
    case('ND')
      do l=1,n
        lambda(l) = -4._rp*sin((1._rp*(2*l-1))*pi/(4._rp*n))**2
      enddo
    end select   
  end subroutine eigenvalues
end module mod_fft
