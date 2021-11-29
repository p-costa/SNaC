module mod_initflow
  use mpi_f08
  use mod_common_mpi, only: myid_block,ierr,comm_block
  use mod_param     , only: pi
  use mod_types
  implicit none
  private
  public initflow,init_inflow
  contains
  subroutine initflow(inivel,is_wallturb,lo,hi,lo_g,hi_g,l,uref,lref,visc,bforce, &
                      xc,xf,yc,yf,zc,zf,dxc,dxf,dyc,dyf,dzc,dzf,u,v,w,p)
    !
    ! computes initial conditions for the velocity field
    !
    implicit none
    character(len=3), intent(in) :: inivel
    logical , intent(in) :: is_wallturb
    integer , intent(in), dimension(3) :: lo,hi,lo_g,hi_g
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(inout) :: uref
    real(rp), intent(in) :: lref,visc,bforce
    real(rp), intent(in), dimension(lo(1)-1:) :: xc,xf,dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: yc,yf,dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: zc,zf,dzc,dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(out) :: u,v,w,p
    real(rp), allocatable, dimension(:) :: u1d
    integer  :: i,j,k
    logical  :: is_noise,is_mean,is_pair
    real(rp) :: reb,retau
    real(rp) :: xcl,xfl,ycl,yfl,zcl,zfl
    !
    allocate(u1d(lo(3):hi(3)))
    is_noise = .false.
    is_mean  = .false.
    is_pair  = .false.
    select case(trim(inivel))
    case('cou')
      call couette(lo(3),hi(3),zc,l(3),uref,u1d)
    case('poi')
      call poiseuille(lo(3),hi(3),zc,l(3),uref,u1d)
      is_mean=.true.
    case('zer')
      u1d(:) = 0._rp
    case('uni')
      u1d(:) = uref
    case('log')
      call log_profile(lo(3),hi(3),zc,l(3),uref,lref,visc,u1d)
      is_noise = .true.
      is_mean = .true.
    case('hcp')
      call poiseuille(lo(3),hi(3),zc,2.*l(3),uref,u1d)
      is_mean = .true.
    case('hcl')
      call log_profile(lo(3),hi(3),zc,2.*l(3),uref,lref,visc,u1d)
      is_noise = .true.
      is_mean=.true.
    case('tgv')
      do k=lo(3),hi(3)
        zcl = zc(k)/l(3)*2._rp*pi
        do j=lo(2),hi(2)
          ycl = yc(j)/l(2)*2._rp*pi
          yfl = yf(j)/l(2)*2._rp*pi
          do i=lo(1),hi(1)
            xcl = xc(i)/l(1)*2._rp*pi
            xfl = xf(i)/l(1)*2._rp*pi
            u(i,j,k) =  sin(xfl)*cos(ycl)*cos(zcl)
            v(i,j,k) = -cos(xcl)*sin(yfl)*cos(zcl)
            w(i,j,k) = 0._rp
            p(i,j,k) = 0._rp!(cos(2._rp*xc)+cos(2._rp*yc))*(cos(2._rp*zc)+2._rp)/16._rp
          end do
        end do
      end do
    case('pdc')
      if(is_wallturb) then ! turbulent flow
        retau  = sqrt(bforce*lref)*uref/visc
        reb    = (retau/.09_rp)**(1._rp/.88_rp)
        uref   = (reb/2._rp)/retau
      else                 ! laminar flow
        uref = (bforce*lref**2/(3._rp*visc))
      end if
      call poiseuille(lo(3),hi(3),zc,l(3),uref,u1d)
      is_mean=.true.
    case default
      if(myid_block == 0) write(stderr,*) 'ERROR: invalid name for initial velocity field'
      if(myid_block == 0) write(stderr,*) ''
      if(myid_block == 0) write(stderr,*) '*** Simulation aborted due to errors in the case file ***'
      if(myid_block == 0) write(stderr,*) '    check INFO_INPUT.md'
      call MPI_FINALIZE()
      error stop
    end select
    if(inivel /= 'tgv') then
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            u(i,j,k) = u1d(k)
            v(i,j,k) = 0._rp
            w(i,j,k) = 0._rp
            p(i,j,k) = 0._rp
          end do
        end do
      end do
    end if
    if(is_noise) then
      call add_noise(lo,hi,lo_g,hi_g,123,.5_rp,u)
      call add_noise(lo,hi,lo_g,hi_g,456,.5_rp,v)
      call add_noise(lo,hi,lo_g,hi_g,789,.5_rp,w)
    end if
    if(is_mean) then
      call set_mean(lo,hi,l,dxc,dyf,dzf,uref,u)
    end if
    if(is_wallturb) is_pair = .true.
    if(is_pair) then
      !
      ! initialize a streamwise vortex pair for a fast transition
      ! to turbulence in a pressure-driven channel:
      !        psi(x,y,z)  = f(z)*g(x,y), with
      !        f(z)        = (1-z**2)**2, and
      !        g(x,y)      = y*exp[-(16x**2-4y**2)]
      ! (x,y,z) --> (streamwise, spanwise, wall-normal) directions
      !
      ! see Henningson and Kim, JFM 1991
      !
      do k=lo(3),hi(3)
        zcl = 2._rp*zc(k)/l(3) - 1._rp ! z rescaled to be between -1 and +1
        zfl = 2._rp*zf(k)/l(3) - 1._rp
        do j=lo(2),hi(2)
          ycl = 2._rp*yc(j)/l(2) - 1._rp ! y rescaled to be between -1 and +1
          yfl = 2._rp*yf(j)/l(2) - 1._rp
          do i=lo(1),hi(1)
            xcl = 2._rp*xc(i)/l(1) - 1._rp ! x rescaled to be between -1 and +1
            xfl = 2._rp*xf(i)/l(1) - 1._rp
            v(i,j,k) = -1._rp*gxy(yfl,xcl)*dfz(zcl)*uref
            w(i,j,k) =  1._rp*fz(zfl)*dgxy(ycl,xcl)*uref
            p(i,j,k) =  0._rp
          end do
        end do
      end do
      !
      ! alternatively, using a Taylor-Green vortex
      ! for the cross-stream velocity components
      ! (commented below)
      !
      !do k=lo(3),hi(3)
      !  zcl = zc(k)/l(3)*2._rp*pi
      !  zfl = zf(k)/l(3)*2._rp*pi
      !  do j=lo(2),hi(2)
      !    ycl = yc(j)/l(2)*2._rp*pi
      !    yfl = yf(j)/l(2)*2._rp*pi
      !    do i=lo(1),hi(1)
      !      xcl = xc(i)/l(1)*2._rp*pi
      !      xfl = xf(i)/l(1)*2._rp*pi
      !      v(i,j,k) =  sin(xcl)*cos(yfl)*cos(zcl)
      !      w(i,j,k) = -cos(xcl)*sin(ycl)*cos(zfl)
      !      p(i,j,k) = 0._rp!(cos(2._rp*xcl)+cos(2._rp*ycl))*(cos(2._rp*zcl)+2._rp)/16._rp
      !    end do
      !  end do
      !end do
    end if
    deallocate(u1d)
  end subroutine initflow
  !
  subroutine add_noise(lo,hi,lo_g,hi_g,iseed,norm,p)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi,lo_g,hi_g
    integer , intent(in) :: iseed
    real(rp), intent(in) :: norm
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    integer(4), allocatable, dimension(:) :: seed
    real(rp) :: rn
    integer :: i,j,k
    allocate(seed(64))
    seed(:) = iseed
    call random_seed( put = seed )
    do k=lo_g(3),hi_g(3)
      do j=lo_g(2),hi_g(2)
        do i=lo_g(1),hi_g(1)
          call random_number(rn)
          if(i>=lo(1).and.i<=hi(1) .and. &
             j>=lo(2).and.j<=hi(2) .and. &
             k>=lo(3).and.k<=hi(3) ) then
             p(i,j,k) = p(i,j,k) + 2._rp*(rn-.5_rp)*norm
          end if
        end do
      end do
    end do
  end subroutine add_noise
  !
  subroutine set_mean(lo,hi,l,dx,dy,dz,mean,p)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(3) :: l
    real(rp), intent(in), dimension(lo(1)-1:) :: dx
    real(rp), intent(in), dimension(lo(2)-1:) :: dy
    real(rp), intent(in), dimension(lo(3)-1:) :: dz
    real(rp), intent(in) :: mean
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp) :: meanold
    integer :: i,j,k
    meanold = 0._rp
    !
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP SHARED(lo,hi,l,dx,dy,dz,p) &
    !$OMP PRIVATE(i,j,k) &
    !$OMP REDUCTION(+:meanold)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          meanold = meanold + p(i,j,k)*dx(i)*dy(j)*dz(k)/(l(1)*l(2)*l(3))
        end do
      end do
    end do
    !$OMP END PARALLEL DO
    call mpi_allreduce(MPI_IN_PLACE,meanold,1,MPI_REAL_RP,MPI_SUM,comm_block)
    if(meanold /= 0._rp) then
      !$OMP WORKSHARE
      p(:,:,:) = p(:,:,:)/meanold*mean
      !$OMP END WORKSHARE
    end if
  end subroutine set_mean
  !
  subroutine couette(lo,hi,zc,l,norm,p)
    !
    ! plane couette profile normalized by the wall velocity difference
    !
    implicit none
    integer , intent(in)   :: lo,hi
    real(rp), intent(in), dimension(lo-1:) :: zc
    real(rp), intent(in)   :: l,norm
    real(rp), intent(out), dimension(lo:) :: p
    integer :: k
    real(rp) :: z
    do k=lo,hi
      z    = zc(k)/l
      p(k) = .5_rp*(1._rp-2._rp*z)*norm
    end do
  end subroutine couette
  !
  subroutine poiseuille(lo,hi,zc,l,norm,p)
    implicit none
    integer , intent(in)   :: lo,hi
    real(rp), intent(in), dimension(lo-1:) :: zc
    real(rp), intent(in)   :: l,norm
    real(rp), intent(out), dimension(lo:) :: p
    integer :: k
    real(rp) :: z
    !
    ! plane poiseuille profile normalized by the bulk velocity
    !
    do k=lo,hi
      z    = zc(k)/l
      p(k) = 6._rp*z*(1._rp-z)*norm
    end do
  end subroutine poiseuille
  !
  function poiseuille_square(r,lmin,lmax) result(vel)
    !
    ! below a poiseuille inflow is calculated
    !   it is convinient that the integral a function of this kind
    !   is equal to one, such that it can be easily recycled
    !   for the inflow of a square duct, as done in
    !   in init_inflow below
    !
    real(rp), intent(in) :: r,lmin,lmax ! position where velocity is calculated
                                        ! and locations lower and upper boundary
    real(rp) :: vel,rr
    !
    rr = (r-lmin)/(lmax-lmin)
    vel = 6._rp*(rr*(1._rp-rr))
  end function poiseuille_square
  !
  subroutine init_inflow(periods,lo,hi,lmin,lmax,x1c,x2c,uref,vel)
    !
    ! below an inflow is calculated for a 2d plane
    ! later this subroutine can be generalized with other shapes of
    ! inflow velocity; by construction, the profile will not vary along
    ! a direction that does not have Dirichlet BCs; note that
    ! since the inflow is evaluated at the center of a face,
    ! there is no danger of experiencing a singularity '0._rp**0'.
    !
    implicit none
    integer , intent(in ), dimension(2       ) :: periods
    integer , intent(in ), dimension(2       ) :: lo,hi
    real(rp), intent(in ), dimension(2       ) :: lmin,lmax
    real(rp), intent(in ), dimension(lo(1)-1:) :: x1c
    real(rp), intent(in ), dimension(lo(2)-1:) :: x2c
    real(rp), intent(in ) :: uref ! target bulk velocity
    real(rp), intent(out), dimension(lo(1)-1:,lo(2)-1:) :: vel
    integer, dimension(2) :: q
    integer :: i1,i2
    ! real(rp), external :: poiseuille_square
    ! procedure (), pointer :: inflow_function => null()
    ! if(inflow_type == POISEUILLE) inflow_function => poiseuille_square
    !
    q(:) = 1
    where(periods(:) /= 0) q(:) = 0
    do i2 = lo(2)-1,hi(2)+1
      do i1 = lo(1)-1,hi(1)+1
        vel(i1,i2) = uref*poiseuille_square(x1c(i1),lmin(1),lmax(1))**q(1) * &
                          poiseuille_square(x2c(i2),lmin(2),lmax(2))**q(2)
      end do
    end do
  end subroutine init_inflow
  !
  subroutine log_profile(lo,hi,zc,l,uref,lref,visc,p)
    implicit none
    integer , intent(in)   :: lo,hi
    real(rp), intent(in), dimension(lo-1:) :: zc
    real(rp), intent(in)   :: l,uref,lref,visc
    real(rp), intent(out), dimension(lo:) :: p
    integer :: k
    real(rp) :: z,reb,retau ! z/lz and bulk Reynolds number
    reb = lref*uref/visc
    retau = 0.09_rp*reb**(0.88_rp) ! from Pope's book
    do k=lo,hi
      z    = zc(k)/l
      if(z>0.5_rp) z = 1._rp-z
      z    = zc(k)*2._rp*retau
      p(k) = 2.5_rp*log(z) + 5.5_rp
      if (z<=11.6_rp) p(k)=z
    end do
  end subroutine log_profile
  !
  ! functions to initialize the streamwise vortex pair
  ! (explained above)
  !
  function fz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: fz
    fz = ((1._rp-zc**2)**2)
  end function
  !
  function dfz(zc)
  real(rp), intent(in) :: zc
  real(rp) :: dfz
    dfz = -4._rp*zc*(1._rp-zc**2)
  end function
  !
  function gxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: gxy
    gxy = yc*exp(-4._rp*(4._rp*xc**2+yc**2))
  end function
  !
  function dgxy(xc,yc)
  real(rp), intent(in) :: xc,yc
  real(rp) :: dgxy
    dgxy = exp(-4._rp*(4._rp*xc**2+yc**2))*(1._rp-8._rp*yc**2)
  end function
end module mod_initflow
