module mod_bound
  use mpi
  use mod_types
  implicit none
  private
  public boundp,bounduvw,updt_rhs_b
  contains
  subroutine bounduvw(cbc,hi,lo,bc,isoutflow,halo,is_bound,nb, &
                      dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in), dimension(3) :: lo,hi
    real(rp), intent(in), dimension(0:1,3,3) :: bc
    logical , intent(in), dimension(0:1,3) :: isoutflow
    integer , intent(in), dimension(3    ) :: halo
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: nb
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    integer :: q,idir,sgn,ioutflowdir
    !
    call updthalo(lo,hi,1,halo(1),nb(:,1),1,u)
    call updthalo(lo,hi,1,halo(2),nb(:,2),2,u)
    call updthalo(lo,hi,1,halo(3),nb(:,3),3,u)
    call updthalo(lo,hi,1,halo(1),nb(:,1),1,v)
    call updthalo(lo,hi,1,halo(2),nb(:,2),2,v)
    call updthalo(lo,hi,1,halo(3),nb(:,3),3,v)
    call updthalo(lo,hi,1,halo(1),nb(:,1),1,w)
    call updthalo(lo,hi,1,halo(2),nb(:,2),2,w)
    call updthalo(lo,hi,1,halo(3),nb(:,3),3,w)
    !
    if(is_bound(0,1)) then
      call set_bc(cbc(0,1,1),0,lo,hi,1,.false.,bc(0,1,1),dxf(lo(1)-1),u)
      call set_bc(cbc(0,1,2),0,lo,hi,1,.true. ,bc(0,1,2),dxc(lo(1)-1),v)
      call set_bc(cbc(0,1,3),0,lo,hi,1,.true. ,bc(0,1,3),dxc(lo(1)-1),w)
    endif
    if(is_bound(1,1)) then
      call set_bc(cbc(1,1,1),1,lo,hi,1,.false.,bc(1,1,1),dxf(hi(1)  ),u)
      call set_bc(cbc(1,1,2),1,lo,hi,1,.true. ,bc(1,1,2),dxc(hi(1)  ),v)
      call set_bc(cbc(1,1,3),1,lo,hi,1,.true. ,bc(1,1,3),dxc(hi(1)  ),w)
    endif
    if(is_bound(0,2)) then
      call set_bc(cbc(0,2,1),0,lo,hi,2,.true. ,bc(0,2,1),dyc(lo(2)-1),u)
      call set_bc(cbc(0,2,2),0,lo,hi,2,.false.,bc(0,2,2),dyf(lo(2)-1),v)
      call set_bc(cbc(0,2,3),0,lo,hi,2,.true. ,bc(0,2,3),dyc(lo(2)-1),w)
     endif
    if(is_bound(1,2)) then
      call set_bc(cbc(1,2,1),1,lo,hi,2,.true. ,bc(1,2,1),dyc(hi(2)  ),u)
      call set_bc(cbc(1,2,2),1,lo,hi,2,.false.,bc(1,2,2),dyf(hi(2)  ),v)
      call set_bc(cbc(1,2,3),1,lo,hi,2,.true. ,bc(1,2,3),dyc(hi(2)  ),w)
    endif
    if(is_bound(0,3)) then
      call set_bc(cbc(0,3,1),0,lo,hi,3,.true. ,bc(0,3,1),dzc(lo(3)-1),u)
      call set_bc(cbc(0,3,2),0,lo,hi,3,.true. ,bc(0,3,2),dzc(lo(3)-1),v)
      call set_bc(cbc(0,3,3),0,lo,hi,3,.false.,bc(0,3,3),dzf(lo(3)-1),w)
    endif
    if(is_bound(1,3)) then
      call set_bc(cbc(1,3,1),1,lo,hi,3,.true. ,bc(1,3,1),dzc(hi(3)  ),u)
      call set_bc(cbc(1,3,2),1,lo,hi,3,.true. ,bc(1,3,2),dzc(hi(3)  ),v)
      call set_bc(cbc(1,3,3),1,lo,hi,3,.false.,bc(1,3,3),dzf(hi(3)  ),w)
    endif
    !
    do q = 1,3
      do idir = 0,1
        if(isoutflow(idir,q).and.is_bound(idir,q)) then
          if(idir.eq.0) sgn = -1
          if(idir.eq.1) sgn = +1
          ioutflowdir = q*sgn
          call outflow(lo,hi,ioutflowdir,dxf,dyf,dzf,u,v,w)
        endif
      enddo
    enddo
    return
  end subroutine bounduvw
  !
  subroutine boundp(cbc,lo,hi,bc,halo,is_bound,nb,dxc,dyc,dzc,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: lo,hi
    real(rp)            , intent(in), dimension(0:1,3) :: bc
    integer , intent(in), dimension(3    ) :: halo
    logical , intent(in), dimension(0:1,3) :: is_bound
    integer , intent(in), dimension(0:1,3) :: nb
    real(rp), intent(in), dimension(lo(1)-1:) :: dxc
    real(rp), intent(in), dimension(lo(2)-1:) :: dyc
    real(rp), intent(in), dimension(lo(3)-1:) :: dzc
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    !
    call updthalo(lo,hi,1,halo(1),nb(:,1),1,p)
    call updthalo(lo,hi,1,halo(2),nb(:,2),2,p)
    call updthalo(lo,hi,1,halo(3),nb(:,3),3,p)
    !
    if(is_bound(0,1)) then
      call set_bc(cbc(0,1),0,lo,hi,1,.true.,bc(0,1),dxc(lo(2)-1),p)
    endif                                                       
    if(is_bound(1,1)) then                                      
      call set_bc(cbc(1,1),1,lo,hi,1,.true.,bc(1,1),dxc(hi(2)  ),p)
    endif                       
    if(is_bound(0,2)) then      
      call set_bc(cbc(0,2),0,lo,hi,2,.true.,bc(0,2),dyc(lo(2)-1),p)
    endif                                                       
    if(is_bound(1,2)) then                                      
      call set_bc(cbc(1,2),1,lo,hi,2,.true.,bc(1,2),dyc(hi(2)  ),p)
    endif                       
    if(is_bound(0,3)) then      
      call set_bc(cbc(0,3),0,lo,hi,3,.true.,bc(0,3),dzc(lo(3)-1),p)
    endif                       
    if(is_bound(1,3)) then      
      call set_bc(cbc(1,3),1,lo,hi,3,.true.,bc(1,3),dzc(hi(3)  ),p)
    endif
    return
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,lo,hi,idir,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in) :: ibound
    integer , intent(in), dimension(3) :: lo,hi
    integer , intent(in) :: idir
    logical , intent(in) :: centered
    real(rp), intent(in) :: rvalue,dr
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp) :: factor,sgn
    !
    factor = rvalue
    if(ctype.eq.'D'.and.centered) then
      factor = 2._rp*factor
      sgn    = -1._rp
    endif
    if(ctype.eq.'N'.and.centered) then
      if(    ibound.eq.0) then
        factor = -dr*factor
      elseif(ibound.eq.1) then
        factor =  dr*factor
      endif
      sgn    = 1._rp
    endif
    !
    select case(ctype)
    case('P')
      select case(idir)
      case(1)
        !$OMP WORKSHARE
        p(lo(idir)-1,:,:) = p(hi(idir),:,:)
        p(hi(idir)+1,:,:) = p(lo(idir),:,:)
        !$OMP END WORKSHARE
      case(2)
        !$OMP WORKSHARE
        p(:,lo(idir)-1,:) = p(:,hi(idir),:)
        p(:,hi(idir)+1,:) = p(:,lo(idir),:)
        !$OMP END WORKSHARE
      case(3)
        !$OMP WORKSHARE
        p(:,:,lo(idir)-1) = p(:,:,hi(idir))
        p(:,:,hi(idir)+1) = p(:,:,lo(idir))
        !$OMP END WORKSHARE
      end select
    case('D','N')
      if(centered) then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(lo(idir)-1,:,:) = factor+sgn*p(lo(idir),:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(hi(idir)+1,:,:) = factor+sgn*p(hi(idir),:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,lo(idir)-1,:) = factor+sgn*p(:,lo(idir),:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,hi(idir)+1,:) = factor+sgn*p(:,hi(idir),:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,lo(idir)-1) = factor+sgn*p(:,:,lo(idir))
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,hi(idir)+1) = factor+sgn*p(:,:,hi(idir))
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'D') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(lo(idir)-1,:,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(hi(idir)  ,:,:) = factor
            p(hi(idir)+1,:,:) = p(hi(idir)-1,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,lo(idir)-1,:) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,hi(idir)  ,:) = factor
            p(:,hi(idir)+1,:) = p(:,hi(idir)-1,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            p(:,:,lo(idir)-1) = factor 
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            p(:,:,hi(idir)  ) = factor
            p(:,:,hi(idir)+1) = p(:,:,hi(idir)-1)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype.eq.'N') then
        select case(idir)
        case(1)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:))
            p(lo(idir)-1,:,:) = factor + p(lo(idir)  ,:,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
            p(hi(idir),:,:) = factor + p(hi(idir)-1,:,:)
            p(hi(idir)+1,:,:) = p(hi(idir),:,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
            p(:,lo(idir)-1,:) = factor + p(:,lo(idir)  ,:)
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
            p(:,hi(idir),:) = factor + p(:,hi(idir)-1,:)
            p(:,hi(idir)+1,:) = p(:,hi(idir),:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound.eq.0) then
            !$OMP WORKSHARE
            !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
            p(:,:,lo(idir)-1) = factor + p(:,:,lo(idir)  )
            !$OMP END WORKSHARE
          elseif(ibound.eq.1) then
            !$OMP WORKSHARE
            !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
            p(:,:,hi(idir)) = factor + p(:,:,hi(idir)-1)
            p(:,:,hi(idir)+1) = p(:,:,hi(idir)) ! not needed
            !$OMP END WORKSHARE
          endif
        end select
      endif
    end select
    return
  end subroutine set_bc
  !
  subroutine outflow(lo,hi,idir,dxf,dyf,dzf,u,v,w)
    implicit none
    integer , intent(in), dimension(3) :: lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(lo(1)-1:) :: dxf
    real(rp), intent(in), dimension(lo(2)-1:) :: dyf
    real(rp), intent(in), dimension(lo(3)-1:) :: dzf
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(inout) :: u,v,w
    integer :: i,j,k
    !
    ! determine face velocity from zero divergence
    !
    select case(idir)
    case(1) ! x direction, right
      i = hi(idir)
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(j,k) &
      !$OMP SHARED(n,i,u,v,w,dx,dyi,dzfi)
      do k=lo(3),hi(3)
        do j=hi(2),hi(2)
          u(i  ,j,k) = u(i-1,j,k) - dxf(i)*((v(i,j,k)-v(i,j-1,k))*dyf(j)+(w(i,j,k)-w(i,j,k-1))*dzf(k))
          u(i+1,j,k) = u(i  ,j,k) ! not needed
        enddo
      enddo
      !$OMP END PARALLEL DO
    case(2) ! y direction, back
      j = hi(idir)
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k) &
      !$OMP SHARED(n,j,u,v,w,dy,dxi,dzfi)
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          v(i,j  ,k) = v(i,j-1,k) - dyf(j)*((u(i,j,k)-u(i-1,j,k))*dxf(i)+(w(i,j,k)-w(i,j,k-1))*dzf(k))
          v(i,j+1,k) = v(i,j  ,k) ! not needed
        enddo
      enddo 
      !$OMP END PARALLEL DO
    case(3) ! z direction, top
      k = hi(idir)
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j) &
      !$OMP SHARED(n,k,u,v,w,dzf,dxi,dyi)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          w(i,j,k  ) = w(i,j,k-1) - dzf(k)*((u(i,j,k)-u(i-1,j,k))/dxf(i)+(v(i,j,k)-v(i,j-1,k))/dyf(j))
          w(i,j,k+1) = w(i,j,k  ) ! not needed
        enddo
      enddo 
      !$OMP END PARALLEL DO
    case(-1) ! x direction, left
      i = lo(idir)-1
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(j,k) &
      !$OMP SHARED(n,i,u,v,w,dx,dyi,dzfi)
      do k=lo(3),hi(3)
        do j=lo(2),hi(2)
          u(i,j,k) = u(i+1,j,k) + dxf(i)*((v(i+1,j,k)-v(i+1,j-1,k))/dyf(j)+(w(i+1,j,k)-w(i+1,j,k-1))/dzf(k))
        enddo
      enddo 
      !$OMP END PARALLEL DO
    case(-2) ! y direction, front
      j = lo(idir)-1
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,k) &
      !$OMP SHARED(n,j,u,v,w,dy,dxi,dzfi)
      do k=lo(3),hi(3)
        do i=lo(1),hi(1)
          v(i,j,k) = v(i,j+1,k) + dyf(j)*((u(i,j+1,k)-u(i-1,j+1,k))/dxf(i)+(w(i,j+1,k)-w(i,j+1,k-1))/dzf(k))
        enddo
      enddo 
      !$OMP END PARALLEL DO
    case(-3) ! z direction, bottom
      k = lo(idir)-1
      !$OMP PARALLEL DO DEFAULT(none) &
      !$OMP PRIVATE(i,j) &
      !$OMP SHARED(n,k,u,v,w,dzf,dxi,dyi)
      do j=lo(2),hi(2)
        do i=lo(1),hi(1)
          w(i,j,k) = w(i,j,k+1) + dzf(k)*((u(i,j,k+1)-u(i-1,j,k+1))/dxf(i)+(v(i,j,k+1)-v(i,j-1,k+1))/dyf(j))
        enddo
      enddo 
      !$OMP END PARALLEL DO
    end select
    return
  end subroutine outflow
  !
  subroutine inflow(lo,hi,idir,vel2d,u,v,w)
    !
    ! n.b.: if the call of inflow should be conditioned by is_bound
    !       just like the outflow subroutine
    !
    implicit none
    integer, intent(in), dimension(3) :: lo,hi
    integer, intent(in) :: idir
    real(rp), dimension(0:,0:), intent(in) :: vel2d
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(inout) :: u,v,w
    integer :: i,j,k
    !
    select case(idir)
      case(1) ! x direction
        i = lo(idir)-1
        do k=lo(3),hi(3)
          do j=lo(2),hi(2)
            u(i,j,k) = vel2d(j,k)
          enddo
        enddo 
      case(2) ! y direction
        j = lo(idir)-1
        do k=lo(3),hi(3)
          do i=lo(1),hi(1)
            v(i,j,k) = vel2d(i,k)
          enddo
        enddo 
      case(3) ! z direction
        k = lo(idir)-1
        do j=lo(2),hi(2)
          do i=lo(1),hi(1)
            w(i,j,k) = vel2d(i,j)
          enddo
        enddo 
    end select
    return
  end subroutine inflow
  !
  subroutine updt_rhs_b(c_or_f,cbc,lo,hi,is_bound,rhsbx,rhsby,rhsbz,p)
    implicit none
    character, intent(in), dimension(3) :: c_or_f
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in), dimension(3) :: lo,hi
    logical , intent(in), dimension(0:1,3) :: is_bound
    real(rp), intent(in), dimension(lo(2):,lo(3):,0:) :: rhsbx
    real(rp), intent(in), dimension(lo(1):,lo(3):,0:) :: rhsby
    real(rp), intent(in), dimension(lo(1):,lo(2):,0:) :: rhsbz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    integer , dimension(3) :: q
    integer :: idir
    q(:) = 0
    do idir = 1,3
      if(c_or_f(idir).eq.'f'.and.cbc(1,idir).eq.'D') q(idir) = 1
    enddo
    if(is_bound(0,1)) then
      !$OMP WORKSHARE
      p(lo(1)     ,lo(2):hi(2),lo(3):hi(3)) = p(lo(1)     ,lo(2):hi(2),lo(3):hi(3)) + rhsbx(:,:,0)
      !$OMP END WORKSHARE
    endif  
    if(is_bound(1,1)) then
      !$OMP WORKSHARE
      p(hi(1)-q(1),lo(2):hi(2),lo(3):hi(3)) = p(hi(1)-q(1),lo(2):hi(2),lo(3):hi(3)) + rhsbx(:,:,1)
      !$OMP END WORKSHARE
    endif
    if(is_bound(0,2)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),lo(2)     ,lo(3):hi(3)) = p(lo(1):hi(1),lo(2)     ,lo(3):hi(3)) + rhsby(:,:,0)
      !$OMP END WORKSHARE
    endif
    if(is_bound(1,2)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),hi(2)-q(2),lo(3):hi(3)) = p(lo(1):hi(1),hi(2)-q(2),lo(3):hi(3)) + rhsby(:,:,1)
      !$OMP END WORKSHARE
    endif
    if(is_bound(0,3)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),lo(2):hi(2),lo(3)     ) = p(lo(1):hi(1),lo(2):hi(2),lo(3)     ) + rhsbz(:,:,0)
      !$OMP END WORKSHARE
    endif
    if(is_bound(1,3)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),lo(2):hi(2),hi(3)-q(3)) = p(lo(1):hi(1),lo(2):hi(2),hi(3)-q(3)) + rhsbz(:,:,1)
      !$OMP END WORKSHARE
    endif
    return
  end subroutine updt_rhs_b
  !
  subroutine updthalo(lo,hi,n,halo,nb,idir,p)
    implicit none
    integer , dimension(3), intent(in) :: lo,hi
    integer , intent(in) :: n,halo
    integer , intent(in), dimension(0:1) :: nb
    integer , intent(in) :: idir
    real(rp), dimension(lo(1)-n:,lo(2)-n:,lo(3)-n:), intent(inout) :: p
    integer :: ierr,status(MPI_STATUS_SIZE)
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halos that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(lo(1)  ,lo(2)-n,lo(3)-n),1,halo,nb(0),0, &
                        p(hi(1)+n,lo(2)-n,lo(3)-n),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(p(hi(1)  ,lo(2)-n,lo(3)-n),1,halo,nb(1),0, &
                        p(lo(1)-n,lo(2)-n,lo(3)-n),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,status,ierr)
         !call MPI_IRECV( p(lo(1)-n,lo(2)-n,lo(3)-n),1,halo,nb(0),1, &
         !                MPI_COMM_WORLD,requests(1),error)
         !call MPI_IRECV( p(hi(1)+n,lo(2)-n,lo(3)-n),1,halo,nb(1),0, &
         !                MPI_COMM_WORLD,requests(2),error)
         !call MPI_ISSEND(p(hi(1)  ,lo(2)-n,lo(3)-n),1,halo,nb(1),1, &
         !                MPI_COMM_WORLD,requests(3),error)
         !call MPI_ISSEND(p(lo(1)  ,lo(2)-n,lo(3)-n),1,halo,nb(0),0, &
         !                MPI_COMM_WORLD,requests(4),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
      call MPI_SENDRECV(p(lo(1)-n,lo(2)  ,lo(3)-n),1,halo,nb(0),0, &
                        p(lo(1)-n,hi(2)+n,lo(3)-n),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(p(lo(1)-n,hi(2)  ,lo(3)-n),1,halo,nb(1),0, &
                        p(lo(1)-n,lo(2)-n,lo(3)-n),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,status,ierr)
         !call MPI_IRECV( p(lo(1)-n,hi(2)+n,lo(3)-n),1,halo,nb(1),0, &
         !                MPI_COMM_WORLD,requests(1),error)
         !call MPI_IRECV( p(lo(1)-n,lo(2)-n,lo(3)-n),1,halo,nb(0),1, &
         !                MPI_COMM_WORLD,requests(2),error)
         !call MPI_ISSEND(p(lo(1)-n,lo(2)  ,lo(3)-n),1,halo,nb(0),0, &
         !               MPI_COMM_WORLD,requests(3),error)
         !call MPI_ISSEND(p(lo(1)-n,hi(2)  ,lo(3)-n),1,halo,nb(1),1, &
         !               MPI_COMM_WORLD,requests(4),error)
         !call MPI_WAITALL(4, requests, statuses, error)
    case(3) ! z direction
      call MPI_SENDRECV(p(lo(1)-n,lo(2)-n,lo(3)  ),1,halo,nb(0),0, &
                        p(lo(1)-n,lo(2)-n,hi(3)+n),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(p(lo(1)-n,lo(2)-n,hi(3)  ),1,halo,nb(1),0, &
                        p(lo(1)-n,lo(2)-n,lo(3)-n),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,status,ierr)
      !call MPI_IRECV( p(lo(1)-n,lo(2)-n,hi(3)+n),1,halo,nb(1),0, &
      !                MPI_COMM_WORLD,requests(1),error)
      !call MPI_IRECV( p(lo(1)-n,lo(2)-n,lo(3)-n),1,halo,nb(0),1, &
      !                MPI_COMM_WORLD,requests(2),error)
      !call MPI_ISSEND(p(lo(1)-n,lo(2)-n,lo(3)  ),1,halo,nb(0),0, &
      !               MPI_COMM_WORLD,requests(3),error)
      !call MPI_ISSEND(p(lo(1)-n,lo(2)-n,hi(3)  ),1,halo,nb(1),1, &
      !               MPI_COMM_WORLD,requests(4),error)
      !call MPI_WAITALL(4, requests, statuses, error)
    end select
    return
  end subroutine updthalo
end module mod_bound
