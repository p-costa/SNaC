module mod_bound
  use mpi_f08
  use mod_types
  implicit none
  private
  public boundp,bounduvw,updt_rhs
  contains
  subroutine bounduvw(cbc,lo,hi,bc,is_correc,halos,is_bound,nb, &
                      dxc,dxf,dyc,dyf,dzc,dzf,u,v,w)
    !
    ! imposes velocity boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3,3) :: cbc
    integer , intent(in   ), dimension(3) :: lo,hi
    real(rp), intent(in   ), dimension(0:1,3,3) :: bc
    logical , intent(in   )                   :: is_correc
    type(MPI_DATATYPE), intent(in   ), dimension(3    ) :: halos
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    integer , intent(in   ), dimension(0:1,3) :: nb
    real(rp), intent(in   ), dimension(lo(1)-1:) :: dxc,dxf
    real(rp), intent(in   ), dimension(lo(2)-1:) :: dyc,dyf
    real(rp), intent(in   ), dimension(lo(3)-1:) :: dzc,dzf
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: u,v,w
    !
    call updthalo(lo,hi,1,halos(1),nb(:,1),1,u)
    call updthalo(lo,hi,1,halos(2),nb(:,2),2,u)
    call updthalo(lo,hi,1,halos(3),nb(:,3),3,u)
    call updthalo(lo,hi,1,halos(1),nb(:,1),1,v)
    call updthalo(lo,hi,1,halos(2),nb(:,2),2,v)
    call updthalo(lo,hi,1,halos(3),nb(:,3),3,v)
    call updthalo(lo,hi,1,halos(1),nb(:,1),1,w)
    call updthalo(lo,hi,1,halos(2),nb(:,2),2,w)
    call updthalo(lo,hi,1,halos(3),nb(:,3),3,w)
    !
    if(is_bound(0,1)) then
      if(.not.is_correc) call set_bc(cbc(0,1,1),0,lo,hi,1,.false.,bc(0,1,1),dxf(lo(1)-1),u)
                         call set_bc(cbc(0,1,2),0,lo,hi,1,.true. ,bc(0,1,2),dxc(lo(1)-1),v)
                         call set_bc(cbc(0,1,3),0,lo,hi,1,.true. ,bc(0,1,3),dxc(lo(1)-1),w)
    endif
    if(is_bound(1,1)) then
      if(.not.is_correc) call set_bc(cbc(1,1,1),1,lo,hi,1,.false.,bc(1,1,1),dxf(hi(1)  ),u)
                         call set_bc(cbc(1,1,2),1,lo,hi,1,.true. ,bc(1,1,2),dxc(hi(1)  ),v)
                         call set_bc(cbc(1,1,3),1,lo,hi,1,.true. ,bc(1,1,3),dxc(hi(1)  ),w)
    endif
    if(is_bound(0,2)) then
                         call set_bc(cbc(0,2,1),0,lo,hi,2,.true. ,bc(0,2,1),dyc(lo(2)-1),u)
      if(.not.is_correc) call set_bc(cbc(0,2,2),0,lo,hi,2,.false.,bc(0,2,2),dyf(lo(2)-1),v)
                         call set_bc(cbc(0,2,3),0,lo,hi,2,.true. ,bc(0,2,3),dyc(lo(2)-1),w)
     endif
    if(is_bound(1,2)) then
                         call set_bc(cbc(1,2,1),1,lo,hi,2,.true. ,bc(1,2,1),dyc(hi(2)  ),u)
      if(.not.is_correc) call set_bc(cbc(1,2,2),1,lo,hi,2,.false.,bc(1,2,2),dyf(hi(2)  ),v)
                         call set_bc(cbc(1,2,3),1,lo,hi,2,.true. ,bc(1,2,3),dyc(hi(2)  ),w)
    endif
    if(is_bound(0,3)) then
                         call set_bc(cbc(0,3,1),0,lo,hi,3,.true. ,bc(0,3,1),dzc(lo(3)-1),u)
                         call set_bc(cbc(0,3,2),0,lo,hi,3,.true. ,bc(0,3,2),dzc(lo(3)-1),v)
      if(.not.is_correc) call set_bc(cbc(0,3,3),0,lo,hi,3,.false.,bc(0,3,3),dzf(lo(3)-1),w)
    endif
    if(is_bound(1,3)) then
                         call set_bc(cbc(1,3,1),1,lo,hi,3,.true. ,bc(1,3,1),dzc(hi(3)  ),u)
                         call set_bc(cbc(1,3,2),1,lo,hi,3,.true. ,bc(1,3,2),dzc(hi(3)  ),v)
      if(.not.is_correc) call set_bc(cbc(1,3,3),1,lo,hi,3,.false.,bc(1,3,3),dzf(hi(3)  ),w)
    endif
  end subroutine bounduvw
  !
  subroutine boundp(cbc,lo,hi,bc,halos,is_bound,nb,dxc,dyc,dzc,p)
    !
    ! imposes pressure boundary conditions
    !
    implicit none
    character(len=1), intent(in), dimension(0:1,3) :: cbc
    integer , intent(in   ), dimension(3    )      :: lo,hi
    real(rp), intent(in   ), dimension(0:1,3)      :: bc
    type(MPI_DATATYPE), intent(in   ), dimension(3    )      :: halos
    logical , intent(in   ), dimension(0:1,3)      :: is_bound
    integer , intent(in   ), dimension(0:1,3)      :: nb
    real(rp), intent(in   ), dimension(lo(1)-1:)   :: dxc
    real(rp), intent(in   ), dimension(lo(2)-1:)   :: dyc
    real(rp), intent(in   ), dimension(lo(3)-1:)   :: dzc
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    !
    call updthalo(lo,hi,1,halos(1),nb(:,1),1,p)
    call updthalo(lo,hi,1,halos(2),nb(:,2),2,p)
    call updthalo(lo,hi,1,halos(3),nb(:,3),3,p)
    !
    if(is_bound(0,1)) then
      call set_bc(cbc(0,1),0,lo,hi,1,.true.,bc(0,1),dxc(lo(1)-1),p)
    endif                                                       
    if(is_bound(1,1)) then                                      
      call set_bc(cbc(1,1),1,lo,hi,1,.true.,bc(1,1),dxc(hi(1)  ),p)
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
  end subroutine boundp
  !
  subroutine set_bc(ctype,ibound,lo,hi,idir,centered,rvalue,dr,p)
    implicit none
    character(len=1), intent(in) :: ctype
    integer , intent(in   ) :: ibound
    integer , intent(in   ), dimension(3) :: lo,hi
    integer , intent(in   ) :: idir
    logical , intent(in   ) :: centered
    real(rp), intent(in   ) :: rvalue,dr
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp) :: factor,sgn
    !
    factor = rvalue
    if(ctype == 'D'.and.centered) then
      factor = 2._rp*factor
      sgn    = -1._rp
    endif
    if(ctype == 'N'.and.centered) then
      if(    ibound == 0) then
        factor = -dr*factor
      elseif(ibound == 1) then
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
          if    (ibound == 0) then
            !$OMP WORKSHARE
            p(lo(idir)-1,:,:) = factor+sgn*p(lo(idir),:,:)
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            p(hi(idir)+1,:,:) = factor+sgn*p(hi(idir),:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            p(:,lo(idir)-1,:) = factor+sgn*p(:,lo(idir),:)
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            p(:,hi(idir)+1,:) = factor+sgn*p(:,hi(idir),:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            p(:,:,lo(idir)-1) = factor+sgn*p(:,:,lo(idir))
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            p(:,:,hi(idir)+1) = factor+sgn*p(:,:,hi(idir))
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype == 'D') then
        select case(idir)
        case(1)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            p(lo(idir)-1,:,:) = factor
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            p(hi(idir)  ,:,:) = factor
            p(hi(idir)+1,:,:) = p(hi(idir)-1,:,:)
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            p(:,lo(idir)-1,:) = factor
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            p(:,hi(idir)  ,:) = factor
            p(:,hi(idir)+1,:) = p(:,hi(idir)-1,:)
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            p(:,:,lo(idir)-1) = factor
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            p(:,:,hi(idir)  ) = factor
            p(:,:,hi(idir)+1) = p(:,:,hi(idir)-1)
            !$OMP END WORKSHARE
          endif
        end select
      elseif(.not.centered.and.ctype == 'N') then
        select case(idir)
        case(1)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            !p(0,:,:) = 1./3.*(-2.*factor+4.*p(1  ,:,:)-p(2  ,:,:))
            p(lo(idir)-1,:,:) = factor + p(lo(idir)  ,:,:)
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            !p(n,:,:) = 1./3.*(-2.*factor+4.*p(n-1,:,:)-p(n-2,:,:))
            p(hi(idir)  ,:,:) = factor + p(hi(idir)-1,:,:)
            p(hi(idir)+1,:,:) =          p(hi(idir)  ,:,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(2)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            !p(:,0  ,:) = 1./3.*(-2.*factor+4.*p(:,1,:)-p(:,2  ,:))
            p(:,lo(idir)-1,:) = factor + p(:,lo(idir)  ,:)
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            !p(:,n,:) = 1./3.*(-2.*factor+4.*p(:,n-1,:)-p(:,n-2,:))
            p(:,hi(idir)  ,:) = factor + p(:,hi(idir)-1,:)
            p(:,hi(idir)+1,:) =          p(:,hi(idir)  ,:) ! not needed
            !$OMP END WORKSHARE
          endif
        case(3)
          if    (ibound == 0) then
            !$OMP WORKSHARE
            !p(:,:,0) = 1./3.*(-2.*factor+4.*p(:,:,1  )-p(:,:,2  ))
            p(:,:,lo(idir)-1) = factor + p(:,:,lo(idir)  )
            !$OMP END WORKSHARE
          elseif(ibound == 1) then
            !$OMP WORKSHARE
            !p(:,:,n) = 1./3.*(-2.*factor+4.*p(:,:,n-1)-p(:,:,n-2))
            p(:,:,hi(idir)  ) = factor + p(:,:,hi(idir)-1)
            p(:,:,hi(idir)+1) =          p(:,:,hi(idir)  ) ! not needed
            !$OMP END WORKSHARE
          endif
        end select
      endif
    end select
  end subroutine set_bc
  !
  subroutine inflow(lo,hi,lo_2d,ibound,idir,vel2d,p)
    !
    implicit none
    integer, intent(in), dimension(3) :: lo,hi
    integer, intent(in), dimension(2) :: lo_2d
    integer, intent(in) :: ibound,idir
    real(rp), dimension(lo_2d(1)-1:,lo_2d(2)-1:), intent(in) :: vel2d
    real(rp), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:), intent(inout) :: p
    integer :: q
    !
    q = (ibound-1)*(lo(idir)-1)+ibound*hi(idir)
    select case(idir)
      case(1)
        p(q,:,:) = vel2d(:,:)
      case(2)
        p(:,q,:) = vel2d(:,:)
      case(3)
        p(:,:,q) = vel2d(:,:)
    end select
  end subroutine inflow
  !
  subroutine updt_rhs(lo,hi,is_bound,rhsbx,rhsby,rhsbz,p)
    !
    ! updates the right-hand-side of the Helmholtz/Poisson equation
    ! with the appropriate boundary conditions
    !
    implicit none
    integer , intent(in   ), dimension(3) :: lo,hi
    logical , intent(in   ), dimension(0:1,3) :: is_bound
    real(rp), intent(in   ), dimension(lo(2):,lo(3):,0:) :: rhsbx
    real(rp), intent(in   ), dimension(lo(1):,lo(3):,0:) :: rhsby
    real(rp), intent(in   ), dimension(lo(1):,lo(2):,0:) :: rhsbz
    real(rp), intent(inout), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    if(is_bound(0,1)) then
      !$OMP WORKSHARE
      p(lo(1),lo(2):hi(2),lo(3):hi(3)) = p(lo(1),lo(2):hi(2),lo(3):hi(3)) + rhsbx(lo(2):hi(2),lo(3):hi(3),0)
      !$OMP END WORKSHARE
    endif  
    if(is_bound(1,1)) then
      !$OMP WORKSHARE
      p(hi(1),lo(2):hi(2),lo(3):hi(3)) = p(hi(1),lo(2):hi(2),lo(3):hi(3)) + rhsbx(lo(2):hi(2),lo(3):hi(3),1)
      !$OMP END WORKSHARE
    endif
    if(is_bound(0,2)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),lo(2),lo(3):hi(3)) = p(lo(1):hi(1),lo(2),lo(3):hi(3)) + rhsby(lo(1):hi(1),lo(3):hi(3),0)
      !$OMP END WORKSHARE
    endif
    if(is_bound(1,2)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),hi(2),lo(3):hi(3)) = p(lo(1):hi(1),hi(2),lo(3):hi(3)) + rhsby(lo(1):hi(1),lo(3):hi(3),1)
      !$OMP END WORKSHARE
    endif
    if(is_bound(0,3)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),lo(2):hi(2),lo(3)) = p(lo(1):hi(1),lo(2):hi(2),lo(3)) + rhsbz(lo(1):hi(1),lo(2):hi(2),0)
      !$OMP END WORKSHARE
    endif
    if(is_bound(1,3)) then
      !$OMP WORKSHARE
      p(lo(1):hi(1),lo(2):hi(2),hi(3)) = p(lo(1):hi(1),lo(2):hi(2),hi(3)) + rhsbz(lo(1):hi(1),lo(2):hi(2),1)
      !$OMP END WORKSHARE
    endif
  end subroutine updt_rhs
  !
  subroutine updthalo(lo,hi,nh,halo,nb,idir,p)
    implicit none
    integer , dimension(3), intent(in) :: lo,hi
    integer , intent(in) :: nh ! number of ghost points
    type(MPI_DATATYPE), intent(in) :: halo
    integer , intent(in), dimension(0:1) :: nb
    integer , intent(in) :: idir
    real(rp), dimension(lo(1)-nh:,lo(2)-nh:,lo(3)-nh:), intent(inout) :: p
    !integer :: requests(4), statuses(MPI_STATUS_SIZE,4)
    !
    !  this subroutine updates the halo that store info
    !  from the neighboring computational sub-domain
    !
    select case(idir)
    case(1) ! x direction
      call MPI_SENDRECV(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call MPI_SENDRECV(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !call MPI_IRECV( p(lo(1)-nh  ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),1, &
      !                MPI_COMM_WORLD,requests(1),error)
      !call MPI_IRECV( p(hi(1)+1   ,lo(2)-nh,lo(3)-nh),1,halo,nb(1),0, &
      !                MPI_COMM_WORLD,requests(2),error)
      !call MPI_ISSEND(p(hi(1)-nh+1,lo(2)-nh,lo(3)-nh),1,halo,nb(1),1, &
      !                MPI_COMM_WORLD,requests(3),error)
      !call MPI_ISSEND(p(lo(1)     ,lo(2)-nh,lo(3)-nh),1,halo,nb(0),0, &
      !                MPI_COMM_WORLD,requests(4),error)
      !call MPI_WAITALL(4, requests, statuses, error)
    case(2) ! y direction
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
                        p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call MPI_SENDRECV(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !call MPI_IRECV( p(lo(1)-nh,hi(2)+1   ,lo(3)-nh),1,halo,nb(1),0, &
      !                MPI_COMM_WORLD,requests(1),error)
      !call MPI_IRECV( p(lo(1)-nh,lo(2)-nh  ,lo(3)-nh),1,halo,nb(0),1, &
      !                MPI_COMM_WORLD,requests(2),error)
      !call MPI_ISSEND(p(lo(1)-nh,lo(2)     ,lo(3)-nh),1,halo,nb(0),0, &
      !               MPI_COMM_WORLD,requests(3),error)
      !call MPI_ISSEND(p(lo(1)-nh,hi(2)-nh+1,lo(3)-nh),1,halo,nb(1),1, &
      !               MPI_COMM_WORLD,requests(4),error)
      !call MPI_WAITALL(4, requests, statuses, error)
    case(3) ! z direction
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
                        p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      call MPI_SENDRECV(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),0, &
                        p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),0, &
                        MPI_COMM_WORLD,MPI_STATUS_IGNORE)
      !call MPI_IRECV( p(lo(1)-nh,lo(2)-nh,hi(3)+1   ),1,halo,nb(1),0, &
      !                MPI_COMM_WORLD,requests(1),error)
      !call MPI_IRECV( p(lo(1)-nh,lo(2)-nh,lo(3)-nh  ),1,halo,nb(0),1, &
      !                MPI_COMM_WORLD,requests(2),error)
      !call MPI_ISSEND(p(lo(1)-nh,lo(2)-nh,lo(3)     ),1,halo,nb(0),0, &
      !               MPI_COMM_WORLD,requests(3),error)
      !call MPI_ISSEND(p(lo(1)-nh,lo(2)-nh,hi(3)-nh+1),1,halo,nb(1),1, &
      !               MPI_COMM_WORLD,requests(4),error)
      !call MPI_WAITALL(4, requests, statuses, error)
    end select
  end subroutine updthalo
end module mod_bound
