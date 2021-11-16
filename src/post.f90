module mod_post
  use mod_types
  implicit none
  private
  public cmpt_vorticity,cmpt_rotation_rate,cmpt_strain_rate,cmpt_q_criterion, &
         cmpt_wall_forces, updt_wall_forces
  contains
  subroutine cmpt_vorticity(n,dxc,dyc,dzc,ux,uy,uz,vox,voy,voz)
    !
    ! computes the vorticity field
    !
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:)       :: dxc,dyc,dzc
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux ,uy ,uz
    real(rp), intent(out), dimension( :, :, :) :: vox,voy,voz
    integer :: i,j,k
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(n,dxc,dyc,dzc,ux,uy,uz,vox,voy,voz) &
    !$OMP PRIVATE(i,j,k)
    !$OMP DO
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
           !
           ! x component of the vorticity at cell center
           !
           vox(i,j,k) = 0.25_rp*( &
                                 (uz(i,j+1,k  )-uz(i,j  ,k  ))/dyc(j  ) - (uy(i,j  ,k+1)-uy(i,j  ,k  ))/dzc(k  ) + &
                                 (uz(i,j+1,k-1)-uz(i,j  ,k-1))/dyc(j  ) - (uy(i,j  ,k  )-uy(i,j  ,k-1))/dzc(k-1) + &
                                 (uz(i,j  ,k  )-uz(i,j-1,k  ))/dyc(j-1) - (uy(i,j-1,k+1)-uy(i,j-1,k  ))/dzc(k  ) + &
                                 (uz(i,j  ,k-1)-uz(i,j-1,k-1))/dyc(j-1) - (uy(i,j-1,k  )-uy(i,j-1,k-1))/dzc(k-1) &
                                )
           !
           ! y component of the vorticity at cell center
           !
           voy(i,j,k) = 0.25_rp*( &
                                 (ux(i  ,j,k+1)-ux(i  ,j,k  ))/dzc(k  ) - (uz(i+1,j,k  )-uz(i  ,j,k  ))/dxc(i  ) + &
                                 (ux(i  ,j,k  )-ux(i  ,j,k-1))/dzc(k-1) - (uz(i+1,j,k-1)-uz(i  ,j,k-1))/dxc(i  ) + &
                                 (ux(i-1,j,k+1)-ux(i-1,j,k  ))/dzc(k  ) - (uz(i  ,j,k  )-uz(i-1,j,k  ))/dxc(i-1) + &
                                 (ux(i-1,j,k  )-ux(i-1,j,k-1))/dzc(k-1) - (uz(i  ,j,k-1)-uz(i-1,j,k-1))/dxc(i-1) &
                                )
           !
           ! z component of the vorticity at cell center
           !
           voz(i,j,k) = 0.25_rp*( &
                                 (uy(i+1,j  ,k)-uy(i  ,j  ,k))/dxc(i  ) - (ux(i  ,j+1,k)-ux(i  ,j  ,k))/dyc(j  ) + &
                                 (uy(i+1,j-1,k)-uy(i  ,j-1,k))/dxc(i  ) - (ux(i  ,j  ,k)-ux(i  ,j-1,k))/dyc(j-1) + &
                                 (uy(i  ,j  ,k)-uy(i-1,j  ,k))/dxc(i-1) - (ux(i-1,j+1,k)-ux(i-1,j  ,k))/dyc(j  ) + &
                                 (uy(i  ,j-1,k)-uy(i-1,j-1,k))/dxc(i-1) - (ux(i-1,j  ,k)-ux(i-1,j-1,k))/dyc(j-1) &
                                )
        end do
      end do
    end do
    !$OMP END PARALLEL
  end subroutine cmpt_vorticity
  !
  subroutine cmpt_strain_rate(n,dxc,dyc,dzc,dxf,dyf,dzf,ux,uy,uz,str)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:)       :: dxc,dyc,dzc, &
                                                  dxf,dyf,dzf
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: str
    real(rp) :: s11,s22,s33,s12,s13,s23
    integer :: i,j,k
    !
    ! compute sijsij, where sij = (1/2)(du_i/dx_j + du_j/dx_i)
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(n,dxc,dyc,dzc,dxf,dyf,dzf,ux,uy,uz,str) &
    !$OMP PRIVATE(i,j,k,s11,s12,s13,s22,s23,s33)
    !$OMP DO
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          s11 = ((ux(i,j,k)-ux(i-1,j,k))/dxf(i))**2
          s22 = ((uy(i,j,k)-uy(i,j-1,k))/dyf(j))**2
          s33 = ((uz(i,j,k)-uz(i,j,k-1))/dzf(k))**2
          s12 = .25_rp*( &
                        ((ux(i  ,j+1,k)-ux(i  ,j  ,k))/dyc(j  ) + (uy(i+1,j  ,k)-uy(i  ,j  ,k))/dxc(i  ))**2 + &
                        ((ux(i  ,j  ,k)-ux(i  ,j-1,k))/dyc(j-1) + (uy(i+1,j-1,k)-uy(i  ,j-1,k))/dxc(i  ))**2 + &
                        ((ux(i-1,j+1,k)-ux(i-1,j  ,k))/dyc(j  ) + (uy(i  ,j  ,k)-uy(i-1,j  ,k))/dxc(i-1))**2 + &
                        ((ux(i-1,j  ,k)-ux(i-1,j-1,k))/dyc(j-1) + (uy(i  ,j-1,k)-uy(i-1,j-1,k))/dxc(i-1))**2 &
                       )*.25_rp
          s13 = .25_rp*( &
                        ((ux(i  ,j,k+1)-ux(i  ,j,k  ))/dzc(k  ) + (uz(i+1,j,k  )-uz(i  ,j,k  ))/dxc(i  ))**2 + &
                        ((ux(i  ,j,k  )-ux(i  ,j,k-1))/dzc(k-1) + (uz(i+1,j,k-1)-uz(i  ,j,k-1))/dxc(i  ))**2 + &
                        ((ux(i-1,j,k+1)-ux(i-1,j,k  ))/dzc(k  ) + (uz(i  ,j,k  )-uz(i-1,j,k  ))/dxc(i-1))**2 + &
                        ((ux(i-1,j,k  )-ux(i-1,j,k-1))/dzc(k-1) + (uz(i  ,j,k-1)-uz(i-1,j,k-1))/dxc(i-1))**2 &
                       )*.25_rp
          s23 = .25_rp*( &
                        ((uy(i,j  ,k+1)-uy(i,j  ,k  ))/dzc(k  ) + (uz(i,j+1,k  )-uz(i,j  ,k  ))/dyc(j  ))**2 + &
                        ((uy(i,j  ,k  )-uy(i,j  ,k-1))/dzc(k-1) + (uz(i,j+1,k-1)-uz(i,j  ,k-1))/dyc(j  ))**2 + &
                        ((uy(i,j-1,k+1)-uy(i,j-1,k  ))/dzc(k  ) + (uz(i,j  ,k  )-uz(i,j-1,k  ))/dyc(j-1))**2 + &
                        ((uy(i,j-1,k  )-uy(i,j-1,k-1))/dzc(k-1) + (uz(i,j  ,k-1)-uz(i,j-1,k-1))/dyc(j-1))**2 &
                       )*.25_rp
          str(i,j,k) = s11+s22+s33 + 2*(s12+s13+s23)
        end do
      end do
    end do
    !$OMP END PARALLEL
  end subroutine cmpt_strain_rate
  !
  subroutine cmpt_rotation_rate(n,dxc,dyc,dzc,ux,uy,uz,ens)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(0:)       :: dxc,dyc,dzc
    real(rp), intent(in ), dimension(0:,0:,0:) :: ux,uy,uz
    real(rp), intent(out), dimension(1:,1:,1:) :: ens
    real(rp) :: e12,e13,e23
    integer :: i,j,k
    !
    ! compute wijwij, where wij = (1/2)(du_i/dx_j - du_j/dx_i)
    !
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(n,dxc,dyc,dzc,ux,uy,uz,ens) &
    !$OMP PRIVATE(i,j,k,e12,e13,e23)
    !$OMP DO
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          e12 = .25_rp*( &
                        ((ux(i  ,j+1,k)-ux(i  ,j  ,k))/dyc(j  ) - (uy(i+1,j  ,k)-uy(i  ,j  ,k))/dxc(i  ))**2 + &
                        ((ux(i  ,j  ,k)-ux(i  ,j-1,k))/dyc(j-1) - (uy(i+1,j-1,k)-uy(i  ,j-1,k))/dxc(i  ))**2 + &
                        ((ux(i-1,j+1,k)-ux(i-1,j  ,k))/dyc(j  ) - (uy(i  ,j  ,k)-uy(i-1,j  ,k))/dxc(i-1))**2 + &
                        ((ux(i-1,j  ,k)-ux(i-1,j-1,k))/dyc(j-1) - (uy(i  ,j-1,k)-uy(i-1,j-1,k))/dxc(i-1))**2 &
                       )*.25_rp
          e13 = .25_rp*( &
                        ((ux(i  ,j,k+1)-ux(i  ,j,k  ))/dzc(k  ) - (uz(i+1,j,k  )-uz(i  ,j,k  ))/dxc(i  ))**2 + &
                        ((ux(i  ,j,k  )-ux(i  ,j,k-1))/dzc(k-1) - (uz(i+1,j,k-1)-uz(i  ,j,k-1))/dxc(i  ))**2 + &
                        ((ux(i-1,j,k+1)-ux(i-1,j,k  ))/dzc(k  ) - (uz(i  ,j,k  )-uz(i-1,j,k  ))/dxc(i-1))**2 + &
                        ((ux(i-1,j,k  )-ux(i-1,j,k-1))/dzc(k-1) - (uz(i  ,j,k-1)-uz(i-1,j,k-1))/dxc(i-1))**2 &
                       )*.25_rp
          e23 = .25_rp*( &
                        ((uy(i,j  ,k+1)-uy(i,j  ,k  ))/dzc(k  ) - (uz(i,j+1,k  )-uz(i,j  ,k  ))/dyc(j  ))**2 + &
                        ((uy(i,j  ,k  )-uy(i,j  ,k-1))/dzc(k-1) - (uz(i,j+1,k-1)-uz(i,j  ,k-1))/dyc(j  ))**2 + &
                        ((uy(i,j-1,k+1)-uy(i,j-1,k  ))/dzc(k  ) - (uz(i,j  ,k  )-uz(i,j-1,k  ))/dyc(j-1))**2 + &
                        ((uy(i,j-1,k  )-uy(i,j-1,k-1))/dzc(k-1) - (uz(i,j  ,k-1)-uz(i,j-1,k-1))/dyc(j-1))**2 &
                       )*.25_rp
          ens(i,j,k) =  2._rp*(e12+e13+e23)
        end do
      end do
    end do
    !$OMP END PARALLEL
  end subroutine cmpt_rotation_rate
  !
  subroutine cmpt_q_criterion(n,ens,str,qcr)
    implicit none
    integer , intent(in ), dimension(3)        :: n
    real(rp), intent(in ), dimension(1:,1:,1:) :: ens,str
    real(rp), intent(out), dimension(1:,1:,1:) :: qcr
    integer  :: i,j,k
    !$OMP PARALLEL DEFAULT(none) &
    !$OMP SHARED(n,ens,str,qcr) &
    !$OMP PRIVATE(i,j,k)
    !$OMP DO
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          qcr(i,j,k) = .5_rp*(ens(i,j,k)-str(i,j,k))
        end do
      end do
    end do
    !$OMP END PARALLEL
  end subroutine cmpt_q_criterion
  !
  subroutine cmpt_wall_forces(n,is_bound,dxc,dxf,dyc,dyf,dzc,dzf,visc,u,v,w,p,bforce, &
                              tau_x,tau_y,tau_z)
    use mpi_f08
    use mod_common_mpi, only: comm_block
    implicit none
    integer , intent(in ), dimension(3)        :: n
    logical , intent(in ), dimension(0:1,3)    :: is_bound
    real(rp), intent(in ), dimension(0:)       :: dxc,dyc,dzc, &
                                                  dxf,dyf,dzf
    real(rp), intent(in )                      :: visc
    real(rp), intent(in ), dimension(0:,0:,0:) :: u,v,w,p
    real(rp), intent(in ), dimension(1:3)      :: bforce
    real(rp), intent(out), dimension(0:1,3)    :: tau_x,tau_y,tau_z
    integer :: i,j,k,idir,ib
    real(rp) :: wall_shear,norm,sgn
    !
    tau_x(:,:) = 0.
    idir = 1
    do ib = 0,1
      if(ib == 0) then
        i    = 1
        norm =  1.
      endif
      if(ib == 1) then
        i    = n(idir) + 1
        norm = -1.
      endif
      if(is_bound(ib,idir)) then
        do k=1,n(3)
          do j=1,n(2)
            tau_x(ib,idir) = tau_x(ib,idir) + norm*( &
                                                    - 0.5*(p(i,j,k)+p(i-1,j,k)) + bforce(idir) &
                                                   )*dyf(j)*dzf(k)
          end do
        end do
      end if
    enddo
    idir = 2
    do ib = 0,1
      if(ib == 0) then
        j   = 1
        sgn =  1.
      endif
      if(ib == 1) then
        j   = n(idir) + 1
        sgn = -1.
      endif
      if(is_bound(ib,idir)) then
        do k=1,n(3)
          do i=1,n(1)
            wall_shear = sgn*(u(i,j,k)-u(i,j-1,k))/dyc(j-1)*visc
            tau_x(ib,idir) = tau_x(ib,idir) + wall_shear*dxc(i)*dzf(k)
          end do
        end do
      end if
    end do
    idir = 3
    do ib = 0,1
      if(ib == 0) then
        k   = 1
        sgn =  1.
      endif
      if(ib == 1) then
        k   = n(idir) + 1
        sgn = -1.
      endif
      if(is_bound(ib,idir)) then
        do j=1,n(2)
          do i=1,n(1)
            wall_shear = sgn*(u(i,j,k)-u(i,j,k-1))/dzc(k-1)*visc
            tau_x(ib,idir) = tau_x(ib,idir) + wall_shear*dxc(i)*dyf(j)
          end do
        end do
      end if
    end do
    call mpi_allreduce(MPI_IN_PLACE,tau_x,2*3,MPI_REAL_RP,MPI_SUM,comm_block)
    !
    tau_y(:,:) = 0.
    idir = 2
    do ib = 0,1
      if(ib == 0) then
        j    = 1
        norm =  1.
      endif
      if(ib == 1) then
        j    = n(idir) + 1
        norm = -1.
      endif
      if(is_bound(ib,idir)) then
        do k=1,n(3)
          do i=1,n(1)
            tau_y(ib,idir) = tau_y(ib,idir) + norm*( &
                                                    - 0.5*(p(i,j,k)+p(i,j-1,k)) + bforce(idir) &
                                                   )*dxf(i)*dzf(k)
          end do
        end do
      end if
    enddo
    idir = 3
    do ib = 0,1
      if(ib == 0) then
        k   = 1
        sgn =  1.
      endif
      if(ib == 1) then
        k   = n(idir) + 1
        sgn = -1.
      endif
      if(is_bound(ib,idir)) then
        do j=1,n(2)
          do i=1,n(1)
            wall_shear = sgn*(v(i,j,k)-v(i,j,k-1))/dzc(k-1)*visc
            tau_y(ib,idir) = tau_y(ib,idir) + wall_shear*dyc(j)*dxf(i)
          end do
        end do
      end if
    end do
    idir = 1
    do ib = 0,1
      if(ib == 0) then
        i   = 1
        sgn =  1.
      endif
      if(ib == 1) then
        i   = n(idir) + 1
        sgn = -1.
      endif
      if(is_bound(ib,idir)) then
        do k=1,n(3)
          do j=1,n(2)
            wall_shear = sgn*(v(i,j,k)-v(i-1,j,k))/dxc(i-1)*visc
            tau_y(ib,idir) = tau_y(ib,idir) + wall_shear*dyc(j)*dzf(k)
          end do
        end do
      end if
    end do
    call mpi_allreduce(MPI_IN_PLACE,tau_y,2*3,MPI_REAL_RP,MPI_SUM,comm_block)
    !
    tau_z(:,:) = 0.
    idir = 3
    do ib = 0,1
      if(ib == 0) then
        k    = 1
        norm =  1.
      endif
      if(ib == 1) then
        k    = n(idir) + 1
        norm = -1.
      endif
      if(is_bound(ib,idir)) then
        do j=1,n(2)
          do i=1,n(1)
            tau_z(ib,idir) = tau_z(ib,idir) + norm*( &
                                                    - 0.5*(p(i,j,k)+p(i,j,k-1)) + bforce(idir) &
                                                   )*dxf(i)*dyf(j)
          end do
        end do
      end if
    enddo
    idir = 1
    do ib = 0,1
      if(ib == 0) then
        i   = 1
        sgn =  1.
      endif
      if(ib == 1) then
        i   = n(idir) + 1
        sgn = -1.
      endif
      if(is_bound(ib,idir)) then
        do k=1,n(3)
          do j=1,n(2)
            wall_shear = sgn*(w(i,j,k)-w(i-1,j,k))/dxc(i-1)*visc
            tau_z(ib,idir) = tau_z(ib,idir) + wall_shear*dzc(k)*dyf(j)
          end do
        end do
      end if
    end do
    idir = 2
    do ib = 0,1
      if(ib == 0) then
        j   = 1
        sgn =  1.
      endif
      if(ib == 1) then
        j   = n(idir) + 1
        sgn = -1.
      endif
      if(is_bound(ib,idir)) then
        do k=1,n(3)
          do i=1,n(1)
            wall_shear = sgn*(w(i,j,k)-w(i,j-1,k))/dyc(j-1)*visc
            tau_z(ib,idir) = tau_z(ib,idir) + wall_shear*dzc(k)*dxf(i)
          end do
        end do
      end if
    end do
    call mpi_allreduce(MPI_IN_PLACE,tau_z,2*3,MPI_REAL_RP,MPI_SUM,comm_block)
  end subroutine cmpt_wall_forces
  !
  subroutine updt_wall_forces(rkpar,tau_x,tau_y,tau_z,tau_x_o,tau_y_o,tau_z_o, &
                              tau_x_acc,tau_y_acc,tau_z_acc)
    implicit none
    real(rp), intent(in   ), dimension(2) :: rkpar
    real(rp), intent(in   ), dimension(0:1,3) :: tau_x,tau_y,tau_z
    real(rp), intent(inout), dimension(0:1,3) :: tau_x_o,tau_y_o,tau_z_o
    real(rp), intent(inout), dimension(0:1,3) :: tau_x_acc,tau_y_acc,tau_z_acc
    integer :: il,iu,is,idir
    !
    ! increment in time shear and normal stresses, consistently with
    ! the RK3 time integration scheme
    !
    ! x momentum terms
    !
    idir = 1; il = 2; iu = 3; is = 1
    tau_x_acc(:,il:iu:is) = tau_x_acc(:,il:iu:is) + rkpar(1)*tau_x(:,il:iu:is) + rkpar(2)*tau_x_o(:,il:iu:is)
    tau_x_acc(:,idir    ) = tau_x_acc(:,idir    ) + (rkpar(1)+rkpar(2))*tau_x(:,idir)
    !
    ! y momentum terms
    !
    idir = 2; il = 1; iu = 3; is = 2
    tau_y_acc(:,il:iu:is) = tau_y_acc(:,il:iu:is) + rkpar(1)*tau_y(:,il:iu:is) + rkpar(2)*tau_y_o(:,il:iu:is)
    tau_y_acc(:,idir    ) = tau_y_acc(:,idir    ) + (rkpar(1)+rkpar(2))*tau_y(:,idir)
    !
    ! z momentum terms
    !
    idir = 3; il = 1; iu = 2; is = 1
    tau_z_acc(:,il:iu:is) = tau_z_acc(:,il:iu:is) + rkpar(1)*tau_z(:,il:iu:is) + rkpar(2)*tau_z_o(:,il:iu:is)
    tau_z_acc(:,idir    ) = tau_z_acc(:,idir    ) + (rkpar(1)+rkpar(2))*tau_z(:,idir)
    !
    ! update previous time step values
    !
    tau_x_o(:,:) = tau_x(:,:)
    tau_y_o(:,:) = tau_y(:,:)
    tau_z_o(:,:) = tau_z(:,:)
  end subroutine updt_wall_forces
end module mod_post
