 subroutine forces 
 
 use shared

force_x(:,:) = 0.0d0
force_y(:,:) = 0.0d0

! First: Update the Temperature via the latent heat release/absorption

 T_prec(:,:) = T(:,:)

 delta_T = T0 - T_in_sol   !T_width*2.d0

! beta = 3.5e-7
 beta = (1 - rho1/rho0) / delta_T

!!!#!!$OMP PARALLEL DEFAULT(shared) &
!!!#!!$OMP PRIVATE(i,j)
!!!#!!$OMP DO

    do j=1,Ny
     do i=1,Nx
         T(i,j)=T_prec(i,j)-(delta_T/stefan_number*(phi(i,j)-phi_prec(i,j))/2.0d0) !*omega_g
     enddo
    enddo
!!!#!!$OMP END DO
!!!#!!
!!!#!!$OMP DO
! Second: compute BUOYANCY FORCE
    do j=1,Ny
     do i=1,Nx
!      force_y(i,j) = rho(i,j) * grav* (-1.0d0  + beta*T(i,j)*(1.0d0+phi(i,j)))   !/2.0d0)
      force_y(i,j) = rho(i,j) * grav * beta*T(i,j) * (1.0d0+phi(i,j)/2.0d0)
     enddo
    enddo
!!!#!!$OMP END DO
!!!#!!
!!!#!!$OMP DO
! Third: Compute SOLID DRAG  (hic sunt leones....)
    do j=1,Ny
     do i=1,Nx
      force_x(i,j) = force_x(i,j) - 1.0d0*u(i,j)*(1.d0-phi(i,j))*0.50d0
      force_y(i,j) = force_y(i,j) - 1.0d0*v(i,j)*(1.d0-phi(i,j))*0.50d0
     enddo
    enddo
!!!#!!$OMP END DO
!!!#!!
!!!#!!$OMP DO 
! Update Velocity Field for equilibrum and then collision:
    do j=1,Ny
     do i=1,Nx
           u(i,j)=u(i,j)+force_x(i,j)/(omega_f*rho(i,j))
           v(i,j)=v(i,j)+force_y(i,j)/(omega_f*rho(i,j))
     enddo
    enddo
!!!#!!$OMP END DO
!!!#!!
!!!#!!$OMP DO
!!!! do j=1,Ny
!!!!  do i=1,Nx
!!!!    write(128,*) i,j,force_x(i,j),force_y(i,j)
!!!!  enddo
!!!!    write(128,'(bn)')
!!!! enddo
!!!!    write(128,'(bn)')
!!!!    write(128,'(bn)')


    do i=1,Nx
      if(phi(i,int(Ny/2)).ge.0.0d0.and.phi(i+1,int(Ny/2)).lt.0.0d0) then
        location = float((i+(i+1))/2)
!!        write(*,*) 'phi = 0, intrface located at ', location, int(Ny/2)
        write(128,*) it, it*alpha/(Nx*Nx), location
      endif
    enddo

!!!#!!$OMP END DO
!!!#!!$OMP END PARALLEL
!!

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine forces"
#endif

 end subroutine forces
