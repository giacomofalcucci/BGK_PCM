 subroutine forces 
 
 use shared

! force_x(:,:) = 0.0d0
! force_y(:,:) = 0.0d0

! First: Update the Temperature via the latent heat release/absorption

! T_prec(:,:) = T(:,:)

 delta_T = T0 - T_in_sol   !T_width*2.d0

! beta = 3.5e-7
 beta = (1 - rho1/rho0) / delta_T

!GA1 !$OMP PARALLEL DEFAULT(NONE)  &
!GA1 !$OMP PRIVATE(i,j)            &
!GA1 !$OMP SHARED(Nx,Ny)           &
!GA1 !$OMP SHARED(delta_T,stefan_number)  &
!GA1 !$OMP SHARED(grav,beta,omega_f)  &
!GA1 !$OMP SHARED(u,v,T,rho)  &
!GA1 !$OMP SHARED(phi,phi_prec)  &
!GA1 !$OMP SHARED(force_y,force_x) 
!GA1 !$OMP DO
!$OMP target teams distribute parallel do collapse(2)
    do j=1,Ny
     do i=1,Nx
         T(i,j)=T(i,j)-(delta_T/stefan_number*(phi(i,j)-phi_prec(i,j))/2.0d0) !*omega_g
     enddo
    enddo


! Second: compute BUOYANCY FORCE
!GA1 !$OMP DO
!$OMP target teams distribute parallel do collapse(2)
    do j=1,Ny
     do i=1,Nx
      force_y(i,j) = rho(i,j) * grav * beta*T(i,j) * (1.0d0+phi(i,j)/2.0d0)
     enddo
    enddo

! Third: Compute SOLID DRAG  (hic sunt leones....)
!GA1 !$OMP DO
!$OMP target teams distribute parallel do collapse(2)
    do j=1,Ny
     do i=1,Nx
      force_x(i,j) =              - 1.0d0*u(i,j)*(1.d0-phi(i,j))*0.50d0
      force_y(i,j) = force_y(i,j) - 1.0d0*v(i,j)*(1.d0-phi(i,j))*0.50d0
     enddo
    enddo

! Update Velocity Field for equilibrum and then collision:
!GA1 !$OMP DO
!$OMP target teams distribute parallel do collapse(2)
    do j=1,Ny
     do i=1,Nx
           u(i,j)=u(i,j)+force_x(i,j)/(omega_f*rho(i,j))
           v(i,j)=v(i,j)+force_y(i,j)/(omega_f*rho(i,j))
     enddo
    enddo
!GA1 !$OMP END PARALLEL 

#ifdef DEBUG_GA    
    do i=1,Nx
      if(phi(i,int(Ny/2)).ge.0.0d0.and.phi(i+1,int(Ny/2)).lt.0.0d0) then
        location = float((i+(i+1))/2)
!!        write(*,*) 'phi = 0, intrface located at ', location, int(Ny/2)
        write(128,*) it, it*alpha/(Nx*Nx), location
      endif
    enddo
#endif


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine forces"
#endif

 end subroutine forces
