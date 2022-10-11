 subroutine forces 
 
 use shared


! First: Update the Temperature via the latent heat release/absorption


!    do concurrent(j=1:Ny,i=1:Nx)
!       T_prec(i,j) = T(i,j)
!    enddo

 delta_T = T0 - T_in_sol   !T_width*2.d0

! beta = 3.5e-7
 beta = (1 - rho1/rho0) / delta_T

    do concurrent(j=1:Ny,i=1:Nx)
       T(i,j)=T(i,j)-(delta_T/stefan_number*(phi(i,j)-phi_prec(i,j))/2.0d0) !*omega_g
!       T(i,j)=T_prec(i,j)-(delta_T/stefan_number*(phi(i,j)-phi_prec(i,j))/2.0d0) !*omega_g
    enddo
!    
! Second: compute BUOYANCY FORCE
    do concurrent(j=1:Ny,i=1:Nx)
      force_y(i,j) = rho(i,j) * grav * beta*T(i,j) * (1.0d0+phi(i,j)/2.0d0)
    enddo
!
! Third: Compute SOLID DRAG  (hic sunt leones....)
    do concurrent(j=1:Ny,i=1:Nx)
      force_x(i,j) =              - 1.0d0*u(i,j)*(1.d0-phi(i,j))*0.50d0
      force_y(i,j) = force_y(i,j) - 1.0d0*v(i,j)*(1.d0-phi(i,j))*0.50d0
    enddo

! Update Velocity Field for equilibrum and then collision:
    do concurrent(j=1:Ny,i=1:Nx)
       u(i,j)=u(i,j)+force_x(i,j)/(omega_f*rho(i,j))
       v(i,j)=v(i,j)+force_y(i,j)/(omega_f*rho(i,j))
    enddo

#ifdef DEBUG_GA
    do i=1,Nx
      if(phi(i,int(Ny/2)).ge.0.0d0.and.phi(i+1,int(Ny/2)).lt.0.0d0) then
        location = float((i+(i+1))/2)
        write(128,*) it, it*alpha/(Nx*Nx), location
      endif
    enddo
#endif

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine forces"
#endif

 end subroutine forces
