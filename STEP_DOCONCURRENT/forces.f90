 subroutine forces 
 
 use shared

 real(kind=mykind)  :: forceX, forceY

! First: Update the Temperature via the latent heat release/absorption
! Second: compute BUOYANCY FORCE
! Third: Compute SOLID DRAG  (hic sunt leones....)
! Fourth: Update Velocity Field for equilibrum and then collision:

 delta_T = T0 - T_in_sol   !T_width*2.d0
 beta = (c10 - rho1/rho0) / delta_T

 do concurrent(j=1:Ny,i=1:Nx)
!        
!1
        T(i,j)=T(i,j)-(delta_T/stefan_number*(phi(i,j)-phi_prec(i,j))*c05) 
!        
!2        
        forceY = rho(i,j)*grav*beta*T(i,j)*(c10+phi(i,j))*c05
!        
!3         
        forceX =              - c10*u(i,j)*(c10-phi(i,j))*c05
        forceY = forceY       - c10*v(i,j)*(c10-phi(i,j))*c05
!        
!4      
        u(i,j)=u(i,j)+forceX/(omega_f*rho(i,j))
        v(i,j)=v(i,j)+forceY/(omega_f*rho(i,j))
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
