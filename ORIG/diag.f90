 subroutine diag
 
 use shared

 real(kind=8)    :: Ray, Gras, Four, Pr, beta_gras, &
                    Delta_Temp,vert_flux_t(0:Ny+1),horiz_flux_t(0:Nx+1), &
                    Nus, Fou, theta, Nus_comp, grad_T, Nus_comp_2, melt, &
                    melt_front, melt_front_expected, theta1, theta2
    
   

 integer(kind=8) :: i_liq, i_sol, Ny_2, i_fix

 Ny_2  = int(Ny/2)
 i_liq = int(Nx/20)
 i_sol = Nx - int(Nx/20)

 Delta_Temp = T0 - T_in_Sol  !  T(i_liq,Ny_2) - T(i_sol,Ny_2)

 if(phi(i_liq,Ny_2).lt.1.0) then
   i_liq = 2
 endif
 if(phi(i_sol,Ny_2).gt.0.0) then
   i_sol = Nx-1
 endif

 ! First, let's compute the Grashoph Number in the two phases

! beta_gras = (1-rho(i_liq,Ny_2)/rho(i_sol,Ny_2)) / &
 beta_gras = 0.0d0
 beta_gras = (rho0/rho1 - 1) / Delta_Temp   !(1-rho1/rho0) / Delta_Temp  

 Gras = grav * beta_gras * Delta_Temp * Ny**3 / kvisc**2
 
 Pr   = (2.*tau_f -1.) / (2.*tau_g - 1.)

 Ray  = Pr * Gras

 Fou  = alpha*it / (Ny**2)

 theta= stefan_number * Fou

 theta1 = Ray**(-0.5d0)
 theta2 = (Nx/Ny) * Ray**(-1./4.) 

  write(*,*) 'Diagnostics CHECK:'
  write(*,*) 'beta_gras, grav, Gras, Delta_Temp, Pr, Ray, Fou, theta'
  write(*,*) beta_gras, grav, Gras, Delta_Temp, Pr, Ray, Fou, theta
  write(*,*) '---------------------------------------------'

! if(theta.le.theta2) then
    Nus  = ((2*theta)**(-0.5d0)) + &
           (0.35d0 * (Ray**(0.25d0)) - ((2*theta)**(-0.5d0))) * &
           ( 1. + (0.0175d0 * (Ray**(3./4.))*(theta**(3./2.)))**(-2.d0))**(1./(-2.))
! else
!    Nus  = ((1. - (Ray**(1./4.))*(theta - theta2)))**(3./5.)
! endif

Nus_comp = 0.0d0

do j=1,Ny
 Nus_comp = Nus_comp + Ny * ((Tin - T(2,j))/2.)/Delta_Temp
enddo

 Nus_comp = Nus_comp / Ny

 Nus_comp_2 = 0.0d0

 do j=1,Ny
    grad_T = 0.0d0
    grad_T = abs((- 3.*T(1,j)+4*T(2,j)-T(3,j))/2.)  ! One-sided, 2nd Order F.D. for Temperature

    Nus_comp_2 = Nus_comp_2 + Ny * grad_T/Delta_Temp
 enddo
 
    Nus_comp_2 = Nus_comp_2 / Ny


! write(*,*) 'Gras =',Gras,'Ray = ',Ray
 if(mod(it,100).eq.0) then

    melt = 0.0d0
    
     do j=1,Ny
      do i=1,Nx
    
       if(phi(i,j).gt.0) then
         melt = melt + 1
       endif
    
      enddo
     enddo
    
     melt = melt/(Nx*Ny)

     if(theta.le.theta1) then
       melt_front = Ny/((2*theta)**(-0.5d0))
       melt_front_expected = Ny / Nus
     elseif(theta.gt.theta1.and.theta.lt.theta2) then
       melt_front=0.0d0
       do j = 1,Ny
        do i = 1,Nx
           if(phi(i,j).eq.0)  then 
             melt_front = melt_front + i
           endif 
        enddo
       enddo
       melt_front = melt_front / Ny
       melt_front_expected = Ny * theta * Ray**(1./4.)
     endif
 
   write(*,*) ' #@#@#@   THETA CHECK :   #@#@#@'
   write(*,*) 'theta=',theta, 'theta1=',theta1, 'theta2=',theta2
   write(164,640) it, theta, Nus, Nus_comp, Nus_comp_2, melt, melt_front, &
                  melt_front_expected, Gras, Pr, Ray, Fou
 endif

640 format(1I8,11(1x,e13.6))

 if(mod(it,500).eq.0) then
 do i=1,Nx
   write(165,*) i, T(i,int(Ny/2))
 enddo
 write(165,'(bn)')
 write(165,'(bn)')
 endif


 do i=0,Nx+1
   if(phi(i,int(0.98*Ny)).eq.0) then
     i_liq = i
   endif
 enddo

! write(*,*) 'i_lq = ', i_liq
!!!
!!! do j=int(Ny/2)+1,Ny+1
!!!   vert_flux_t(j) = v(i_liq,j) * T(i_liq,j) + alpha*(T(i_liq,j-1)-T(i_liq,j)) 
!!! enddo
!!! do j=0,int(Ny/2)
!!!   vert_flux_t(j) = v(i_liq,j) * T(i_liq,j) + alpha*(T(i_liq,j+1)-T(i_liq,j))  
!!! enddo
!!!
!!! do j=0,Ny+1
!!!   write(468,*) i_liq,j,vert_flux_t(j)
!!! enddo
!!!   write(468,'(bn)')
!!!   write(468,'(bn)')
!!!


 end subroutine diag
