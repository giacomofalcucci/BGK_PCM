 subroutine phase_field
 
 use shared
! real(kind=mykind) :: react_plus(1:Nx,1:Ny),react_minus(1:Nx,1:Ny), & 
!                  k_plus(1:Nx,1:Ny),k_minus(1:Nx,1:Ny)
!                  provo a spostarle nel COMMON, per evitare casini con OpenMP

real(kind=mykind) :: react_plus,react_minus 
real(kind=mykind) :: k_plus,k_minus

T_act  = 0.0d0
T_crit = 273.0

f_plus  = 0.5d0
f_minus = 0.5d0 

! Compute Phi Selector (Succi, PRL 2001)

phi_prec(:,:)    = phi(:,:)

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j,k_plus,k_minus)
!!!!$OMP PRIVATE(i,j,react_plus,react_minus,k_plus,k_minus,react_tot,phi)
!!!!$OMP SECTIONS
!!!!$OMP SECTION
!$OMP DO
 do j = 1,Ny
  do i = 1,Nx

     react_plus  = 0.5d0*(1.+tanh(((T(i,j)-T_crit)-T_act)/T_width))
     react_minus = 0.5d0*(1.-tanh(((T(i,j)-T_crit)+T_act)/T_width))

     k_plus = + (0.5d0*(1.+phi(i  ,j  )))**2.0d0 &
                   + (0.5d0*(1.+phi(i+1,j  )))**2.0d0 & 
                   + (0.5d0*(1.+phi(i  ,j+1)))**2.0d0 & 
                   + (0.5d0*(1.+phi(i-1,j  )))**2.0d0 & 
                   + (0.5d0*(1.+phi(i  ,j-1)))**2.0d0 & 
                   + (0.5d0*(1.+phi(i+1,j+1)))**2.0d0 & 
                   + (0.5d0*(1.+phi(i-1,j+1)))**2.0d0 & 
                   + (0.5d0*(1.+phi(i-1,j-1)))**2.0d0 & 
                   + (0.5d0*(1.+phi(i+1,j-1)))**2.0d0   

     k_minus =+ (0.5d0*(1.-phi(i  ,j  )))**2.0d0 &
                   + (0.5d0*(1.-phi(i+1,j  )))**2.0d0 &
                   + (0.5d0*(1.-phi(i  ,j+1)))**2.0d0 &
                   + (0.5d0*(1.-phi(i-1,j  )))**2.0d0 &
                   + (0.5d0*(1.-phi(i  ,j-1)))**2.0d0 &
                   + (0.5d0*(1.-phi(i+1,j+1)))**2.0d0 &
                   + (0.5d0*(1.-phi(i-1,j+1)))**2.0d0 &
                   + (0.5d0*(1.-phi(i-1,j-1)))**2.0d0 &
                   + (0.5d0*(1.-phi(i+1,j-1)))**2.0d0  

      k_plus  = (k_plus  / 8.0d0) * react_plus
      k_minus = (k_minus / 8.0d0) * react_minus

    phi(i,j)       = phi_prec(i,j)  + (  & 
                     +  f_plus * k_plus  * (1. - phi_prec(i,j))  &
                     -  f_minus* k_minus * (1. + phi_prec(i,j)) )

   enddo
 enddo
!OMP END DO
!!!!!!!!!!$OMP END SECTIONS
!$OMP END PARALLEL

! write(*,*) 'f_plus=',f_plus,'f_minus=',f_minus

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine phase_field"
#endif

 end subroutine phase_field
