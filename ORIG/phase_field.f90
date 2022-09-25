 subroutine phase_field
 
 use shared
! real(kind=dp) :: react_plus(1:Nx,1:Ny),react_minus(1:Nx,1:Ny), & 
!                  k_plus(1:Nx,1:Ny),k_minus(1:Nx,1:Ny)
!                  provo a spostarle nel COMMON, per evitare casini con OpenMP

T_act  = 0.0d0
T_crit = 273.0

f_plus  = 0.5d0
f_minus = 0.5d0 

! Compute Phi Selector (Succi, PRL 2001)

phi_prec(:,:)    = phi(:,:)
k_plus(:,:)      = 0.0d0
k_minus(:,:)     = 0.0d0
react_plus(:,:)  = 0.0d0
react_minus(:,:) = 0.0d0
react_tot(:,:)   = 0.0d0

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j,k_plus,k_minus)
!!!!$OMP PRIVATE(i,j,react_plus,react_minus,k_plus,k_minus,react_tot,phi)
!!!!$OMP SECTIONS
!!!!$OMP SECTION
!$OMP DO
 do j = 1,Ny
  do i = 1,Nx

   react_plus(i,j)  = 0.5d0*(1.+tanh(((T(i,j)-T_crit)-T_act)/T_width))
   react_minus(i,j) = 0.5d0*(1.-tanh(((T(i,j)-T_crit)+T_act)/T_width))

  enddo
 enddo 
!$OMP END DO
k_plus(:,:)  = 0.0d0
k_minus(:,:) = 0.0d0
!!!!$OMP SECTION
!$OMP DO
 do j = 1,Ny
  do i = 1,Nx
   do k=0,npop-1
       k_plus(i,j)  = k_plus(i,j) + (0.5d0*(1.+phi(i+cx(k),j+cy(k))))**2.0d0
       k_minus(i,j) = k_minus(i,j)+ (0.5d0*(1.-phi(i+cx(k),j+cy(k))))**2.0d0
   enddo
       k_plus(i,j)  = (k_plus(i,j)  / 8.0d0) * react_plus(i,j)
       k_minus(i,j) = (k_minus(i,j) / 8.0d0) * react_minus(i,j)
  enddo
 enddo
!$OMP END DO
!!!!$OMP SECTION
!$OMP DO
 do j=1,Ny
   do i=1,Nx
    react_tot (i,j) = f_plus * k_plus(i,j)  * (1. - phi_prec(i,j)) - &
                      f_minus* k_minus(i,j) * (1. + phi_prec(i,j)) 
  
    phi (i,j)       = phi_prec(i,j) + react_tot(i,j)

   enddo
 enddo
!OMP END DO
!!!!!!!!!!$OMP END SECTIONS
!$OMP END PARALLEL

! write(*,*) 'f_plus=',f_plus,'f_minus=',f_minus

 end subroutine phase_field
