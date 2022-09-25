 subroutine BC_g

 use shared

 !======================================== WEST

 do j = 1,Ny
    uu = u(0,j)*u(0,j)
    vv = v(0,j)*v(0,j)
    uv = 1.5d0*(uu+vv)
    upv = u(0,j)+v(0,j)
    umv = u(0,j)-v(0,j)
    g(1,0,j) = w_eq(1)*Tin*(1.d0+3.d0*u(0,j)+4.5d0*uu-uv)
    g(5,0,j) = w_eq(5)*Tin*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
    g(8,0,j) = w_eq(8)*Tin*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    T(0,j) = Tin

!!!    g(2,0,j) = T0*(w_eq(2)+w_eq(4))-g(4,0,j)
!!!    g(6,0,j) = T0*(w_eq(6)+w_eq(8))-g(8,0,j)
!!!    g(9,0,j) = T0*(w_eq(9)+w_eq(7))-g(7,0,j)

!    g(2,0,j) = (g(2,0,j)+u(1,j)*g(2,1,j))/(1.d0+u(1,j))
!    g(6,0,j) = (g(6,0,j)+u(1,j)*g(6,1,j))/(1.d0+u(1,j))
!    g(9,0,j) = (g(9,0,j)+u(1,j)*g(9,1,j))/(1.d0+u(1,j))
!    T(0,j) = T(1,j)

 end do

 !======================================== EAST

 do j = 1,Ny
    uu = u(Nx+1,j)*u(Nx+1,j)
    vv = v(Nx+1,j)*v(Nx+1,j)
    uv = 1.5d0*(uu+vv)
    upv = u(Nx+1,j)+v(Nx+1,j)
    umv = u(Nx+1,j)-v(Nx+1,j)
    g(3,Nx+1,j) = w_eq(3)*T(Nx,j)*(1.d0-3.d0*u(Nx+1,j)+4.5d0*uu-uv)
    g(6,Nx+1,j) = w_eq(6)*T(Nx,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
    g(7,Nx+1,j) = w_eq(7)*T(Nx,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
    T(Nx+1,j) = T(Nx,j)

!    g(4,Nx+1,j) = (g(4,Nx+1,j)+u(Nx,j)*g(4,Nx,j))/(1.d0+u(Nx,j))
!    g(7,Nx+1,j) = (g(7,Nx+1,j)+u(Nx,j)*g(7,Nx,j))/(1.d0+u(Nx,j))
!    g(8,Nx+1,j) = (g(8,Nx+1,j)+u(Nx,j)*g(8,Nx,j))/(1.d0+u(Nx,j))
!    T(Nx+1,j) = T(Nx,j)

 end do

 !======================================== NORTH

!!! do i = 0,Nx+1
!!!    T(i,Ny+1)   = T(i,Ny) !!T0
!!!    g(4,i,Ny+1) = w_eq(4)*T(i,Ny+1)
!!!    g(7,i,Ny+1) = w_eq(7)*T(i,Ny+1)
!!!    g(8,i,Ny+1) = w_eq(8)*T(i,Ny+1)
!!!!    T(i,Ny+1) = T(i,Ny) !T0
!!!    phi(i,Ny+1) = phi(i,Ny)
!!! end do



!!!! do i = 0,int(Nx/2)-1
!!!!    T(i,Ny+1)   = T(i,Ny) !!T0
!!!!    g(4,i,Ny+1) = w_eq(4)*T(i,Ny+1)
!!!!    g(7,i,Ny+1) = w_eq(7)*T(i,Ny+1)
!!!!    g(8,i,Ny+1) = w_eq(8)*T(i,Ny+1)
!!!!!    T(i,Ny+1) = T(i,Ny) !T0
!!!!    phi(i,Ny+1) = phi(i,Ny)
!!!! end do
!!!! do i = int(Nx/2),Nx+1
!!!!    T(i,Ny+1)   = T(i,Ny) !T_in_sol
!!!!    g(4,i,Ny+1) = w_eq(4)*T(i,Ny+1)
!!!!    g(7,i,Ny+1) = w_eq(7)*T(i,Ny+1)
!!!!    g(8,i,Ny+1) = w_eq(8)*T(i,Ny+1)
!!!!!    T(i,Ny+1) = T(i,Ny) !T0
!!!!    phi(i,Ny+1) = phi(i,Ny)
!!!! end do

!!!!  do i=1,Nx     ! PURE "HALF-WAY NO-SLIP" for the TEMPERATURE
!!!!     g(4,i,Ny+1) = g(2,i  ,Ny)
!!!!     g(7,i,Ny+1) = g(5,i-1,Ny)
!!!!     g(8,i,Ny+1) = g(6,i+1,Ny)
!!!!!
!!!!     T(i,Ny+1)   = T(i,Ny)
!!!!     phi(i,Ny+1) = phi(i,Ny)
!!!!  enddo


!!!!  do i=1,Nx     ! PURE "FREE-SLIP" for the TEMPERATURE
!!!!     g(4,i,Ny+1) = g(2,i ,Ny)
!!!!     g(7,i,Ny+1) = g(6,i ,Ny)
!!!!     g(8,i,Ny+1) = g(5,i ,Ny)
!!!!!
!!!!     T(i,Ny+1)   = T(i,Ny)
!!!!     phi(i,Ny+1) = phi(i,Ny)
!!!!  enddo


  do i=0,Nx+1     ! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
   do k=0,npop-1
     g(k,i,Ny+1) = g(k,i ,Ny)
   enddo
     T(i,Ny+1)   = T(i,Ny)
     phi(i,Ny+1) = phi(i,Ny)
  enddo

 !======================================== SOUTH

!!! do i = 0,Nx+1
!!!    T(i,0) = T(i,1) !!T0
!!!    g(2,i,0) = w_eq(2)*T(i,0)
!!!    g(5,i,0) = w_eq(5)*T(i,0)
!!!    g(6,i,0) = w_eq(6)*T(i,0)
!!!!    T(i,0) = T(i,1)  !T0
!!!    phi(i,0) = phi(i,1)
!!! end do

!!! do i = 0,int(Nx/2)-1
!!!    T(i,0) = T(i,1) !!T0
!!!    g(2,i,0) = w_eq(2)*T(i,0)
!!!    g(5,i,0) = w_eq(5)*T(i,0)
!!!    g(6,i,0) = w_eq(6)*T(i,0)
!!!!    T(i,0) = T(i,1)  !T0
!!!    phi(i,0) = phi(i,1)
!!! end do
!!!!
!!! do i = int(Nx/2),Nx+1
!!!    T(i,0) = T(i,1) ! T_in_sol
!!!    g(2,i,0) = w_eq(2)*T(i,0)
!!!    g(5,i,0) = w_eq(5)*T(i,0)
!!!    g(6,i,0) = w_eq(6)*T(i,0)
!!!!    T(i,0) = T(i,1)  !T0
!!!    phi(i,0) = phi(i,1)
!!! end do

!!!! NO-SLIP 
!!!    do i = 1, Nx
!!!      g(2,i,0) = g(4,i  ,1)
!!!      g(5,i,0) = g(7,i+1,1)
!!!      g(6,i,0) = g(8,i-1,1)
!!!!
!!!      T(i,0) = T(i,1)
!!!      phi(i,0) = phi(i,1)
!!!    enddo
!!!
!!!!! FREE SLIP
!!!!    do i = 1, Nx
!!!!      g(2,i,0) = g(4,i ,1)
!!!!      g(5,i,0) = g(8,\i ,1)
!!!!      g(6,i,0) = g(7,i ,1)
!!!!!
!!!!      T(i,0) = T(i,1)
!!!!      phi(i,0) = phi(i,1)
!!!!    enddo

  do i=0,Nx+1     ! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
   do k=0,npop-1
     g(k,i,0) = g(k,i ,1)
   enddo
     T(i,0)   = T(i,1)
     phi(i,0) = phi(i,1)
  enddo


!!!
!!!
!!!!!!!VERTEXES  -> Helf-Way No-Slip, as well
!!!!!!!
!!!    g(5,0,0)       = g(7,1,1)
!!!    g(8,0,Ny+1)    = g(6,1,Ny)
!!!    g(6,Nx+1,0)    = g(8,Nx,1)
!!!    g(7,Nx+1,Ny+1) = g(5,Nx,Ny)
!!!!!! 

 end subroutine BC_g
