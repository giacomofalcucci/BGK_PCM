 subroutine BC_g

 use shared

 !======================================== WEST

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,j)            &
!$OMP PRIVATE(uu,vv,uv,upv,umv)  &
!$OMP SHARED(Nx,Ny,Tin)          &
!$OMP SHARED(w_eq1,w_eq2,w_eq3,w_eq4)     &
!$OMP SHARED(w_eq5,w_eq6,w_eq7,w_eq8)     &
!$OMP SHARED(g0,g1,g2,g3,g4)  &
!$OMP SHARED(g5,g6,g7,g8)  &
!$OMP SHARED(u,v,T,phi)
!$OMP DO
 do j = 1,Ny
    uu = u(0,j)*u(0,j)
    vv = v(0,j)*v(0,j)
    uv = 1.5d0*(uu+vv)
    upv = u(0,j)+v(0,j)
    umv = u(0,j)-v(0,j)
    g1(0,j) = w_eq1*Tin*(1.d0+3.d0*u(0,j)+4.5d0*uu-uv)
    g5(0,j) = w_eq5*Tin*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
    g8(0,j) = w_eq8*Tin*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    T(0,j) = Tin
 end do

 !======================================== EAST

!$OMP DO
 do j = 1,Ny
    uu = u(Nx+1,j)*u(Nx+1,j)
    vv = v(Nx+1,j)*v(Nx+1,j)
    uv = 1.5d0*(uu+vv)
    upv = u(Nx+1,j)+v(Nx+1,j)
    umv = u(Nx+1,j)-v(Nx+1,j)
    g3(Nx+1,j) = w_eq3*T(Nx,j)*(1.d0-3.d0*u(Nx+1,j)+4.5d0*uu-uv)
    g6(Nx+1,j) = w_eq6*T(Nx,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
    g7(Nx+1,j) = w_eq7*T(Nx,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
    T(Nx+1,j) = T(Nx,j)
 end do

 !======================================== NORTH


! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
!$OMP DO
  do i=0,Nx+1     
     g0(i,Ny+1) = g0(i ,Ny)
     g1(i,Ny+1) = g1(i ,Ny)
     g2(i,Ny+1) = g2(i ,Ny)
     g3(i,Ny+1) = g3(i ,Ny)
     g4(i,Ny+1) = g4(i ,Ny)
     g5(i,Ny+1) = g5(i ,Ny)
     g6(i,Ny+1) = g6(i ,Ny)
     g7(i,Ny+1) = g7(i ,Ny)
     g8(i,Ny+1) = g8(i ,Ny)
     T(i,Ny+1)   = T(i,Ny)
     phi(i,Ny+1) = phi(i,Ny)
  enddo

 !======================================== SOUTH

! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
!$OMP DO
  do i=0,Nx+1     
     g0(i,0) = g0(i ,1)
     g1(i,0) = g1(i ,1)
     g2(i,0) = g2(i ,1)
     g3(i,0) = g3(i ,1)
     g4(i,0) = g4(i ,1)
     g5(i,0) = g5(i ,1)
     g6(i,0) = g6(i ,1)
     g7(i,0) = g7(i ,1)
     g8(i,0) = g8(i ,1)
     T(i,0)   = T(i,1)
     phi(i,0) = phi(i,1)
  enddo
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_g"
#endif


 end subroutine BC_g
