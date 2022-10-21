 subroutine BC_g

 use shared

 !======================================== WEST

 !$cuf kernel  do <<<*,*>>>
 do j = 1,Ny
    uu = u_gpu(0,j)*u_gpu(0,j)
    vv = v_gpu(0,j)*v_gpu(0,j)
    uv = 1.5d0*(uu+vv)
    upv = u_gpu(0,j)+v_gpu(0,j)
    umv = u_gpu(0,j)-v_gpu(0,j)
    g1_gpu(0,j) = w_eq1*Tin*(1.d0+3.d0*u_gpu(0,j)+4.5d0*uu-uv)
    g5_gpu(0,j) = w_eq5*Tin*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
    g8_gpu(0,j) = w_eq8*Tin*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    T_gpu(0,j) = Tin
 end do

 !======================================== EAST

 !$cuf kernel  do <<<*,*>>>
 do j = 1,Ny
    uu = u_gpu(Nx+1,j)*u_gpu(Nx+1,j)
    vv = v_gpu(Nx+1,j)*v_gpu(Nx+1,j)
    uv = 1.5d0*(uu+vv)
    upv = u_gpu(Nx+1,j)+v_gpu(Nx+1,j)
    umv = u_gpu(Nx+1,j)-v_gpu(Nx+1,j)
    g3_gpu(Nx+1,j) = w_eq3*T_gpu(Nx,j)*(1.d0-3.d0*u_gpu(Nx+1,j)+4.5d0*uu-uv)
    g6_gpu(Nx+1,j) = w_eq6*T_gpu(Nx,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
    g7_gpu(Nx+1,j) = w_eq7*T_gpu(Nx,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
    T_gpu(Nx+1,j) = T_gpu(Nx,j)
 end do

 !======================================== NORTH


 !$cuf kernel  do <<<*,*>>>
  do i=0,Nx+1     ! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
     g0_gpu(i,Ny+1) = g0_gpu(i ,Ny)
     g1_gpu(i,Ny+1) = g1_gpu(i ,Ny)
     g2_gpu(i,Ny+1) = g2_gpu(i ,Ny)
     g3_gpu(i,Ny+1) = g3_gpu(i ,Ny)
     g4_gpu(i,Ny+1) = g4_gpu(i ,Ny)
     g5_gpu(i,Ny+1) = g5_gpu(i ,Ny)
     g6_gpu(i,Ny+1) = g6_gpu(i ,Ny)
     g7_gpu(i,Ny+1) = g7_gpu(i ,Ny)
     g8_gpu(i,Ny+1) = g8_gpu(i ,Ny)
     T_gpu(i,Ny+1)   = T_gpu(i,Ny)
     phi_gpu(i,Ny+1) = phi_gpu(i,Ny)
  enddo

 !======================================== SOUTH

 !$cuf kernel  do <<<*,*>>>
  do i=0,Nx+1     ! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
     g0_gpu(i,0) = g0_gpu(i ,1)
     g1_gpu(i,0) = g1_gpu(i ,1)
     g2_gpu(i,0) = g2_gpu(i ,1)
     g3_gpu(i,0) = g3_gpu(i ,1)
     g4_gpu(i,0) = g4_gpu(i ,1)
     g5_gpu(i,0) = g5_gpu(i ,1)
     g6_gpu(i,0) = g6_gpu(i ,1)
     g7_gpu(i,0) = g7_gpu(i ,1)
     g8_gpu(i,0) = g8_gpu(i ,1)
     T_gpu(i,0)   = T_gpu(i,1)
     phi_gpu(i,0) = phi_gpu(i,1)
  enddo

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_g"
#endif


 end subroutine BC_g
