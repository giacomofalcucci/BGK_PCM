 subroutine BC_f

 use shared

 !======================================== WEST

 if     (westBC_st .eq. 0) then
    ! periodic
      !$cuf kernel  do <<<*,*>>>
    do j = 0,Ny+1
         f0_gpu(0,j) = f0_gpu(Nx,j)
         f1_gpu(0,j) = f1_gpu(Nx,j)
         f2_gpu(0,j) = f2_gpu(Nx,j)
         f3_gpu(0,j) = f3_gpu(Nx,j)
         f4_gpu(0,j) = f4_gpu(Nx,j)
         f5_gpu(0,j) = f5_gpu(Nx,j)
         f6_gpu(0,j) = f6_gpu(Nx,j)
         f7_gpu(0,j) = f7_gpu(Nx,j)
         f8_gpu(0,j) = f8_gpu(Nx,j)
    end do
 elseif (westBC_st .eq. 1) then
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
    !$cuf kernel  do <<<*,*>>>
    do j = 1,Ny
       f1_gpu(0,j) = w_eq1*rho_gpu(1,j)*(1.d0+3.d0*u0+4.5d0*uu-uv)
       f5_gpu(0,j) = w_eq5*rho_gpu(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f8_gpu(0,j) = w_eq8*rho_gpu(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u_gpu(0,j) = u0
       v_gpu(0,j) = v0
       rho_gpu(0,j) = rho_gpu(1,j)
    end do
 elseif (westBC_st .eq. 2) then 
    ! fixed parabolic velocity, pressure gradient = 0
      !$cuf kernel  do <<<*,*>>>
    do j = 0,Ny+1
       uu = u_gpu(0,j)*u_gpu(0,j)
       vv = v0*v0
       uv = 1.5d0*(uu+vv)
       upv = u_gpu(0,j)+v0
       umv = u_gpu(0,j)-v0
       f1_gpu(0,j) = w_eq1*rho_gpu(1,j)*(1.d0+3.d0*u_gpu(0,j)+4.5d0*uu-uv)
       f5_gpu(0,j) = w_eq5*rho_gpu(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f8_gpu(0,j) = w_eq8*rho_gpu(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 elseif (westBC_st .eq. 3) then
    ! fixed constant pressure, velocity gradient = 0
      !$cuf kernel  do <<<*,*>>>
    do j = 0,Ny+1
       uu = u_gpu(1,j)*u_gpu(1,j)
       vv = v_gpu(1,j)*v_gpu(1,j)
       uv = 1.5d0*(uu+vv)
       upv = u_gpu(1,j)+v_gpu(1,j)
       umv = u_gpu(1,j)-v_gpu(1,j)
       f1_gpu(0,j) = w_eq1*rho1*(1.d0+3.d0*u_gpu(1,j)+4.5d0*uu-uv)
       f5_gpu(0,j) = w_eq5*rho1*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f8_gpu(0,j) = w_eq8*rho1*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u_gpu(0,j) = u_gpu(1,j)
       v_gpu(0,j) = v_gpu(1,j)
       rho_gpu(0,j) = rho1
    end do
 elseif (westBC_st .eq. 4) then
    ! bounceback
      !$cuf kernel  do <<<*,*>>>
    do j = 1,Ny
       f1_gpu(0,j  ) = f3_gpu(1,j)
       f8_gpu(0,j  ) = f6_gpu(1,j-1)
       f5_gpu(0,j  ) = f7_gpu(1,j+1)
    end do
 end if

 !======================================== EAST

 if     (eastBC_st .eq. 0) then
    ! periodic
      !$cuf kernel  do <<<*,*>>>
    do j = 0,Ny+1
       f0_gpu(Nx+1,j) = f0_gpu(1,j)
       f1_gpu(Nx+1,j) = f1_gpu(1,j)
       f2_gpu(Nx+1,j) = f2_gpu(1,j)
       f3_gpu(Nx+1,j) = f3_gpu(1,j)
       f4_gpu(Nx+1,j) = f4_gpu(1,j)
       f5_gpu(Nx+1,j) = f5_gpu(1,j)
       f6_gpu(Nx+1,j) = f6_gpu(1,j)
       f7_gpu(Nx+1,j) = f7_gpu(1,j)
       f8_gpu(Nx+1,j) = f8_gpu(1,j)
    end do
 elseif (eastBC_st .eq. 3) then
    ! fixed constant pressure, velocity gradient = 0
      !$cuf kernel  do <<<*,*>>>
    do j = 1,Ny
       uu = u_gpu(Nx,j)*u_gpu(Nx,j)
       vv = v_gpu(Nx,j)*v_gpu(Nx,j)
       uv = 1.5d0*(uu+vv)
       upv = u_gpu(Nx,j)+v_gpu(Nx,j)
       umv = u_gpu(Nx,j)-v_gpu(Nx,j)
       f3_gpu(Nx+1,j) = w_eq3*rho0*(1.d0-3.d0*u_gpu(Nx,j)+4.5d0*uu-uv)
       f6_gpu(Nx+1,j) = w_eq6*rho0*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
       f7_gpu(Nx+1,j) = w_eq7*rho0*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       u_gpu(Nx+1,j) = u_gpu(Nx,j)
       v_gpu(Nx+1,j) = v_gpu(Nx,j)
       rho_gpu(Nx+1,j) = rho0
    end do
 elseif (eastBC_st .eq. 4) then
    ! bounceback
      !$cuf kernel  do <<<*,*>>>
    do j = 1,Ny
       f3_gpu(Nx+1,j  ) = f1_gpu(Nx,j)
       f6_gpu(Nx+1,j  ) = f8_gpu(Nx,j+1)
       f7_gpu(Nx+1,j  ) = f5_gpu(Nx,j-1)
    end do
 end if

 !======================================== NORTH

 if     (northBC_st .eq. 0) then
    ! periodic
      !$cuf kernel  do <<<*,*>>>
    do i = 0,Nx+1
       f0_gpu(i,Ny+1) = f0_gpu(i,1)
       f1_gpu(i,Ny+1) = f1_gpu(i,1)
       f2_gpu(i,Ny+1) = f2_gpu(i,1)
       f3_gpu(i,Ny+1) = f3_gpu(i,1)
       f4_gpu(i,Ny+1) = f4_gpu(i,1)
       f5_gpu(i,Ny+1) = f5_gpu(i,1)
       f6_gpu(i,Ny+1) = f6_gpu(i,1)
       f7_gpu(i,Ny+1) = f7_gpu(i,1)
       f8_gpu(i,Ny+1) = f8_gpu(i,1)
    end do
 elseif (northBC_st .eq. 4) then
    ! bounceback
      !$cuf kernel  do <<<*,*>>>
    do i = 1,Nx
       f4_gpu(i,Ny+1) = f2_gpu(i,Ny)
       f7_gpu(i,Ny+1) = f5_gpu(i-1,Ny)
       f8_gpu(i,Ny+1) = f6_gpu(i+1,Ny)
       rho_gpu(i,Ny+1) = rho_gpu(i,Ny)
       u_gpu(i,Ny+1) = 0.d0
       v_gpu(i,Ny+1) = 0.d0
    end do

! To understand to keep or not
!GAAA       f8_gpu(0,Ny+1) = f6_gpu(1,Ny) 

elseif (northBC_st .eq. 5) then
    ! symmetric
      !$cuf kernel  do <<<*,*>>>
    do i = 0,Nx+1
       f4_gpu(i,Ny+1) = f2_gpu(i,Ny-1)
       f7_gpu(i,Ny+1) = f6_gpu(i,Ny-1)
       f8_gpu(i,Ny+1) = f5_gpu(i,Ny-1)
    end do
 elseif (northBC_st .eq. 6) then
    ! open-boundary
      !$cuf kernel  do <<<*,*>>>
    do i = 0,Nx+1
       f4_gpu(i,Ny+1) = f4_gpu(i,Ny)
       f7_gpu(i,Ny+1) = f7_gpu(i,Ny)
       f8_gpu(i,Ny+1) = f8_gpu(i,Ny)
    end do
 elseif (northBC_st .eq. 1) then
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
      !$cuf kernel  do <<<*,*>>>
    do i = 0,Nx+1
       f4_gpu(i,Ny+1) = w_eq4*rho_gpu(i,Ny)*(1.d0-3.d0*v0+4.5d0*vv-uv)
       f7_gpu(i,Ny+1) = w_eq7*rho_gpu(i,Ny)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       f8_gpu(i,Ny+1) = w_eq8*rho_gpu(i,Ny)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 end if

 !======================================== SOUTH

 if     (southBC_st .eq. 0) then
    ! periodic
      !$cuf kernel  do <<<*,*>>>
    do i = 0,Nx+1
       f0_gpu(i,0) = f0_gpu(i,Ny)
       f1_gpu(i,0) = f1_gpu(i,Ny)
       f2_gpu(i,0) = f2_gpu(i,Ny)
       f3_gpu(i,0) = f3_gpu(i,Ny)
       f4_gpu(i,0) = f4_gpu(i,Ny)
       f5_gpu(i,0) = f5_gpu(i,Ny)
       f6_gpu(i,0) = f6_gpu(i,Ny)
       f7_gpu(i,0) = f7_gpu(i,Ny)
       f8_gpu(i,0) = f8_gpu(i,Ny)
    end do
 elseif (southBC_st .eq. 4) then
    ! bounceback
      !$cuf kernel  do <<<*,*>>>
    do i = 1,Nx
       f2_gpu(i,0) = f4_gpu(i  ,1)
       f5_gpu(i,0) = f7_gpu(i+1,1)
       f6_gpu(i,0) = f8_gpu(i-1,1)
       rho_gpu(i,0) = rho_gpu(i,1)
       u_gpu(i,0) = 0.d0
       v_gpu(i,0) = 0.d0
    end do

! To understand to keep or not
!GAAA      f5_gpu(0,0)  = f7_gpu(1,1)

 elseif (southBC_st .eq. 6) then
    ! open-boundary
      !$cuf kernel  do <<<*,*>>>
    do i = 0,Nx+1
       f2_gpu(i,0) = f2_gpu(i,1)
       f5_gpu(i,0) = f5_gpu(i,1)
       f6_gpu(i,0) = f6_gpu(i,1)
    end do
 end if

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_f"
#endif


 end subroutine BC_f
