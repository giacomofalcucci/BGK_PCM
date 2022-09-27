 subroutine BC_f

 use shared

 !======================================== WEST

 if     (westBC_st .eq. 0) then
    ! periodic
    do k = 0,npop-1
       do j = 0,Ny+1
          ff(0,j,k) = ff(Nx,j,k)
       end do
    end do
 elseif (westBC_st .eq. 1) then
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
    do j = 1,Ny
       ff(0,j,1) = w_eq(1)*rho(1,j)*(1.d0+3.d0*u0+4.5d0*uu-uv)
       ff(0,j,5) = w_eq(5)*rho(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       ff(0,j,8) = w_eq(8)*rho(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u(0,j) = u0
       v(0,j) = v0
       rho(0,j) = rho(1,j)
    end do
 elseif (westBC_st .eq. 2) then 
    ! fixed parabolic velocity, pressure gradient = 0
    do j = 0,Ny+1
       uu = u(0,j)*u(0,j)
       vv = v0*v0
       uv = 1.5d0*(uu+vv)
       upv = u(0,j)+v0
       umv = u(0,j)-v0
       ff(0,j,1) = w_eq(1)*rho(1,j)*(1.d0+3.d0*u(0,j)+4.5d0*uu-uv)
       ff(0,j,5) = w_eq(5)*rho(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       ff(0,j,8) = w_eq(8)*rho(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 elseif (westBC_st .eq. 3) then
    ! fixed constant pressure, velocity gradient = 0
    do j = 0,Ny+1
       uu = u(1,j)*u(1,j)
       vv = v(1,j)*v(1,j)
       uv = 1.5d0*(uu+vv)
       upv = u(1,j)+v(1,j)
       umv = u(1,j)-v(1,j)
       ff(0,j,1) = w_eq(1)*rho1*(1.d0+3.d0*u(1,j)+4.5d0*uu-uv)
       ff(0,j,5) = w_eq(5)*rho1*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       ff(0,j,8) = w_eq(8)*rho1*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u(0,j) = u(1,j)
       v(0,j) = v(1,j)
       rho(0,j) = rho1
    end do
 elseif (westBC_st .eq. 4) then
    ! bounceback
    do j = 1,Ny
       ff(0,j  ,1) = ff(1,j  ,3)
       ff(0,j  ,8) = ff(1,j-1,6)
       ff(0,j  ,5) = ff(1,j+1,7)
    end do
 end if

 !======================================== EAST

 if     (eastBC_st .eq. 0) then
    ! periodic
    do k = 0,npop-1
       do j = 0,Ny+1
          ff(Nx+1,j,k) = ff(1,j,k)
       end do
    end do
 elseif (eastBC_st .eq. 3) then
    ! fixed constant pressure, velocity gradient = 0
    do j = 1,Ny
       uu = u(Nx,j)*u(Nx,j)
       vv = v(Nx,j)*v(Nx,j)
       uv = 1.5d0*(uu+vv)
       upv = u(Nx,j)+v(Nx,j)
       umv = u(Nx,j)-v(Nx,j)
       ff(Nx+1,j,3) = w_eq(3)*rho0*(1.d0-3.d0*u(Nx,j)+4.5d0*uu-uv)
       ff(Nx+1,j,6) = w_eq(6)*rho0*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
       ff(Nx+1,j,7) = w_eq(7)*rho0*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       u(Nx+1,j) = u(Nx,j)
       v(Nx+1,j) = v(Nx,j)
       rho(Nx+1,j) = rho0
    end do
 elseif (eastBC_st .eq. 4) then
    ! bounceback
    do j = 1,Ny
       ff(Nx+1,j,3) = ff(Nx,j  ,1)
       ff(Nx+1,j,6) = ff(Nx,j+1,8)
       ff(Nx+1,j,7) = ff(Nx,j-1,5)
    end do
 end if

 !======================================== NORTH

 if     (northBC_st .eq. 0) then
    ! periodic
    do k = 0,npop-1
       do i = 0,Nx+1
          ff(i,Ny+1,k) = ff(i,1,k)
       end do
    end do
 elseif (northBC_st .eq. 4) then
    ! bounceback
    do i = 1,Nx
       ff(i,Ny+1,4) = ff(i  ,Ny,2)
       ff(i,Ny+1,7) = ff(i-1,Ny,5)
       ff(i,Ny+1,8) = ff(i+1,Ny,6)
       rho(i,Ny+1) = rho(i,Ny)
       u(i,Ny+1) = 0.d0
       v(i,Ny+1) = 0.d0
    end do

       ff(0,Ny+1,8) = ff(1,Ny,6) 

elseif (northBC_st .eq. 5) then
    ! symmetric
    do i = 0,Nx+1
       ff(i,Ny+1,4) = ff(i,Ny-1,2)
       ff(i,Ny+1,7) = ff(i,Ny-1,6)
       ff(i,Ny+1,8) = ff(i,Ny-1,5)
    end do
 elseif (northBC_st .eq. 6) then
    ! open-boundary
    do i = 0,Nx+1
       ff(i,Ny+1,4) = ff(i,Ny,4)
       ff(i,Ny+1,7) = ff(i,Ny,7)
       ff(i,Ny+1,8) = ff(i,Ny,8)
    end do
 elseif (northBC_st .eq. 1) then
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
    do i = 0,Nx+1
       ff(i,Ny+1,4) = w_eq(4)*rho(i,Ny)*(1.d0-3.d0*v0+4.5d0*vv-uv)
       ff(i,Ny+1,7) = w_eq(7)*rho(i,Ny)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       ff(i,Ny+1,8) = w_eq(8)*rho(i,Ny)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 end if

 !======================================== SOUTH

 if     (southBC_st .eq. 0) then
    ! periodic
    do k = 0,npop-1
       do i = 0,Nx+1
          ff(i,0,k) = ff(i,Ny,k)
       end do
    end do
 elseif (southBC_st .eq. 4) then
    ! bounceback
    do i = 1,Nx
       ff(i,0,2) = ff(i  ,1,4)
       ff(i,0,5) = ff(i+1,1,7)
       ff(i,0,6) = ff(i-1,1,8)
       rho(i,0) = rho(i,1)
       u(i,0) = 0.d0
       v(i,0) = 0.d0
    end do

      ff(0,0,5)  = ff(1,1,7)
 elseif (southBC_st .eq. 6) then
    ! open-boundary
    do i = 0,Nx+1
       ff(i,0,2) = ff(i,1,2)
       ff(i,0,5) = ff(i,1,5)
       ff(i,0,6) = ff(i,1,6)
    end do
 end if

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_f"
#endif


 end subroutine BC_f
