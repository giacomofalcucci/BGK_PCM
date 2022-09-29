 subroutine BC_f

 use shared

 !======================================== WEST

 if     (westBC_st .eq. 0) then
    ! periodic
    do j = 0,Ny+1
       do k = 0,npop-1
          f(k,0,j) = f(k,Nx,j)
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
       f(1,0,j) = w_eq(1)*rho(1,j)*(1.d0+3.d0*u0+4.5d0*uu-uv)
       f(5,0,j) = w_eq(5)*rho(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f(8,0,j) = w_eq(8)*rho(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
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
       f(1,0,j) = w_eq(1)*rho(1,j)*(1.d0+3.d0*u(0,j)+4.5d0*uu-uv)
       f(5,0,j) = w_eq(5)*rho(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f(8,0,j) = w_eq(8)*rho(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 elseif (westBC_st .eq. 3) then
    ! fixed constant pressure, velocity gradient = 0
    do j = 0,Ny+1
       uu = u(1,j)*u(1,j)
       vv = v(1,j)*v(1,j)
       uv = 1.5d0*(uu+vv)
       upv = u(1,j)+v(1,j)
       umv = u(1,j)-v(1,j)
       f(1,0,j) = w_eq(1)*rho1*(1.d0+3.d0*u(1,j)+4.5d0*uu-uv)
       f(5,0,j) = w_eq(5)*rho1*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f(8,0,j) = w_eq(8)*rho1*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u(0,j) = u(1,j)
       v(0,j) = v(1,j)
       rho(0,j) = rho1
    end do
 elseif (westBC_st .eq. 4) then
    ! bounceback
    do j = 1,Ny
       f(1,0,j  ) = f(3,1,j)
       f(8,0,j  ) = f(6,1,j-1)
       f(5,0,j  ) = f(7,1,j+1)
    end do
 end if

 !======================================== EAST

 if     (eastBC_st .eq. 0) then
    ! periodic
    do j = 0,Ny+1
       do k = 0,npop-1
          f(k,Nx+1,j) = f(k,1,j)
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
       f(3,Nx+1,j) = w_eq(3)*rho0*(1.d0-3.d0*u(Nx,j)+4.5d0*uu-uv)
       f(6,Nx+1,j) = w_eq(6)*rho0*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
       f(7,Nx+1,j) = w_eq(7)*rho0*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       u(Nx+1,j) = u(Nx,j)
       v(Nx+1,j) = v(Nx,j)
       rho(Nx+1,j) = rho0
    end do
 elseif (eastBC_st .eq. 4) then
    ! bounceback
    do j = 1,Ny
       f(3,Nx+1,j  ) = f(1,Nx,j)
       f(6,Nx+1,j  ) = f(8,Nx,j+1)
       f(7,Nx+1,j  ) = f(5,Nx,j-1)
    end do
 end if

 !======================================== NORTH

 if     (northBC_st .eq. 0) then
    ! periodic
    do i = 0,Nx+1
       do k = 0,npop-1
          f(k,i,Ny+1) = f(k,i,1)
       end do
    end do
 elseif (northBC_st .eq. 4) then
    ! bounceback
    do i = 1,Nx
       f(4,i,Ny+1) = f(2,i,Ny)
       f(7,i,Ny+1) = f(5,i-1,Ny)
       f(8,i,Ny+1) = f(6,i+1,Ny)
       rho(i,Ny+1) = rho(i,Ny)
       u(i,Ny+1) = 0.d0
       v(i,Ny+1) = 0.d0
    end do

       f(8,0,Ny+1) = f(6,1,Ny) 

elseif (northBC_st .eq. 5) then
    ! symmetric
    do i = 0,Nx+1
       f(4,i,Ny+1) = f(2,i,Ny-1)
       f(7,i,Ny+1) = f(6,i,Ny-1)
       f(8,i,Ny+1) = f(5,i,Ny-1)
    end do
 elseif (northBC_st .eq. 6) then
    ! open-boundary
    do i = 0,Nx+1
       f(4,i,Ny+1) = f(4,i,Ny)
       f(7,i,Ny+1) = f(7,i,Ny)
       f(8,i,Ny+1) = f(8,i,Ny)
    end do
 elseif (northBC_st .eq. 1) then
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
    do i = 0,Nx+1
       f(4,i,Ny+1) = w_eq(4)*rho(i,Ny)*(1.d0-3.d0*v0+4.5d0*vv-uv)
       f(7,i,Ny+1) = w_eq(7)*rho(i,Ny)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       f(8,i,Ny+1) = w_eq(8)*rho(i,Ny)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 end if

 !======================================== SOUTH

 if     (southBC_st .eq. 0) then
    ! periodic
    do i = 0,Nx+1
       do k = 0,npop-1
          f(k,i,0) = f(k,i,Ny)
       end do
    end do
 elseif (southBC_st .eq. 4) then
    ! bounceback
    do i = 1,Nx
       f(2,i,0) = f(4,i,1)
       f(5,i,0) = f(7,i+1,1)
       f(6,i,0) = f(8,i-1,1)
       rho(i,0) = rho(i,1)
       u(i,0) = 0.d0
       v(i,0) = 0.d0
    end do

      f(5,0,0)  = f(7,1,1)
 elseif (southBC_st .eq. 6) then
    ! open-boundary
    do i = 0,Nx+1
       f(2,i,0) = f(2,i,1)
       f(5,i,0) = f(5,i,1)
       f(6,i,0) = f(6,i,1)
    end do
 end if

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_f"
#endif


 end subroutine BC_f
