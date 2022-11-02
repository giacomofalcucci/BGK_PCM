 subroutine BC_f

 use shared

 !======================================== WEST

 if     (westBC_st .eq. 0) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! periodic
    do concurrent(j=0:Ny+1)
         f0(0,j) = f0(Nx,j)
         f1(0,j) = f1(Nx,j)
         f2(0,j) = f2(Nx,j)
         f3(0,j) = f3(Nx,j)
         f4(0,j) = f4(Nx,j)
         f5(0,j) = f5(Nx,j)
         f6(0,j) = f6(Nx,j)
         f7(0,j) = f7(Nx,j)
         f8(0,j) = f8(Nx,j)
    end do
!    
 elseif (westBC_st .eq. 1) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
    do concurrent(j=1:Ny)
       f1(0,j) = w_eq1*rho(1,j)*(1.d0+3.d0*u0+4.5d0*uu-uv)
       f5(0,j) = w_eq5*rho(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f8(0,j) = w_eq8*rho(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u(0,j) = u0
       v(0,j) = v0
       rho(0,j) = rho(1,j)
    end do
!    
 elseif (westBC_st .eq. 2) then 
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! fixed parabolic velocity, pressure gradient = 0
    do concurrent(j=0:Ny+1)
       uu = u(0,j)*u(0,j)
       vv = v0*v0
       uv = 1.5d0*(uu+vv)
       upv = u(0,j)+v0
       umv = u(0,j)-v0
       f1(0,j) = w_eq1*rho(1,j)*(1.d0+3.d0*u(0,j)+4.5d0*uu-uv)
       f5(0,j) = w_eq5*rho(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f8(0,j) = w_eq8*rho(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
!    
 elseif (westBC_st .eq. 3) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! fixed constant pressure, velocity gradient = 0
    do concurrent(j=0:Ny+1)
       uu = u(1,j)*u(1,j)
       vv = v(1,j)*v(1,j)
       uv = 1.5d0*(uu+vv)
       upv = u(1,j)+v(1,j)
       umv = u(1,j)-v(1,j)
       f1(0,j) = w_eq1*rho1*(1.d0+3.d0*u(1,j)+4.5d0*uu-uv)
       f5(0,j) = w_eq5*rho1*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f8(0,j) = w_eq8*rho1*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
       u(0,j) = u(1,j)
       v(0,j) = v(1,j)
       rho(0,j) = rho1
    end do
!    
 elseif (westBC_st .eq. 4) then
    ! bounceback
    do concurrent(j=1:Ny)
       f1(0,j  ) = f3(1,j)
       f8(0,j  ) = f6(1,j-1)
       f5(0,j  ) = f7(1,j+1)
       rho(0,j) = rho(1,j)
       u(0,j) = 0.d0
       v(0,j) = 0.d0
    end do
!
 end if

 !======================================== EAST

 if     (eastBC_st .eq. 0) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! periodic
    do concurrent(j=0:Ny+1)
       f0(Nx+1,j) = f0(1,j)
       f1(Nx+1,j) = f1(1,j)
       f2(Nx+1,j) = f2(1,j)
       f3(Nx+1,j) = f3(1,j)
       f4(Nx+1,j) = f4(1,j)
       f5(Nx+1,j) = f5(1,j)
       f6(Nx+1,j) = f6(1,j)
       f7(Nx+1,j) = f7(1,j)
       f8(Nx+1,j) = f8(1,j)
    end do
!    
 elseif (eastBC_st .eq. 3) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! fixed constant pressure, velocity gradient = 0
    do concurrent(j=1:Ny)
       uu = u(Nx,j)*u(Nx,j)
       vv = v(Nx,j)*v(Nx,j)
       uv = 1.5d0*(uu+vv)
       upv = u(Nx,j)+v(Nx,j)
       umv = u(Nx,j)-v(Nx,j)
       f3(Nx+1,j) = w_eq3*rho0*(1.d0-3.d0*u(Nx,j)+4.5d0*uu-uv)
       f6(Nx+1,j) = w_eq6*rho0*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
       f7(Nx+1,j) = w_eq7*rho0*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       u(Nx+1,j) = u(Nx,j)
       v(Nx+1,j) = v(Nx,j)
       rho(Nx+1,j) = rho0
    end do
!    
 elseif (eastBC_st .eq. 4) then
    ! bounceback
    do concurrent(j=1:Ny)
       f3(Nx+1,j  ) = f1(Nx,j)
       f6(Nx+1,j  ) = f8(Nx,j+1)
       f7(Nx+1,j  ) = f5(Nx,j-1)
       rho(Nx+1,j) = rho(Nx,j)
       u(Nx+1,j) = 0.d0
       v(Nx+1,j) = 0.d0
    end do
!
!   Corner fix
    f6(Nx+1,0   ) = f8(Nx,1)
    f7(Nx+1,Ny+1) = f5(Nx,Ny)
!
 end if

 !======================================== NORTH

 if     (northBC_st .eq. 0) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! periodic
    do concurrent(i=0:Nx+1)
       f0(i,Ny+1) = f0(i,1)
       f1(i,Ny+1) = f1(i,1)
       f2(i,Ny+1) = f2(i,1)
       f3(i,Ny+1) = f3(i,1)
       f4(i,Ny+1) = f4(i,1)
       f5(i,Ny+1) = f5(i,1)
       f6(i,Ny+1) = f6(i,1)
       f7(i,Ny+1) = f7(i,1)
       f8(i,Ny+1) = f8(i,1)
    end do
!         
 elseif (northBC_st .eq. 4) then
    ! bounceback
    do concurrent(i=1:Nx)
       f4(i,Ny+1) = f2(i,Ny)
       f7(i,Ny+1) = f5(i-1,Ny)
       f8(i,Ny+1) = f6(i+1,Ny)
       rho(i,Ny+1) = rho(i,Ny)
       u(i,Ny+1) = 0.d0
       v(i,Ny+1) = 0.d0
    end do

!   Corner fix
    f8(0,Ny+1) = f6(1,Ny) 
!
elseif (northBC_st .eq. 5) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! symmetric
    do concurrent(i=0:Nx+1)
       f4(i,Ny+1) = f2(i,Ny-1)
       f7(i,Ny+1) = f6(i,Ny-1)
       f8(i,Ny+1) = f5(i,Ny-1)
    end do
 elseif (northBC_st .eq. 6) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! open-boundary
    do concurrent(i=0:Nx+1)
       f4(i,Ny+1) = f4(i,Ny)
       f7(i,Ny+1) = f7(i,Ny)
       f8(i,Ny+1) = f8(i,Ny)
    end do
 elseif (northBC_st .eq. 1) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! fixed constant velocity, pressure gradient = 0
    uu = u0*u0
    vv = v0*v0
    uv = 1.5d0*(uu+vv)
    upv = u0+v0
    umv = u0-v0
    do concurrent(i=0:Nx+1)
       f4(i,Ny+1) = w_eq4*rho(i,Ny)*(1.d0-3.d0*v0+4.5d0*vv-uv)
       f7(i,Ny+1) = w_eq7*rho(i,Ny)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       f8(i,Ny+1) = w_eq8*rho(i,Ny)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    end do
 end if

 !======================================== SOUTH

 if     (southBC_st .eq. 0) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! periodic
    do concurrent(i=0:Nx+1)
       f0(i,0) = f0(i,Ny)
       f1(i,0) = f1(i,Ny)
       f2(i,0) = f2(i,Ny)
       f3(i,0) = f3(i,Ny)
       f4(i,0) = f4(i,Ny)
       f5(i,0) = f5(i,Ny)
       f6(i,0) = f6(i,Ny)
       f7(i,0) = f7(i,Ny)
       f8(i,0) = f8(i,Ny)
    end do
!         
 elseif (southBC_st .eq. 4) then
!         
    ! bounceback
    do concurrent(i=1:Nx)
       f2(i,0) = f4(i  ,1)
       f5(i,0) = f7(i+1,1)
       f6(i,0) = f8(i-1,1)
       rho(i,0) = rho(i,1)
       u(i,0) = 0.d0
       v(i,0) = 0.d0
    end do

!   Corner fix
    f5(0,0)  = f7(1,1)
!
 elseif (southBC_st .eq. 6) then
!         
    write(6,*) "ERROR: BC still to  validate"
    stop
!    
    ! open-boundary
    do concurrent(i=0:Nx+1)
       f2(i,0) = f2(i,1)
       f5(i,0) = f5(i,1)
       f6(i,0) = f6(i,1)
    end do
 end if

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_f"
#endif


 end subroutine BC_f
