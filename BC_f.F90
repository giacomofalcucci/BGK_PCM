!-----------------------------------------------------------
        subroutine BC_f
!-----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        uu = u0*u0
        vv = v0*v0
        uv = 1.5d0*(uu+vv)
        upv = u0+v0
        umv = u0-v0

!!! AT THE MOMENT   ->  ONLY v_INLET e p_OUTLET

!$acc kernels
!$acc loop independent
        do j = 1, ny
!  WEST
           f1(0,j) = cte09*rhod1(1,j)*(1.d0+3.d0*u0+4.5d0*uu-uv)
           f5(0,j) = cte36*rhod1(1,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
           f8(0,j) = cte36*rhod1(1,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
           u1(0,j) = u0
           v1(0,j) = v0
           rhod1(0,j) = rhod1(1,j)
!  EAST
           uu = u1(nx,j)*u1(nx,j)
           vv = v1(nx,j)*v1(nx,j)
           uv = 1.5d0*(uu+vv)
           upv = u1(nx,j)+v1(nx,j)
           umv = u1(nx,j)-v1(nx,j)
           f3(nx+1,j) = cte09*rho0*(1.d0-3.d0*u1(Nx,j)+4.5d0*uu-uv)
           f6(nx+1,j) = cte36*rho0*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
           f7(nx+1,j) = cte36*rho0*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
           u1(nx+1,j) = u1(nx,j)
           v1(nx+1,j) = v1(nx,j)
           rhod1(nx+1,j) = rho0
        enddo
!$acc end kernels


!!!  AT THE MOMENT:   ONLY  WALL NO SLIP

!$acc kernels
!$acc loop independent
        do i = 1, nx
!   NORTH
           f4(i,ny+1)  = f2(i,ny)
           f7(i,ny+1)  = f5(i-1,ny)
           f8(i,ny+1)  = f6(i+1,ny)
           rhod1(i,ny+1) = rhod1(i,ny)
           u1(i,ny+1)   = 0.d0
           v1(i,ny+1)   = 0.d0
!   SOUTH
           f2(i,0) = f4(i,1)
           f5(i,0) = f7(i+1,1)
           f6(i,0) = f8(i-1,1)
           rhod1(i,0) = rhod1(i,1)
           u1(i,0) = 0.d0
           v1(i,0) = 0.d0
        enddo
!$acc end kernels

        rhod1(0,0)       = rhod1(1,0)
        rhod1(nx+1,0)    = rho0
        rhod1(nx+1,ny+1) = rho0
        rhod1(0,ny+1)    = rhod1(1,ny+1)

        f5(0,0)          = f7(1,1)
        f8(0,ny+1)       = f6(1,ny)
        f6(nx+1,0)       = f8(nx,1)
        f7(nx+1,ny+1)    = f5(nx,ny)

        end
