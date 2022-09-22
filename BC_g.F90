!-----------------------------------------------------------
        subroutine BC_g
!-----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

!!! AT THE MOMENT   ->  ONLY v_INLET e p_OUTLET

!$acc kernels
!$acc loop independent
        do j = 1, ny
! WEST
          uu = u1(0,j)*u1(0,j)
          vv = v1(0,j)*v1(0,j)
          uv = 1.5d0*(uu+vv)
          upv = u1(0,j)+v1(0,j)
          umv = u1(0,j)-v1(0,j)
          g1(0,j) = cte09*Tin*(1.d0+3.d0*u1(0,j)+4.5d0*uu-uv)
          g5(0,j) = cte36*Tin*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          g8(0,j) = cte36*Tin*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
          T(0,j) = Tin
!
!  EAST
         uu = u1(nx+1,j)*u1(nx+1,j)
         vv = v1(nx+1,j)*v1(nx+1,j)
         uv = 1.5d0*(uu+vv)
         upv = u1(nx+1,j)+v1(nx+1,j)
         umv = u1(nx+1,j)-v1(nx+1,j)
         g3(nx+1,j) = cte09*T(nx,j)*(1.d0-3.d0*u1(nx+1,j)+4.5d0*uu-uv)
         g6(nx+1,j) = cte36*T(nx,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
         g7(nx+1,j) = cte36*T(nx,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
         T(nx+1,j) = T(nx,j)
        enddo
!$acc end kernels


!!!  AT THE MOMENT:   ONLY  WALL NO SLIP

!$acc kernels
!$acc loop independent
        do i = 1, nx
!   NORTH
           g4(i,ny+1)  = g2(i,ny)
           g7(i,ny+1)  = g5(i-1,ny)
           g8(i,ny+1)  = g6(i+1,ny)
           T(i,ny+1)   = T(i,ny)
           phi(i,ny+1) = phi(i,ny)
!   SOUTH
           g2(i,0) = g4(i,1)
           g5(i,0) = g7(i+1,1)
           g6(i,0) = g8(i-1,1)
           T(i,0)  = T(i,1)
           phi(i,0)= phi(i,1)
        enddo
!$acc end kernels

        T(0,0)         = T0
        T(nx+1,0)      = T_in_sol
        T(nx+1,ny+1)   = T_in_sol
        T(0,ny+1)      = T0

        g5(0,0)        = g7(1,1)
        g8(0,ny+1)     = g6(1,ny)
        g6(nx+1,0)     = g8(nx,1)
        g7(nx+1,ny+1)  = g5(nx,ny)

        end
