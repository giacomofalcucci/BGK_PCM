 subroutine equilibrium

 use shared

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j)     !!!!,uu,vv,uv,upv,umv)
!$OMP SECTIONS
!$OMP SECTION
!!!$OMP DO
 do j = 1,Ny
    do i = 1,Nx

       if (flag(i,j) .eq. 0) then

          uu = u(i,j)*u(i,j)
          vv = v(i,j)*v(i,j)
          uv = 1.5d0*(uu+vv)
          upv = u(i,j)+v(i,j)
          umv = u(i,j)-v(i,j)
   
          ffeq(i,j,0) = w_eq(0)*rho(i,j)*(1.d0-uv)
          ffeq(i,j,1) = w_eq(1)*rho(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          ffeq(i,j,2) = w_eq(2)*rho(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          ffeq(i,j,3) = w_eq(3)*rho(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          ffeq(i,j,4) = w_eq(4)*rho(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          ffeq(i,j,5) = w_eq(5)*rho(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          ffeq(i,j,6) = w_eq(6)*rho(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          ffeq(i,j,7) = w_eq(7)*rho(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          ffeq(i,j,8) = w_eq(8)*rho(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
          
          ggeq(i,j,0) = w_eq(0)*T(i,j)*(1.d0-uv)
          ggeq(i,j,1) = w_eq(1)*T(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          ggeq(i,j,2) = w_eq(2)*T(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          ggeq(i,j,3) = w_eq(3)*T(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          ggeq(i,j,4) = w_eq(4)*T(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          ggeq(i,j,5) = w_eq(5)*T(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          ggeq(i,j,6) = w_eq(6)*T(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          ggeq(i,j,7) = w_eq(7)*T(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          ggeq(i,j,8) = w_eq(8)*T(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)

       end if

    end do
 end do
!$OMP END SECTIONS
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine equilibrium"
#endif

 end subroutine equilibrium
