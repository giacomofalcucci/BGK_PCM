 subroutine collision
 
 use shared

 real(kind=mykind)   ::  feq0(0:npop-1)
 real(kind=mykind)   ::  geq0(0:npop-1)

 do j = 1,Ny
    do i = 1,Nx

          uu = u(i,j)*u(i,j)
          vv = v(i,j)*v(i,j)
          uv = 1.5d0*(uu+vv)
          upv = u(i,j)+v(i,j)
          umv = u(i,j)-v(i,j)
   
          feq0(0) = w_eq(0)*rho(i,j)*(1.d0-uv)
          feq0(1) = w_eq(1)*rho(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          feq0(2) = w_eq(2)*rho(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          feq0(3) = w_eq(3)*rho(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          feq0(4) = w_eq(4)*rho(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          feq0(5) = w_eq(5)*rho(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          feq0(6) = w_eq(6)*rho(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          feq0(7) = w_eq(7)*rho(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          feq0(8) = w_eq(8)*rho(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
          
          geq0(0) = w_eq(0)*T(i,j)*(1.d0-uv)
          geq0(1) = w_eq(1)*T(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          geq0(2) = w_eq(2)*T(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          geq0(3) = w_eq(3)*T(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          geq0(4) = w_eq(4)*T(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          geq0(5) = w_eq(5)*T(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          geq0(6) = w_eq(6)*T(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          geq0(7) = w_eq(7)*T(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          geq0(8) = w_eq(8)*T(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)

          f(0,i,j) = f(0,i,j)-omega_f*(f(0,i,j)-feq0(0))
          f(1,i,j) = f(1,i,j)-omega_f*(f(1,i,j)-feq0(1))
          f(2,i,j) = f(2,i,j)-omega_f*(f(2,i,j)-feq0(2))
          f(3,i,j) = f(3,i,j)-omega_f*(f(3,i,j)-feq0(3))
          f(4,i,j) = f(4,i,j)-omega_f*(f(4,i,j)-feq0(4))
          f(5,i,j) = f(5,i,j)-omega_f*(f(5,i,j)-feq0(5))
          f(6,i,j) = f(6,i,j)-omega_f*(f(6,i,j)-feq0(6))
          f(7,i,j) = f(7,i,j)-omega_f*(f(7,i,j)-feq0(7))
          f(8,i,j) = f(8,i,j)-omega_f*(f(8,i,j)-feq0(8))

          g(0,i,j) = g(0,i,j)-omega_g*(g(0,i,j)-geq0(0))
          g(1,i,j) = g(1,i,j)-omega_g*(g(1,i,j)-geq0(1))
          g(2,i,j) = g(2,i,j)-omega_g*(g(2,i,j)-geq0(2))
          g(3,i,j) = g(3,i,j)-omega_g*(g(3,i,j)-geq0(3))
          g(4,i,j) = g(4,i,j)-omega_g*(g(4,i,j)-geq0(4))
          g(5,i,j) = g(5,i,j)-omega_g*(g(5,i,j)-geq0(5))
          g(6,i,j) = g(6,i,j)-omega_g*(g(6,i,j)-geq0(6))
          g(7,i,j) = g(7,i,j)-omega_g*(g(7,i,j)-geq0(7))
          g(8,i,j) = g(8,i,j)-omega_g*(g(8,i,j)-geq0(8))
    end do
 end do
!$OMP END DO
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine collision"
#endif

 end subroutine collision
