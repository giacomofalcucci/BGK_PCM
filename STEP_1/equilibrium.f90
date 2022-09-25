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
   
          feq(0,i,j) = w_eq(0)*rho(i,j)*(1.d0-uv)
          feq(1,i,j) = w_eq(1)*rho(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          feq(2,i,j) = w_eq(2)*rho(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          feq(3,i,j) = w_eq(3)*rho(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          feq(4,i,j) = w_eq(4)*rho(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          feq(5,i,j) = w_eq(5)*rho(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          feq(6,i,j) = w_eq(6)*rho(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          feq(7,i,j) = w_eq(7)*rho(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          feq(8,i,j) = w_eq(8)*rho(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
          
          geq(0,i,j) = w_eq(0)*T(i,j)*(1.d0-uv)
          geq(1,i,j) = w_eq(1)*T(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          geq(2,i,j) = w_eq(2)*T(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          geq(3,i,j) = w_eq(3)*T(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          geq(4,i,j) = w_eq(4)*T(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          geq(5,i,j) = w_eq(5)*T(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          geq(6,i,j) = w_eq(6)*T(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          geq(7,i,j) = w_eq(7)*T(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          geq(8,i,j) = w_eq(8)*T(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)

       end if

    end do
 end do
!$OMP END SECTIONS
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine equilibrium"
#endif

 end subroutine equilibrium
