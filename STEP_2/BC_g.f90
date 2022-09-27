 subroutine BC_g

 use shared

 !======================================== WEST

 do j = 1,Ny
    uu = u(0,j)*u(0,j)
    vv = v(0,j)*v(0,j)
    uv = 1.5d0*(uu+vv)
    upv = u(0,j)+v(0,j)
    umv = u(0,j)-v(0,j)
    gg(0,j,1) = w_eq(1)*Tin*(1.d0+3.d0*u(0,j)+4.5d0*uu-uv)
    gg(0,j,5) = w_eq(5)*Tin*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
    gg(0,j,8) = w_eq(8)*Tin*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
    T(0,j) = Tin


 end do

 !======================================== EAST

 do j = 1,Ny
    uu = u(Nx+1,j)*u(Nx+1,j)
    vv = v(Nx+1,j)*v(Nx+1,j)
    uv = 1.5d0*(uu+vv)
    upv = u(Nx+1,j)+v(Nx+1,j)
    umv = u(Nx+1,j)-v(Nx+1,j)
    gg(Nx+1,j,3) = w_eq(3)*T(Nx,j)*(1.d0-3.d0*u(Nx+1,j)+4.5d0*uu-uv)
    gg(Nx+1,j,6) = w_eq(6)*T(Nx,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
    gg(Nx+1,j,7) = w_eq(7)*T(Nx,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
    T(Nx+1,j) = T(Nx,j)

 end do

 !======================================== NORTH


  do k=0,npop-1
     do i=0,Nx+1     ! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
        gg(i,Ny+1,k) = gg(i ,Ny,k)
     enddo
  enddo
  do i=0,Nx+1  
     T(i,Ny+1)   = T(i,Ny)
     phi(i,Ny+1) = phi(i,Ny)
  enddo

 !======================================== SOUTH


  do k=0,npop-1
     do i=0,Nx+1     ! Zero-gradient for the Temperature -> no thermal flux (<->adiabatic?!?!?)
        gg(i,0,k) = gg(i ,1,k)
     enddo
  enddo
  do i=0,Nx+1  
     T(i,0)   = T(i,1)
     phi(i,0) = phi(i,1)
  enddo

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine BC_g"
#endif


 end subroutine BC_g
