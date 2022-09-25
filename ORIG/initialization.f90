 subroutine initialization

 use shared

 interf = int(Nx/40)

 !rho = rho0
! T = T0
 Tin = T0 !5.d0*T0

 u0 = Re*kvisc/length
 v0 = 0.d0
 Mach = u0/cs

do j=0,Ny+1
!  do i=0,int(Nx/2)-1
 do i =0, interf
    T(i,j)   = T0
    rho(i,j) = rho1 
  enddo
!  do i=int(Nx/2),Nx+1
  do i=interf+1,Nx+1
    T(i,j)   = T_in_sol
    rho(i,j) = rho0
  enddo
enddo


 write(*,*) '****************************************************'
 write(*,*) 'alpha= ' , alpha
 write(*,*) 'kvisc= ' , kvisc
 write(*,*) 'u0= ' , u0
 write(*,*) 'T0= ' , T0
 write(*,*) '****************************************************'

 if (u0 .ge. cs/3.d0) then
    write(*,*) '=============================='
    write(*,*) 'WARNING! Mach number too high!'
    write(*,*) '=============================='
 end if

 u = u0
 v = v0

 cx(0) = 0.0d0
 cx(1) = 1.0d0
 cx(2) = 0.0d0
 cx(3) = -1.0d0
 cx(4) = 0.0d0
 cx(5) = 1.0d0
 cx(6) = -1.0d0
 cx(7) = -1.0d0
 cx(8) = 1.0d0

 cy(0) = 0.0d0
 cy(1) = 0.0d0
 cy(2) = 1.0d0
 cy(3) = 0.0d0
 cy(4) = -1.0d0
 cy(5) = 1.0d0
 cy(6) = 1.0d0
 cy(7) = -1.0d0
 cy(8) = -1.0d0

 do j = 0,Ny+1
    T(0,j) = Tin
 end do

 do i = 0,Nx+1
    u(i,Ny+1) = 0.d0
    v(i,Ny+1) = 0.d0
    u(i,0) = 0.d0
    v(i,0) = 0.d0
 end do

 if (westBC_st .eq. 3) then
    u = 0.d0
    v = 0.d0
    do j = 0,Ny+1
       rho(0,j) = rho1
    end do
 end if

 do j = 1,Ny
    if     (westBC_st .eq. 1) then
       u(0,j) = u0
       v(0,j) = v0
    elseif (westBC_st .eq. 2) then
       u(0,j) = -(6.d0*u0/(Ny+1)**2.d0)*(y(j)**2.d0-(Ny+1)*y(j))
       v(0,j) = v0
    end if
 end do

 do j = 0,Ny+1
    do i = 0,Nx+1

       uu = u(i,j)*u(i,j)
       vv = v(i,j)*v(i,j)
       uv = 1.5d0*(uu+vv)
       upv = u(i,j)+v(i,j)
       umv = u(i,j)-v(i,j)

       f(0,i,j) = w_eq(0)*rho(i,j)*(1.d0-uv)
       f(1,i,j) = w_eq(1)*rho(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
       f(2,i,j) = w_eq(2)*rho(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
       f(3,i,j) = w_eq(3)*rho(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
       f(4,i,j) = w_eq(4)*rho(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
       f(5,i,j) = w_eq(5)*rho(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       f(6,i,j) = w_eq(6)*rho(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
       f(7,i,j) = w_eq(7)*rho(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       f(8,i,j) = w_eq(8)*rho(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)

       g(0,i,j) = w_eq(0)*T(i,j)*(1.d0-uv)
       g(1,i,j) = w_eq(1)*T(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
       g(2,i,j) = w_eq(2)*T(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
       g(3,i,j) = w_eq(3)*T(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
       g(4,i,j) = w_eq(4)*T(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
       g(5,i,j) = w_eq(5)*T(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
       g(6,i,j) = w_eq(6)*T(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
       g(7,i,j) = w_eq(7)*T(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
       g(8,i,j) = w_eq(8)*T(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)

   end do
 end do

 do j=0,Ny+1
  do i=0,Nx+1

! Phase Field initialization: "-1" -> SOLID, "+1" -> LIQUID
  phi(i,j)      = -1.0d0
  phi_prec(i,j) = -1.0d0 

  enddo
 enddo

 do j=0,Ny+1
!  do i=0,int(Nx/2)
  do i=0,interf

! Phase Field initialization: "-1" -> SOLID, "+1" -> LIQUID
  phi(i,j)      = 1.0d0
  phi_prec(i,j) = 1.0d0

  enddo
 enddo


 end subroutine initialization
