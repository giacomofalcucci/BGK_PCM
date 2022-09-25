 subroutine setup

 use shared

 open(unit=112,file='input.dat',status='old')
    ! simulation parameters
    read(112,*)    id(1),  iter
    read(112,*)    id(2),  noutput
    read(112,*)    id(3),  dump
    ! number of nodes
    read(112,*)    id(4),  Nx
    read(112,*)    id(5),  Ny
    ! fluid flow parameters
    read(112,*)    id(6),  Re
    read(112,*)    id(7),  rho0
    read(112,*)    id(8),  rho1
    read(112,*)    id(9),  tau_f
    read(112,*)    id(10), tau_g
    read(112,*)    id(11), T0
    read(112,*)    id(12), T_in_sol
    read(112,*)    id(13), grav
    read(112,*)    id(14), eps
    read(112,*)    id(15), stefan_number
    read(112,*)    id(16), T_width
    ! characteristic length
    read(112,*)    id(17), length
    ! structured boundary conditions
    read(112,*)    id(18), westBC_st
    read(112,*)    id(19), eastBC_st
    read(112,*)    id(20), northBC_st
    read(112,*)    id(21), southBC_st
 close(112)

 allocate(u(0:Nx+1,0:Ny+1),v(0:Nx+1,0:Ny+1),rho(0:Nx+1,0:Ny+1))
 allocate(u2(0:Nx+1,0:Ny+1),v2(0:Nx+1,0:Ny+1),rho2(0:Nx+1,0:Ny+1))
 allocate(flag(0:Nx+1,0:Ny+1))
 allocate(x(0:Nx+1),y(0:Ny+1))
 allocate(feq(0:npop-1,0:Nx+1,0:Ny+1),f(0:npop-1,0:Nx+1,0:Ny+1))
 allocate(ftemp(0:npop-1,0:Nx+1,0:Ny+1))
 allocate(geq(0:npop-1,0:Nx+1,0:Ny+1),g(0:npop-1,0:Nx+1,0:Ny+1))
 allocate(T(0:Nx+1,0:Ny+1),T_prec(0:Nx+1,0:Ny+1))
 allocate(T2(0:Nx+1,0:Ny+1))
 allocate(phi(0:Nx+1,0:Ny+1),phi_prec(0:Nx+1,0:Ny+1))
 allocate(force_x(0:Nx+1,0:Ny+1),force_y(0:Nx+1,0:Ny+1))
 allocate(react_tot(1:Nx,1:Ny))

 w_eq(0) = 4.d0/9.d0
 w_eq(1) = 1.d0/9.d0
 w_eq(2) = 1.d0/9.d0
 w_eq(3) = 1.d0/9.d0
 w_eq(4) = 1.d0/9.d0
 w_eq(5) = 1.d0/36.d0
 w_eq(6) = 1.d0/36.d0
 w_eq(7) = 1.d0/36.d0
 w_eq(8) = 1.d0/36.d0

 x(0) = 0.d0
 y(0) = 0.d0
 do i = 1,Nx+1
    x(i) = x(i-1) + 1.d0
 end do
 do j = 1,Ny+1
    y(j) = y(j-1) + 1.d0
 end do

 cs = sqrt(1.d0/3.d0)
 cs2 = 1.d0/3.d0
 cs4 = 1.d0/9.d0

 flag = 0

 omega_f = 1.d0/tau_f
 omega_g = 1.d0/tau_g
 kvisc = cs2*(tau_f-0.5d0)
 alpha = cs2*(tau_g-0.5d0)

 samplex = 40
 sampley = 20
 
 a = Nx/samplex
 b = Ny/sampley

 allocate(react_plus(1:Nx,1:Ny),react_minus(1:Nx,1:Ny), &
          k_plus(1:Nx,1:Ny),k_minus(1:Nx,1:Ny)) 
 allocate(Nu(1:samplex),Tavg(1:samplex),gradTwall(1:samplex))
 allocate(Nu2(1:samplex),Tavg2(1:samplex),gradTwall2(1:samplex))
 allocate(Nu3(1:samplex))


 end subroutine setup
