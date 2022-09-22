!--------------------------------------------------
        subroutine alloca

        use storage
        implicit none
!
! -------vector allocation
!
        allocate( fp0(0:nx+1,0:ny+1))
        allocate( fp1(0:nx+1,0:ny+1))
        allocate( fp2(0:nx+1,0:ny+1))
        allocate( fp3(0:nx+1,0:ny+1))
        allocate( fp4(0:nx+1,0:ny+1))
        allocate( fp5(0:nx+1,0:ny+1))
        allocate( fp6(0:nx+1,0:ny+1))
        allocate( fp7(0:nx+1,0:ny+1))
        allocate( fp8(0:nx+1,0:ny+1))
!        
        allocate(  f0(0:nx+1,0:ny+1))
        allocate(  f1(0:nx+1,0:ny+1))
        allocate(  f2(0:nx+1,0:ny+1))
        allocate(  f3(0:nx+1,0:ny+1))
        allocate(  f4(0:nx+1,0:ny+1))
        allocate(  f5(0:nx+1,0:ny+1))
        allocate(  f6(0:nx+1,0:ny+1))
        allocate(  f7(0:nx+1,0:ny+1))
        allocate(  f8(0:nx+1,0:ny+1))
!
        allocate( gp0(0:nx+1,0:ny+1))
        allocate( gp1(0:nx+1,0:ny+1))
        allocate( gp2(0:nx+1,0:ny+1))
        allocate( gp3(0:nx+1,0:ny+1))
        allocate( gp4(0:nx+1,0:ny+1))
        allocate( gp5(0:nx+1,0:ny+1))
        allocate( gp6(0:nx+1,0:ny+1))
        allocate( gp7(0:nx+1,0:ny+1))
        allocate( gp8(0:nx+1,0:ny+1))
!        
        allocate(  g0(0:nx+1,0:ny+1))
        allocate(  g1(0:nx+1,0:ny+1))
        allocate(  g2(0:nx+1,0:ny+1))
        allocate(  g3(0:nx+1,0:ny+1))
        allocate(  g4(0:nx+1,0:ny+1))
        allocate(  g5(0:nx+1,0:ny+1))
        allocate(  g6(0:nx+1,0:ny+1))
        allocate(  g7(0:nx+1,0:ny+1))
        allocate(  g8(0:nx+1,0:ny+1))
!
        allocate(u1(0:nx+2,0:ny+2))
        allocate(v1(0:nx+2,0:ny+2))
        allocate(u2(0:nx+2,0:ny+2))
        allocate(v2(0:nx+2,0:ny+2))
        allocate(psi(-1:nx+2,-1:ny+2))
        allocate(rhod1(-1:nx+2,-1:ny+2))
        allocate(rhod2(-1:nx+2,-1:ny+2))
        allocate(p(0:nx+1,0:ny+1))

        allocate(T(0:Nx+1,0:Ny+1))
        allocate(T2(0:Nx+1,0:Ny+1))
        allocate(T_prec(0:Nx+1,0:Ny+1))
        allocate(phi(0:Nx+1,0:Ny+1))
        allocate(phi_prec(0:Nx+1,0:Ny+1))
        allocate(react_tot(1:Nx,1:Ny))

        allocate(iflag(0:nx+1,0:ny+1))
        allocate(param(1:nx,1:ny))

#ifdef DEBUG
        write(6,*) "Completed subroutine alloca"
#endif

        return
        end subroutine alloca

