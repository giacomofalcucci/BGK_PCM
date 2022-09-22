!
        module storage
!
        implicit double precision(a-h,o-z)
!
        integer, parameter:: sp=kind(1.0)
        integer, parameter:: dp=selected_real_kind(2*precision(1.0_sp))
!
        integer :: nx
        integer :: ny
        integer, parameter :: npop = 9
!
        character*5 fileout
        integer::  nrhout
        logical iforce,iobst
!/phys/ 
        real(dp) u0,v0,uf,fom,length  
        real(dp) f_guo
        real(dp) kvisc, alpha, T0, Tin, rho0, rho1, eps, T_in_sol, grav
        real(dp) Re, Nus, Nus_comp, Nus_comp_2, Gras, Pr, Ray, Fou
!/constants/
        real(dp) cs2,cs22,cssq,rhoin,omega_f,tau_f,fpois,den,visc
        real(dp) omega_g,tau_g,stefan_number
        real(dp) w0,w1,w2,w4,w5,w8,gnn,gnnn,rhoaver,dinvrho
        real(dp) rhopsi,dt,dx,dump,c1_2,c2_2,c4_2,c5_2,c8_2 
        real(dp), parameter :: cte04 = (4.d0/ 9.d0)
        real(dp), parameter :: cte09 = (1.d0/ 9.d0)
        real(dp), parameter :: cte36 = (1.d0/36.d0)
!/count/ 
        integer istep,nout,ndiag,nsteps,nobst,icond
!/arrays to allocate/
        real(dp), dimension (:,:), allocatable :: f0
        real(dp), dimension (:,:), allocatable :: f1
        real(dp), dimension (:,:), allocatable :: f2
        real(dp), dimension (:,:), allocatable :: f3
        real(dp), dimension (:,:), allocatable :: f4
        real(dp), dimension (:,:), allocatable :: f5
        real(dp), dimension (:,:), allocatable :: f6
        real(dp), dimension (:,:), allocatable :: f7
        real(dp), dimension (:,:), allocatable :: f8
!
        real(dp), dimension (:,:), allocatable :: fp0
        real(dp), dimension (:,:), allocatable :: fp1
        real(dp), dimension (:,:), allocatable :: fp2
        real(dp), dimension (:,:), allocatable :: fp3
        real(dp), dimension (:,:), allocatable :: fp4
        real(dp), dimension (:,:), allocatable :: fp5
        real(dp), dimension (:,:), allocatable :: fp6
        real(dp), dimension (:,:), allocatable :: fp7
        real(dp), dimension (:,:), allocatable :: fp8
!
        real(dp), dimension (:,:), allocatable :: g0
        real(dp), dimension (:,:), allocatable :: g1
        real(dp), dimension (:,:), allocatable :: g2
        real(dp), dimension (:,:), allocatable :: g3
        real(dp), dimension (:,:), allocatable :: g4
        real(dp), dimension (:,:), allocatable :: g5
        real(dp), dimension (:,:), allocatable :: g6
        real(dp), dimension (:,:), allocatable :: g7
        real(dp), dimension (:,:), allocatable :: g8
!
        real(dp), dimension (:,:), allocatable :: gp0
        real(dp), dimension (:,:), allocatable :: gp1
        real(dp), dimension (:,:), allocatable :: gp2
        real(dp), dimension (:,:), allocatable :: gp3
        real(dp), dimension (:,:), allocatable :: gp4
        real(dp), dimension (:,:), allocatable :: gp5
        real(dp), dimension (:,:), allocatable :: gp6
        real(dp), dimension (:,:), allocatable :: gp7
        real(dp), dimension (:,:), allocatable :: gp8
!
        real(dp), dimension (:,:), allocatable :: u1
        real(dp), dimension (:,:), allocatable :: v1
        real(dp), dimension (:,:), allocatable :: u2
        real(dp), dimension (:,:), allocatable :: v2
        real(dp), dimension (:,:), allocatable :: psi
        real(dp), dimension (:,:), allocatable :: rhod1
        real(dp), dimension (:,:), allocatable :: rhod2
        real(dp), dimension (:,:), allocatable :: p
        real(dp), dimension (:,:), allocatable :: T
        real(dp), dimension (:,:), allocatable :: T_prec
        real(dp), dimension (:,:), allocatable :: phi
        real(dp), dimension (:,:), allocatable :: phi_prec
        real(dp), dimension (:,:), allocatable :: react_tot
        real(dp), dimension (:,:), allocatable :: param
!       
        real(dp), dimension (:), allocatable   :: Tavg
        real(dp), dimension (:), allocatable   :: Tavg2
        real(dp), dimension (:), allocatable   :: gradTwall
        real(dp), dimension (:), allocatable   :: gradTwall2
!
        integer, dimension (:,:), allocatable ::  iflag
!

!/arrays/
        real(dp), dimension (0:npop-1) ::     w
        real(dp), dimension (0:npop-1) ::     u_ci
!                
!        
        end module  storage

