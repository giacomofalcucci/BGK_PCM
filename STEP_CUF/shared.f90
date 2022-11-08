module shared

 implicit none
! 
 integer, parameter           :: dp = selected_real_kind(15,307)
 integer, parameter           :: sp = kind(1.)
#ifdef SINGLEPRECISION
 integer, parameter           :: mykind = sp
#else
 integer, parameter           :: mykind = dp
#endif
 integer, parameter           :: n_th = 1
 integer, parameter           :: npop = 9
! 
! real constant
 real(kind=mykind), parameter     :: c00 = 0.0000000000000000000
 real(kind=mykind), parameter     :: c05 = 0.5000000000000000000
 real(kind=mykind), parameter     :: c10 = 1.0000000000000000000
 real(kind=mykind), parameter     :: c15 = 1.5000000000000000000
 real(kind=mykind), parameter     :: c20 = 2.0000000000000000000
 real(kind=mykind), parameter     :: c30 = 3.0000000000000000000
 real(kind=mykind), parameter     :: c45 = 4.5000000000000000000
 real(kind=mykind), parameter     :: pi = acos(-1.d0)
! real scalars 
 real(kind=mykind)                :: kvisc, alpha, T0, Tin, rho0, rho1
 real(kind=mykind)                :: eps, T_in_sol
 real(kind=mykind)                :: w_eq0,w_eq1,w_eq2,w_eq3,w_eq4,w_eq5,w_eq6,w_eq7,w_eq8
 real(kind=mykind)                :: sum1, sum2, sum3, sum4
 real(kind=mykind)                :: omega_g, omega_f, tau_g, tau_f, grav
 real(kind=mykind)                :: mass, length, temp_aver
 real(kind=mykind)                :: T_crit, T_act, T_width,f_plus,f_minus
 real(kind=mykind)            :: umv, upv, uv, uu, vv, u0, v0, irho, dsdh, dhdh
 real(kind=mykind)                :: cs, cs2, cs4
 real(kind=mykind)                :: dist, radius, x_center, y_center
 real(kind=mykind)                :: delta_T,stefan_number
 real(kind=mykind)            :: Re, Mach
! vector
 real(kind=mykind), allocatable   :: cx(:),cy(:)
 real(kind=mykind), allocatable   :: Nu(:), Nu2(:), Nu3(:)
 real(kind=mykind), allocatable   :: Tavg(:), Tavg2(:)
 real(kind=mykind), allocatable   :: gradTwall(:), gradTwall2(:)
 real(kind=mykind), allocatable   :: x(:), y(:)
! 2 indexes matrices 
 real(kind=mykind), allocatable   :: phi(:,:),phi_prec(:,:)
 real(kind=mykind), device, allocatable   :: phi_gpu(:,:),phi_prec_gpu(:,:)
 real(kind=mykind), allocatable   :: T_prec(:,:) 
 real(kind=mykind), allocatable   :: T(:,:),u(:,:), v(:,:), rho(:,:)
 real(kind=mykind), allocatable   :: T2(:,:), u2(:,:), v2(:,:), rho2(:,:)
 real(kind=mykind), device, allocatable  :: T_gpu(:,:)
 real(kind=mykind), device, allocatable  :: u_gpu(:,:)
 real(kind=mykind), device, allocatable  :: v_gpu(:,:)
 real(kind=mykind), device, allocatable  :: rho_gpu(:,:)
 real(kind=mykind), device, allocatable  :: T2_gpu(:,:), u2_gpu(:,:), v2_gpu(:,:), rho2_gpu(:,:)
 real(kind=mykind), allocatable   :: flag(:,:)
! f populations 
 real(kind=mykind), allocatable   :: f0(:,:)
 real(kind=mykind), allocatable   :: f1(:,:)
 real(kind=mykind), allocatable   :: f2(:,:)
 real(kind=mykind), allocatable   :: f3(:,:)
 real(kind=mykind), allocatable   :: f4(:,:)
 real(kind=mykind), allocatable   :: f5(:,:)
 real(kind=mykind), allocatable   :: f6(:,:)
 real(kind=mykind), allocatable   :: f7(:,:)
 real(kind=mykind), allocatable   :: f8(:,:)
! f populations  (gpu)
 real(kind=mykind), device, allocatable   :: f0_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f1_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f2_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f3_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f4_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f5_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f6_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f7_gpu(:,:)
 real(kind=mykind), device, allocatable   :: f8_gpu(:,:)
! fp populations  (post streaming)
 real(kind=mykind), allocatable   :: fp0(:,:)
 real(kind=mykind), allocatable   :: fp1(:,:)
 real(kind=mykind), allocatable   :: fp2(:,:)
 real(kind=mykind), allocatable   :: fp3(:,:)
 real(kind=mykind), allocatable   :: fp4(:,:)
 real(kind=mykind), allocatable   :: fp5(:,:)
 real(kind=mykind), allocatable   :: fp6(:,:)
 real(kind=mykind), allocatable   :: fp7(:,:)
 real(kind=mykind), allocatable   :: fp8(:,:)
! fp populations  (post streaming) gpu
 real(kind=mykind), device, allocatable   :: fp0_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp1_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp2_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp3_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp4_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp5_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp6_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp7_gpu(:,:)
 real(kind=mykind), device, allocatable   :: fp8_gpu(:,:)
! g populations 
 real(kind=mykind), allocatable   :: g0(:,:)
 real(kind=mykind), allocatable   :: g1(:,:)
 real(kind=mykind), allocatable   :: g2(:,:)
 real(kind=mykind), allocatable   :: g3(:,:)
 real(kind=mykind), allocatable   :: g4(:,:)
 real(kind=mykind), allocatable   :: g5(:,:)
 real(kind=mykind), allocatable   :: g6(:,:)
 real(kind=mykind), allocatable   :: g7(:,:)
 real(kind=mykind), allocatable   :: g8(:,:)
! g populations (gpu)
 real(kind=mykind), device, allocatable   :: g0_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g1_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g2_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g3_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g4_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g5_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g6_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g7_gpu(:,:)
 real(kind=mykind), device, allocatable   :: g8_gpu(:,:)
! gp populations  (post streaming)
 real(kind=mykind), allocatable   :: gp0(:,:)
 real(kind=mykind), allocatable   :: gp1(:,:)
 real(kind=mykind), allocatable   :: gp2(:,:)
 real(kind=mykind), allocatable   :: gp3(:,:)
 real(kind=mykind), allocatable   :: gp4(:,:)
 real(kind=mykind), allocatable   :: gp5(:,:)
 real(kind=mykind), allocatable   :: gp6(:,:)
 real(kind=mykind), allocatable   :: gp7(:,:)
 real(kind=mykind), allocatable   :: gp8(:,:)
! gp populations (post streaming) gpu
 real(kind=mykind), device, allocatable   :: gp0_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp1_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp2_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp3_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp4_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp5_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp6_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp7_gpu(:,:)
 real(kind=mykind), device, allocatable   :: gp8_gpu(:,:)
! 
 integer                      :: Nx, Ny, it, iter, it0, a, b, samplex, sampley
 integer                      :: i, j, k, m
 character*15                 :: id
 integer                      :: eastBC_st, westBC_st, northBC_st, southBC_st
 integer                      :: dump, noutput
 integer                      :: interf

!======================================================================
! OUTPUT_PARAVIEW variables
!======================================================================
  character(280)              :: Nxtot, Nytot
  character(280)              :: string, command
  character(280)              :: fn = ' '
  character(280)              :: filename_1 = ' '
  character(280)              :: filename
  character(280)              :: velocity, density, coordinates
  character(280)              :: temperature, ph_field
  real(kind=mykind)               :: sizereal
  integer                     :: sizeint
  integer                     :: ss, ss1, ss2, deltaOut, itOut
  integer                     :: offset_velocity, offset_density, offset_temp
  integer                     :: offset_coordinates, offset_phasef

end module shared
