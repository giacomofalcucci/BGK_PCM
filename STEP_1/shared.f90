module shared

 implicit none
! 
 integer, parameter           :: dp = selected_real_kind(15,307)
 integer, parameter           :: mykind = dp
 integer, parameter           :: n_th = 1
 integer, parameter           :: npop = 9
! 
! real scalars 
 real(kind=mykind), parameter     :: pi = acos(-1.d0)
 real(kind=mykind)                :: kvisc, alpha, T0, Tin, rho0, rho1
 real(kind=mykind)                :: eps, T_in_sol
 real(kind=mykind)                :: w_eq(0:8), sum1, sum2, sum3, sum4
 real(kind=mykind)                :: cx(0:8),cy(0:8)
 real(kind=mykind)                :: omega_g, omega_f, tau_g, tau_f, grav
 real(kind=mykind)                :: mass, length, temp_aver
 real(kind=mykind)                :: T_crit, T_act, T_width,f_plus,f_minus
 real(kind=mykind)            :: umv, upv, uv, uu, vv, u0, v0, irho, dsdh, dhdh
 real(kind=mykind)                :: cs, cs2, cs4
 real(kind=mykind)                :: dist, radius, x_center, y_center
 real(kind=mykind)                :: delta_T,stefan_number
 real(kind=mykind)            :: Re, Mach
! vector
 real(kind=mykind), allocatable   :: Nu(:), Nu2(:), Nu3(:)
 real(kind=mykind), allocatable   :: Tavg(:), Tavg2(:)
 real(kind=mykind), allocatable   :: gradTwall(:), gradTwall2(:)
 real(kind=mykind), allocatable   :: x(:), y(:)
! 2 indexes matrices 
 real(kind=mykind), allocatable   :: phi(:,:),phi_prec(:,:),react_tot(:,:)
 real(kind=mykind), allocatable   :: react_plus(:,:),react_minus(:,:)
 real(kind=mykind), allocatable   :: T(:,:), T_prec(:,:), T2(:,:)
 real(kind=mykind), allocatable   :: u(:,:), v(:,:), rho(:,:)
 real(kind=mykind), allocatable   :: u2(:,:), v2(:,:), rho2(:,:)
 real(kind=mykind), allocatable   :: force_x(:,:), force_y(:,:)
 real(kind=mykind), allocatable   :: k_plus(:,:),k_minus(:,:)
 real(kind=mykind), allocatable   :: flag(:,:)
! 3 indexes matrices 
 real(kind=mykind), allocatable   :: f(:,:,:), feq(:,:,:), ftemp(:,:,:)
 real(kind=mykind), allocatable   :: g(:,:,:), geq(:,:,:)
 integer                      :: Nx, Ny, it, iter, it0, a, b, samplex, sampley
 integer                      :: i, j, k, m
 character                    :: id(50)
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
