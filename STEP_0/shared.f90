module shared

 implicit none
 integer, parameter           :: dp = selected_real_kind(15,307)
 integer, parameter           :: n_th = 1
 real(kind=dp), parameter     :: pi = acos(-1.d0)
 integer, parameter           :: npop = 9
 real(kind=dp)                :: kvisc, alpha, T0, Tin, rho0, rho1
 real(kind=dp)                :: eps, T_in_sol
 real(kind=dp)                :: w_eq(0:8), sum1, sum2, sum3, sum4
 real(kind=dp)                :: cx(0:8),cy(0:8)
 real(kind=dp)                :: omega_g, omega_f, tau_g, tau_f, grav
 real(kind=dp), allocatable   :: Nu(:), Nu2(:), Nu3(:)
 real(kind=dp), allocatable   :: Tavg(:), Tavg2(:), gradTwall(:), gradTwall2(:)
 real(kind=dp), allocatable   :: phi(:,:),phi_prec(:,:),react_tot(:,:)
 real(kind=dp), allocatable   :: f(:,:,:), feq(:,:,:), ftemp(:,:,:)
 real(kind=dp), allocatable   :: g(:,:,:), geq(:,:,:), T(:,:), T_prec(:,:)
 real(kind=dp), allocatable   :: u(:,:), v(:,:), rho(:,:)
 real(kind=dp), allocatable   :: u2(:,:), v2(:,:), rho2(:,:),T2(:,:)
 real(kind=dp), allocatable   :: force_x(:,:), force_y(:,:)
 real(kind=dp), allocatable   :: react_plus(:,:),react_minus(:,:)
 real(kind=dp), allocatable   :: k_plus(:,:),k_minus(:,:)
 real(kind=dp)                :: mass, length, temp_aver
 real(kind=dp)                :: T_crit, T_act, T_width,f_plus,f_minus
 integer                      :: Nx, Ny, it, iter, it0, a, b, samplex, sampley
 real(kind=dp), allocatable   :: x(:), y(:)
 real(kind=dp), allocatable   :: flag(:,:)
 real(kind=dp)                :: umv, upv, uv, uu, vv, u0, v0, irho, dsdh, dhdh
 real(kind=dp)                :: Re, Mach
 integer                      :: i, j, k, m
 character                    :: id(50)
 integer                      :: eastBC_st, westBC_st, northBC_st, southBC_st
 integer                      :: dump, noutput
 real(kind=dp)                :: cs, cs2, cs4
 real(kind=dp)                :: dist, radius, x_center, y_center
 real(kind=dp)                :: delta_T,stefan_number
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
  real(kind=dp)               :: sizereal
  integer                     :: sizeint
  integer                     :: ss, ss1, ss2, deltaOut, itOut
  integer                     :: offset_velocity, offset_density, offset_temp
  integer                     :: offset_coordinates, offset_phasef

end module shared
