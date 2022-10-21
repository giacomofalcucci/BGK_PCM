 subroutine collision
 
 use shared

 real(kind=mykind)   ::  feq00, feq01, feq02, feq03, feq04, feq05, feq06, feq07, feq08
 real(kind=mykind)   ::  geq00, geq01, geq02, geq03, geq04, geq05, geq06, geq07, geq08
 ! temporary scalar
 real(kind=mykind)   ::  ft00, ft01, ft02, ft03, ft04, ft05, ft06, ft07, ft08
 real(kind=mykind)   ::  gt00, gt01, gt02, gt03, gt04, gt05, gt06, gt07, gt08

 real(kind=mykind)   ::  pxx,pyy,pxy,fneq
 real(kind=mykind)   ::  Txx,Tyy,Txy,gneq

 !$cuf kernel  do(2) <<<*,(128,4)>>>
 do j = 1,Ny
    do i = 1,Nx

          uu = u_gpu(i,j)*u_gpu(i,j)
          vv = v_gpu(i,j)*v_gpu(i,j)
          uv = 1.5d0*(uu+vv)
          upv = u_gpu(i,j)+v_gpu(i,j)
          umv = u_gpu(i,j)-v_gpu(i,j)

          feq00 = w_eq0*rho_gpu(i,j)*(1.d0-uv)
          feq01 = w_eq1*rho_gpu(i,j)*(1.d0+3.d0*u_gpu(i,j)+4.5d0*uu-uv)
          feq02 = w_eq2*rho_gpu(i,j)*(1.d0+3.d0*v_gpu(i,j)+4.5d0*vv-uv)
          feq03 = w_eq3*rho_gpu(i,j)*(1.d0-3.d0*u_gpu(i,j)+4.5d0*uu-uv)
          feq04 = w_eq4*rho_gpu(i,j)*(1.d0-3.d0*v_gpu(i,j)+4.5d0*vv-uv)
          feq05 = w_eq5*rho_gpu(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          feq06 = w_eq6*rho_gpu(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          feq07 = w_eq7*rho_gpu(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          feq08 = w_eq8*rho_gpu(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
          
          geq00 = w_eq0*T_gpu(i,j)*(1.d0-uv)
          geq01 = w_eq1*T_gpu(i,j)*(1.d0+3.d0*u_gpu(i,j)+4.5d0*uu-uv)
          geq02 = w_eq2*T_gpu(i,j)*(1.d0+3.d0*v_gpu(i,j)+4.5d0*vv-uv)
          geq03 = w_eq3*T_gpu(i,j)*(1.d0-3.d0*u_gpu(i,j)+4.5d0*uu-uv)
          geq04 = w_eq4*T_gpu(i,j)*(1.d0-3.d0*v_gpu(i,j)+4.5d0*vv-uv)
          geq05 = w_eq5*T_gpu(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          geq06 = w_eq6*T_gpu(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          geq07 = w_eq7*T_gpu(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          geq08 = w_eq8*T_gpu(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
!
             !non equilibrium part of the momentum flux tensor
           pxx=    + (            - cs2) *  (fp0_gpu(i,j) - feq00) & 
                   + (1           - cs2) *  (fp1_gpu(i,j) - feq01) & 
                   + (            - cs2) *  (fp2_gpu(i,j) - feq02) & 
                   + (1           - cs2) *  (fp3_gpu(i,j) - feq03) & 
                   + (            - cs2) *  (fp4_gpu(i,j) - feq04) & 
                   + (1           - cs2) *  (fp5_gpu(i,j) - feq05) & 
                   + (1           - cs2) *  (fp6_gpu(i,j) - feq06) & 
                   + (1           - cs2) *  (fp7_gpu(i,j) - feq07) & 
                   + (1           - cs2) *  (fp8_gpu(i,j) - feq08)   

           pyy=    + (            - cs2) *  (fp0_gpu(i,j) - feq00) & 
                   + (            - cs2) *  (fp1_gpu(i,j) - feq01) & 
                   + (1           - cs2) *  (fp2_gpu(i,j) - feq02) & 
                   + (            - cs2) *  (fp3_gpu(i,j) - feq03) & 
                   + (1           - cs2) *  (fp4_gpu(i,j) - feq04) & 
                   + (1           - cs2) *  (fp5_gpu(i,j) - feq05) & 
                   + (1           - cs2) *  (fp6_gpu(i,j) - feq06) & 
                   + (1           - cs2) *  (fp7_gpu(i,j) - feq07) & 
                   + (1           - cs2) *  (fp8_gpu(i,j) - feq08)   

           pxy=    + (1                ) *  (fp5_gpu(i,j) - feq05) & 
                   - (1                ) *  (fp6_gpu(i,j) - feq06) & 
                   + (1                ) *  (fp7_gpu(i,j) - feq07) & 
                   - (1                ) *  (fp8_gpu(i,j) - feq08)   
!
           Txx=    + (            - cs2) *  (gp0_gpu(i,j) - geq00) &
                   + (1           - cs2) *  (gp1_gpu(i,j) - geq01) &
                   + (            - cs2) *  (gp2_gpu(i,j) - geq02) &
                   + (1           - cs2) *  (gp3_gpu(i,j) - geq03) &
                   + (            - cs2) *  (gp4_gpu(i,j) - geq04) &
                   + (1           - cs2) *  (gp5_gpu(i,j) - geq05) &
                   + (1           - cs2) *  (gp6_gpu(i,j) - geq06) &
                   + (1           - cs2) *  (gp7_gpu(i,j) - geq07) &
                   + (1           - cs2) *  (gp8_gpu(i,j) - geq08)  

           Tyy=    + (            - cs2) *  (gp0_gpu(i,j) - geq00) & 
                   + (            - cs2) *  (gp1_gpu(i,j) - geq01) &
                   + (1           - cs2) *  (gp2_gpu(i,j) - geq02) &
                   + (            - cs2) *  (gp3_gpu(i,j) - geq03) &
                   + (1           - cs2) *  (gp4_gpu(i,j) - geq04) &
                   + (1           - cs2) *  (gp5_gpu(i,j) - geq05) &
                   + (1           - cs2) *  (gp6_gpu(i,j) - geq06) &
                   + (1           - cs2) *  (gp7_gpu(i,j) - geq07) &
                   + (1           - cs2) *  (gp8_gpu(i,j) - geq08)  

           Txy=    + (1                ) *  (gp5_gpu(i,j) - geq05) &
                   - (1                ) *  (gp6_gpu(i,j) - geq06) &
                   + (1                ) *  (gp7_gpu(i,j) - geq07) &
                   - (1                ) *  (gp8_gpu(i,j) - geq08)  
!
             ft0= feq00 + ((0.5d0 * w_eq0) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (  - cs2)*pyy            )
             ft1= feq01 + ((0.5d0 * w_eq1) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (  - cs2)*pyy            )
             ft2= feq02 + ((0.5d0 * w_eq2) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (1 - cs2)*pyy            )  
             ft3= feq03 + ((0.5d0 * w_eq3) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (  - cs2)*pyy            ) 
             ft4= feq04 + ((0.5d0 * w_eq4) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (1 - cs2)*pyy            )  
             ft5= feq05 + ((0.5d0 * w_eq5) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy +2.0d0*pxy )
             ft6= feq06 + ((0.5d0 * w_eq6) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy -2.0d0*pxy )
             ft7= feq07 + ((0.5d0 * w_eq7) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy +2.0d0*pxy )
             ft8= feq08 + ((0.5d0 * w_eq8) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy -2.0d0*pxy )
!
             gt0= geq00 + ((0.5d0 * w_eq0) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (  - cs2)*Tyy            ) 
             gt1= geq01 + ((0.5d0 * w_eq1) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (  - cs2)*Tyy            )
             gt2= geq02 + ((0.5d0 * w_eq2) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (1 - cs2)*Tyy            )  
             gt3= geq03 + ((0.5d0 * w_eq3) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (  - cs2)*Tyy            ) 
             gt4= geq04 + ((0.5d0 * w_eq4) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (1 - cs2)*Tyy            )  
             gt5= geq05 + ((0.5d0 * w_eq5) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy +2.0d0*Txy )
             gt6= geq06 + ((0.5d0 * w_eq6) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy -2.0d0*Txy )
             gt7= geq07 + ((0.5d0 * w_eq7) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy +2.0d0*Txy )
             gt8= geq08 + ((0.5d0 * w_eq8) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy -2.0d0*Txy )

          f0_gpu(i,j) = ft0-omega_f*(ft0-feq00)
          f1_gpu(i,j) = ft1-omega_f*(ft1-feq01)
          f2_gpu(i,j) = ft2-omega_f*(ft2-feq02)
          f3_gpu(i,j) = ft3-omega_f*(ft3-feq03)
          f4_gpu(i,j) = ft4-omega_f*(ft4-feq04)
          f5_gpu(i,j) = ft5-omega_f*(ft5-feq05)
          f6_gpu(i,j) = ft6-omega_f*(ft6-feq06)
          f7_gpu(i,j) = ft7-omega_f*(ft7-feq07)
          f8_gpu(i,j) = ft8-omega_f*(ft8-feq08)
!
          g0_gpu(i,j) = gt0-omega_g*(gt0-geq00)
          g1_gpu(i,j) = gt1-omega_g*(gt1-geq01)
          g2_gpu(i,j) = gt2-omega_g*(gt2-geq02)
          g3_gpu(i,j) = gt3-omega_g*(gt3-geq03)
          g4_gpu(i,j) = gt4-omega_g*(gt4-geq04)
          g5_gpu(i,j) = gt5-omega_g*(gt5-geq05)
          g6_gpu(i,j) = gt6-omega_g*(gt6-geq06)
          g7_gpu(i,j) = gt7-omega_g*(gt7-geq07)
          g8_gpu(i,j) = gt8-omega_g*(gt8-geq08)
    end do
 end do
!$acc end kernels

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine collision"
#endif

 end subroutine collision
