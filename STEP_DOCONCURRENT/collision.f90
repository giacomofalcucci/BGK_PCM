 subroutine collision
 
 use shared
 implicit none

 real(kind=mykind)   ::  feq00, feq01, feq02, feq03, feq04, feq05, feq06, feq07, feq08
 real(kind=mykind)   ::  geq00, geq01, geq02, geq03, geq04, geq05, geq06, geq07, geq08
 ! temporary scalar
 real(kind=mykind)   ::  ft00, ft01, ft02, ft03, ft04, ft05, ft06, ft07, ft08
 real(kind=mykind)   ::  gt00, gt01, gt02, gt03, gt04, gt05, gt06, gt07, gt08

 real(kind=mykind)   ::  pxx,pyy,pxy,fneq
 real(kind=mykind)   ::  Txx,Tyy,Txy,gneq

 do concurrent (j = 1:Ny,i = 1:Nx)

          uu = u(i,j)*u(i,j)
          vv = v(i,j)*v(i,j)
          uv = c15*(uu+vv)
          upv = u(i,j)+v(i,j)
          umv = u(i,j)-v(i,j)

          feq00 = w_eq0*rho(i,j)*(c10-uv)
          feq01 = w_eq1*rho(i,j)*(c10+c30*u(i,j)+c45*uu-uv)
          feq02 = w_eq2*rho(i,j)*(c10+c30*v(i,j)+c45*vv-uv)
          feq03 = w_eq3*rho(i,j)*(c10-c30*u(i,j)+c45*uu-uv)
          feq04 = w_eq4*rho(i,j)*(c10-c30*v(i,j)+c45*vv-uv)
          feq05 = w_eq5*rho(i,j)*(c10+c30*upv+c45*upv*upv-uv)
          feq06 = w_eq6*rho(i,j)*(c10-c30*umv+c45*umv*umv-uv)
          feq07 = w_eq7*rho(i,j)*(c10-c30*upv+c45*upv*upv-uv)
          feq08 = w_eq8*rho(i,j)*(c10+c30*umv+c45*umv*umv-uv)
          
          geq00 = w_eq0*T(i,j)*(c10-uv)
          geq01 = w_eq1*T(i,j)*(c10+c30*u(i,j)+c45*uu-uv)
          geq02 = w_eq2*T(i,j)*(c10+c30*v(i,j)+c45*vv-uv)
          geq03 = w_eq3*T(i,j)*(c10-c30*u(i,j)+c45*uu-uv)
          geq04 = w_eq4*T(i,j)*(c10-c30*v(i,j)+c45*vv-uv)
          geq05 = w_eq5*T(i,j)*(c10+c30*upv+c45*upv*upv-uv)
          geq06 = w_eq6*T(i,j)*(c10-c30*umv+c45*umv*umv-uv)
          geq07 = w_eq7*T(i,j)*(c10-c30*upv+c45*upv*upv-uv)
          geq08 = w_eq8*T(i,j)*(c10+c30*umv+c45*umv*umv-uv)
!
             !non equilibrium part of the momentum flux tensor
           pxx=    + (            - cs2) *  (fp0(i,j) - feq00) & 
                   + (c10           - cs2) *  (fp1(i,j) - feq01) & 
                   + (            - cs2) *  (fp2(i,j) - feq02) & 
                   + (c10           - cs2) *  (fp3(i,j) - feq03) & 
                   + (            - cs2) *  (fp4(i,j) - feq04) & 
                   + (c10           - cs2) *  (fp5(i,j) - feq05) & 
                   + (c10           - cs2) *  (fp6(i,j) - feq06) & 
                   + (c10           - cs2) *  (fp7(i,j) - feq07) & 
                   + (c10           - cs2) *  (fp8(i,j) - feq08)   

           pyy=    + (            - cs2) *  (fp0(i,j) - feq00) & 
                   + (            - cs2) *  (fp1(i,j) - feq01) & 
                   + (c10           - cs2) *  (fp2(i,j) - feq02) & 
                   + (            - cs2) *  (fp3(i,j) - feq03) & 
                   + (c10           - cs2) *  (fp4(i,j) - feq04) & 
                   + (c10           - cs2) *  (fp5(i,j) - feq05) & 
                   + (c10           - cs2) *  (fp6(i,j) - feq06) & 
                   + (c10           - cs2) *  (fp7(i,j) - feq07) & 
                   + (c10           - cs2) *  (fp8(i,j) - feq08)   

           pxy=    + (c10                ) *  (fp5(i,j) - feq05) & 
                   - (c10                ) *  (fp6(i,j) - feq06) & 
                   + (c10                ) *  (fp7(i,j) - feq07) & 
                   - (c10                ) *  (fp8(i,j) - feq08)   
!
           Txx=    + (            - cs2) *  (gp0(i,j) - geq00) &
                   + (c10           - cs2) *  (gp1(i,j) - geq01) &
                   + (            - cs2) *  (gp2(i,j) - geq02) &
                   + (c10           - cs2) *  (gp3(i,j) - geq03) &
                   + (            - cs2) *  (gp4(i,j) - geq04) &
                   + (c10           - cs2) *  (gp5(i,j) - geq05) &
                   + (c10           - cs2) *  (gp6(i,j) - geq06) &
                   + (c10           - cs2) *  (gp7(i,j) - geq07) &
                   + (c10           - cs2) *  (gp8(i,j) - geq08)  

           Tyy=    + (            - cs2) *  (gp0(i,j) - geq00) & 
                   + (            - cs2) *  (gp1(i,j) - geq01) &
                   + (c10           - cs2) *  (gp2(i,j) - geq02) &
                   + (            - cs2) *  (gp3(i,j) - geq03) &
                   + (c10           - cs2) *  (gp4(i,j) - geq04) &
                   + (c10           - cs2) *  (gp5(i,j) - geq05) &
                   + (c10           - cs2) *  (gp6(i,j) - geq06) &
                   + (c10           - cs2) *  (gp7(i,j) - geq07) &
                   + (c10           - cs2) *  (gp8(i,j) - geq08)  

           Txy=    + (c10                ) *  (gp5(i,j) - geq05) &
                   - (c10                ) *  (gp6(i,j) - geq06) &
                   + (c10                ) *  (gp7(i,j) - geq07) &
                   - (c10                ) *  (gp8(i,j) - geq08)  
!
             ft00= feq00 + ((c05 * w_eq0) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (  - cs2)*pyy            )
             ft01= feq01 + ((c05 * w_eq1) / (cs2*cs2)) * &
                         ((c10 - cs2)*pxx + (  - cs2)*pyy            )
             ft02= feq02 + ((c05 * w_eq2) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (c10 - cs2)*pyy            )  
             ft03= feq03 + ((c05 * w_eq3) / (cs2*cs2)) * &
                         ((c10 - cs2)*pxx + (  - cs2)*pyy            ) 
             ft04= feq04 + ((c05 * w_eq4) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (c10 - cs2)*pyy            )  
             ft05= feq05 + ((c05 * w_eq5) / (cs2*cs2)) * &
                         ((c10 - cs2)*pxx + (c10 - cs2)*pyy +c20*pxy )
             ft06= feq06 + ((c05 * w_eq6) / (cs2*cs2)) * &
                         ((c10 - cs2)*pxx + (c10 - cs2)*pyy -c20*pxy )
             ft07= feq07 + ((c05 * w_eq7) / (cs2*cs2)) * &
                         ((c10 - cs2)*pxx + (c10 - cs2)*pyy +c20*pxy )
             ft08= feq08 + ((c05 * w_eq8) / (cs2*cs2)) * &
                         ((c10 - cs2)*pxx + (c10 - cs2)*pyy -c20*pxy )
!
             gt00= geq00 + ((c05 * w_eq0) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (  - cs2)*Tyy            ) 
             gt01= geq01 + ((c05 * w_eq1) / (cs2*cs2)) * &
                         ((c10 - cs2)*Txx + (  - cs2)*Tyy            )
             gt02= geq02 + ((c05 * w_eq2) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (c10 - cs2)*Tyy            )  
             gt03= geq03 + ((c05 * w_eq3) / (cs2*cs2)) * &
                         ((c10 - cs2)*Txx + (  - cs2)*Tyy            ) 
             gt04= geq04 + ((c05 * w_eq4) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (c10 - cs2)*Tyy            )  
             gt05= geq05 + ((c05 * w_eq5) / (cs2*cs2)) * &
                         ((c10 - cs2)*Txx + (c10 - cs2)*Tyy +c20*Txy )
             gt06= geq06 + ((c05 * w_eq6) / (cs2*cs2)) * &
                         ((c10 - cs2)*Txx + (c10 - cs2)*Tyy -c20*Txy )
             gt07= geq07 + ((c05 * w_eq7) / (cs2*cs2)) * &
                         ((c10 - cs2)*Txx + (c10 - cs2)*Tyy +c20*Txy )
             gt08= geq08 + ((c05 * w_eq8) / (cs2*cs2)) * &
                         ((c10 - cs2)*Txx + (c10 - cs2)*Tyy -c20*Txy )

          f0(i,j) = ft00-omega_f*(ft00-feq00)
          f1(i,j) = ft01-omega_f*(ft01-feq01)
          f2(i,j) = ft02-omega_f*(ft02-feq02)
          f3(i,j) = ft03-omega_f*(ft03-feq03)
          f4(i,j) = ft04-omega_f*(ft04-feq04)
          f5(i,j) = ft05-omega_f*(ft05-feq05)
          f6(i,j) = ft06-omega_f*(ft06-feq06)
          f7(i,j) = ft07-omega_f*(ft07-feq07)
          f8(i,j) = ft08-omega_f*(ft08-feq08)
!
          g0(i,j) = gt00-omega_g*(gt00-geq00)
          g1(i,j) = gt01-omega_g*(gt01-geq01)
          g2(i,j) = gt02-omega_g*(gt02-geq02)
          g3(i,j) = gt03-omega_g*(gt03-geq03)
          g4(i,j) = gt04-omega_g*(gt04-geq04)
          g5(i,j) = gt05-omega_g*(gt05-geq05)
          g6(i,j) = gt06-omega_g*(gt06-geq06)
          g7(i,j) = gt07-omega_g*(gt07-geq07)
          g8(i,j) = gt08-omega_g*(gt08-geq08)
 end do

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine collision"
#endif

 end subroutine collision
