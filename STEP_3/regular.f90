 subroutine regular

 use shared

 real(kind=mykind)   ::  pxx,pyy,pxy,fneq
 real(kind=mykind)   ::  Txx,Tyy,Txy,gneq
 real(kind=mykind)   ::  feq0(0:npop-1)
 real(kind=mykind)   ::  geq0(0:npop-1)

 do j = 1,Ny
    do i = 1,Nx

          uu = u(i,j)*u(i,j)
          vv = v(i,j)*v(i,j)
          uv = 1.5d0*(uu+vv)
          upv = u(i,j)+v(i,j)
          umv = u(i,j)-v(i,j)

          feq0(0) = w_eq(0)*rho(i,j)*(1.d0-uv)
          feq0(1) = w_eq(1)*rho(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          feq0(2) = w_eq(2)*rho(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          feq0(3) = w_eq(3)*rho(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          feq0(4) = w_eq(4)*rho(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          feq0(5) = w_eq(5)*rho(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          feq0(6) = w_eq(6)*rho(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          feq0(7) = w_eq(7)*rho(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          feq0(8) = w_eq(8)*rho(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
          
          geq0(0) = w_eq(0)*T(i,j)*(1.d0-uv)
          geq0(1) = w_eq(1)*T(i,j)*(1.d0+3.d0*u(i,j)+4.5d0*uu-uv)
          geq0(2) = w_eq(2)*T(i,j)*(1.d0+3.d0*v(i,j)+4.5d0*vv-uv)
          geq0(3) = w_eq(3)*T(i,j)*(1.d0-3.d0*u(i,j)+4.5d0*uu-uv)
          geq0(4) = w_eq(4)*T(i,j)*(1.d0-3.d0*v(i,j)+4.5d0*vv-uv)
          geq0(5) = w_eq(5)*T(i,j)*(1.d0+3.d0*upv+4.5d0*upv*upv-uv)
          geq0(6) = w_eq(6)*T(i,j)*(1.d0-3.d0*umv+4.5d0*umv*umv-uv)
          geq0(7) = w_eq(7)*T(i,j)*(1.d0-3.d0*upv+4.5d0*upv*upv-uv)
          geq0(8) = w_eq(8)*T(i,j)*(1.d0+3.d0*umv+4.5d0*umv*umv-uv)
!
             !non equilibrium part of the momentum flux tensor
           pxx=    + (            - cs2) *  (fp(0,i,j) - feq0(0)) & 
                   + (1           - cs2) *  (fp(1,i,j) - feq0(1)) & 
                   + (            - cs2) *  (fp(2,i,j) - feq0(2)) & 
                   + (1           - cs2) *  (fp(3,i,j) - feq0(3)) & 
                   + (            - cs2) *  (fp(4,i,j) - feq0(4)) & 
                   + (1           - cs2) *  (fp(5,i,j) - feq0(5)) & 
                   + (1           - cs2) *  (fp(6,i,j) - feq0(6)) & 
                   + (1           - cs2) *  (fp(7,i,j) - feq0(7)) & 
                   + (1           - cs2) *  (fp(8,i,j) - feq0(8))   

           pyy=    + (            - cs2) *  (fp(0,i,j) - feq0(0)) & 
                   + (            - cs2) *  (fp(1,i,j) - feq0(1)) & 
                   + (1           - cs2) *  (fp(2,i,j) - feq0(2)) & 
                   + (            - cs2) *  (fp(3,i,j) - feq0(3)) & 
                   + (1           - cs2) *  (fp(4,i,j) - feq0(4)) & 
                   + (1           - cs2) *  (fp(5,i,j) - feq0(5)) & 
                   + (1           - cs2) *  (fp(6,i,j) - feq0(6)) & 
                   + (1           - cs2) *  (fp(7,i,j) - feq0(7)) & 
                   + (1           - cs2) *  (fp(8,i,j) - feq0(8))   

           pxy=    + (1                ) *  (fp(5,i,j) - feq0(5)) & 
                   - (1                ) *  (fp(6,i,j) - feq0(6)) & 
                   + (1                ) *  (fp(7,i,j) - feq0(7)) & 
                   - (1                ) *  (fp(8,i,j) - feq0(8))   
!
           Txx=    + (            - cs2) *  (gp(0,i,j) - geq0(0)) &
                   + (1           - cs2) *  (gp(1,i,j) - geq0(1)) &
                   + (            - cs2) *  (gp(2,i,j) - geq0(2)) &
                   + (1           - cs2) *  (gp(3,i,j) - geq0(3)) &
                   + (            - cs2) *  (gp(4,i,j) - geq0(4)) &
                   + (1           - cs2) *  (gp(5,i,j) - geq0(5)) &
                   + (1           - cs2) *  (gp(6,i,j) - geq0(6)) &
                   + (1           - cs2) *  (gp(7,i,j) - geq0(7)) &
                   + (1           - cs2) *  (gp(8,i,j) - geq0(8))  

           Tyy=    + (            - cs2) *  (gp(0,i,j) - geq0(0)) & 
                   + (            - cs2) *  (gp(1,i,j) - geq0(1)) &
                   + (1           - cs2) *  (gp(2,i,j) - geq0(2)) &
                   + (            - cs2) *  (gp(3,i,j) - geq0(3)) &
                   + (1           - cs2) *  (gp(4,i,j) - geq0(4)) &
                   + (1           - cs2) *  (gp(5,i,j) - geq0(5)) &
                   + (1           - cs2) *  (gp(6,i,j) - geq0(6)) &
                   + (1           - cs2) *  (gp(7,i,j) - geq0(7)) &
                   + (1           - cs2) *  (gp(8,i,j) - geq0(8))  

           Txy=    + (1                ) *  (gp(5,i,j) - geq0(5)) &
                   - (1                ) *  (gp(6,i,j) - geq0(6)) &
                   + (1                ) *  (gp(7,i,j) - geq0(7)) &
                   - (1                ) *  (gp(8,i,j) - geq0(8))  
!
             f(0,i,j)= feq0(0) + ((0.5d0 * w_eq(0)) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (  - cs2)*pyy            )
             f(1,i,j)= feq0(1) + ((0.5d0 * w_eq(1)) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (  - cs2)*pyy            )
             f(2,i,j)= feq0(2) + ((0.5d0 * w_eq(2)) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (1 - cs2)*pyy            )  
             f(3,i,j)= feq0(3) + ((0.5d0 * w_eq(3)) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (  - cs2)*pyy            ) 
             f(4,i,j)= feq0(4) + ((0.5d0 * w_eq(4)) / (cs2*cs2)) * &
                         ((  - cs2)*pxx + (1 - cs2)*pyy            )  
             f(5,i,j)= feq0(5) + ((0.5d0 * w_eq(5)) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy +2.0d0*pxy )
             f(6,i,j)= feq0(6) + ((0.5d0 * w_eq(6)) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy -2.0d0*pxy )
             f(7,i,j)= feq0(7) + ((0.5d0 * w_eq(7)) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy +2.0d0*pxy )
             f(8,i,j)= feq0(8) + ((0.5d0 * w_eq(8)) / (cs2*cs2)) * &
                         ((1 - cs2)*pxx + (1 - cs2)*pyy -2.0d0*pxy )
!
             g(0,i,j)= geq0(0) + ((0.5d0 * w_eq(0)) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (  - cs2)*Tyy            ) 
             g(1,i,j)= geq0(1) + ((0.5d0 * w_eq(1)) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (  - cs2)*Tyy            )
             g(2,i,j)= geq0(2) + ((0.5d0 * w_eq(2)) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (1 - cs2)*Tyy            )  
             g(3,i,j)= geq0(3) + ((0.5d0 * w_eq(3)) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (  - cs2)*Tyy            ) 
             g(4,i,j)= geq0(4) + ((0.5d0 * w_eq(4)) / (cs2*cs2)) * &
                         ((  - cs2)*Txx + (1 - cs2)*Tyy            )  
             g(5,i,j)= geq0(5) + ((0.5d0 * w_eq(5)) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy +2.0d0*Txy )
             g(6,i,j)= geq0(6) + ((0.5d0 * w_eq(6)) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy -2.0d0*Txy )
             g(7,i,j)= geq0(7) + ((0.5d0 * w_eq(7)) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy +2.0d0*Txy )
             g(8,i,j)= geq0(8) + ((0.5d0 * w_eq(8)) / (cs2*cs2)) * &
                         ((1 - cs2)*Txx + (1 - cs2)*Tyy -2.0d0*Txy )

         enddo
        enddo


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine regular"
#endif

        return 

end subroutine regular
