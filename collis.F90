!----------------------------------------------------------
        subroutine collis
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)


!$OMP PARALLEL DEFAULT(NONE)                                &
!$OMP PRIVATE(i,j)                                          &
!$OMP PRIVATE(fnnx,fnny,fnnnx,fnnny,f2x,f2y)                &
!$OMP PRIVATE(usq,vsq,sumsq,sumsq2,u22,v22,ui,vi,uv,rhoij)  &
!$OMP SHARED(nx,ny,omega)                                   &
!$OMP SHARED(gnnn,w1,gnn,w2,cs2,cs22,cssq)                  &
!$OMP SHARED(psi,u1,v1,rhod1,u2,v2,rho1)                    &
!$OMP SHARED(f0,f1,f2,f3,f4,f5,f6,f7,f8)                    &
!$OMP SHARED(g0,g1,g2,g3,g4,g5,g6,g7,g8)                    &
!$OMP SHARED(fp0,fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8)           &
!$OMP SHARED(gp0,gp1,gp2,gp3,gp4,gp5,gp6,gp7,gp8)           &
!$OMP SHARED(feq0,feq1,feq2,feq3,feq4,feq5,feq6,feq7,feq8)  &
!$OMP SHARED(geq0,geq1,geq2,geq3,geq4,geq5,geq6,geq7,geq8)
!$OMP DO
!$acc kernels
!$acc loop independent collapse(2)
        do j = 1, ny
!DIR$ IVDEP
           do i = 1, nx

! forcing
! non local Lennard-Jones (may put a density functional, not density itself)
              fnnx  = psi(i,j)*(psi(i+1,j)-psi(i-1,j))
              fnny  = psi(i,j)*(psi(i,j+1)-psi(i,j-1))
              fnnnx = psi(i,j)*((psi(i+1,j-1) + psi(i+1,j+1)    &
                               - psi(i-1,j-1) - psi(i-1,j+1)))
              fnnny = psi(i,j)*((psi(i+1,j+1) + psi(i-1,j+1)    &
                               - psi(i-1,j-1) - psi(i+1,j-1)))

! interactions
              f2x =-(gnnn * fnnx * w1 + gnnn * fnnnx * w2)- &  
                    (gnn * fnnx*cte09+ gnn * fnnnx *cte36)  
              f2y =-(gnnn * fnny * w1 + gnnn * fnnny * w2)- & 
                    (gnn * fnny*cte09+ gnn * fnnny *cte36)

! SHIFT Equilibrium
              u1(i,j)=u1(i,j)+f2x/(omega*rhod1(i,j))
              v1(i,j)=v1(i,j)+f2y/(omega*rhod1(i,j))

! compute equilibrium              
              usq = u1(i,j) * u1(i,j) 
              vsq = v1(i,j) * v1(i,j)
              sumsq = (usq + vsq) / cs22
              sumsq2 = sumsq * (1.0d0 - cs2) / cs2
              u22 = usq / cssq 
              v22 = vsq / cssq
              ui = u1(i,j) / cs2
              vi = v1(i,j) / cs2
              uv = ui * vi
              rhoij = rhod1(i,j)

              feq0 = cte04*rhoij*(1.0d0 - sumsq)

              feq1 = cte09*rhoij*(1.0d0 - sumsq + u22 + ui)
              feq2 = cte09*rhoij*(1.0d0 - sumsq + v22 + vi)
              feq3 = cte09*rhoij*(1.0d0 - sumsq + u22 - ui)
              feq4 = cte09*rhoij*(1.0d0 - sumsq + v22 - vi)

              feq5 = cte36*rhoij*(1.0d0 + sumsq2 +ui+vi+uv)
              feq6 = cte36*rhoij*(1.0d0 + sumsq2 -ui+vi-uv)
              feq7 = cte36*rhoij*(1.0d0 + sumsq2 -ui-vi+uv)
              feq8 = cte36*rhoij*(1.0d0 + sumsq2 +ui-vi-uv)

! compute correction              
              f0(i,j) =  f0(i,j) * (1.0d0 - omega) + omega * feq0
!              
              f1(i,j) = fp1(i,j) * (1.0d0 - omega) + omega * feq1
              f2(i,j) = fp2(i,j) * (1.0d0 - omega) + omega * feq2
              f3(i,j) = fp3(i,j) * (1.0d0 - omega) + omega * feq3
              f4(i,j) = fp4(i,j) * (1.0d0 - omega) + omega * feq4
!              
              f5(i,j) = fp5(i,j) * (1.0d0 - omega) + omega * feq5
              f6(i,j) = fp6(i,j) * (1.0d0 - omega) + omega * feq6
              f7(i,j) = fp7(i,j) * (1.0d0 - omega) + omega * feq7
              f8(i,j) = fp8(i,j) * (1.0d0 - omega) + omega * feq8

! moved from hydrovarPOST
              rhod1(i,j) = f0(i,j) + &
                           f1(i,j) + &
                           f2(i,j) + &
                           f3(i,j) + &
                           f4(i,j) + &
                           f5(i,j) + &
                           f6(i,j) + &
                           f7(i,j) + &
                           f8(i,j)
              rho1 = 1.d0 / rhod1(i,j)        
              u1(i,j) = ( f1(i,j) - f3(i,j) + f5(i,j)            &
                        - f6(i,j) - f7(i,j) + f8(i,j) ) * rho1 
              v1(i,j) = ( f5(i,j) + f2(i,j) + f6(i,j)            &
                        - f7(i,j) - f4(i,j) - f8(i,j) ) * rho1

! computing mean velocity (pre/post collision)
              u1(i,j)=(u1(i,j)+u2(i,j))*0.5d0
              v1(i,j)=(v1(i,j)+v2(i,j))*0.5d0

           enddo
        enddo
!$OMP END PARALLEL
!$acc end kernels




#ifdef DEBUG
        write(6,*) "Completed subroutine collis"
#endif
        end subroutine collis

