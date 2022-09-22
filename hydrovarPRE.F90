!---------------------------------------------
        subroutine hydrovarPRE
!----------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

! Calculation of velocities and pseudopotential
!$OMP PARALLEL DEFAULT(NONE)                                &
!$OMP PRIVATE(i,j,rho1)                                     &
!$OMP SHARED(nx,ny,rhopsi)                                  &
!$OMP SHARED(u1,v1,rhod1,u2,v2,psi)                         &
!$OMP SHARED(f0)                                            &
!$OMP SHARED(fp1,fp2,fp3,fp4,fp5,fp6,fp7,fp8)               
!$OMP DO
!        
!$acc kernels
!$acc loop independent collapse(2)
        do j = 0, ny+1
!DIR$ IVDEP
           do i = 0, nx+1
              rhod1(i,j) =  f0(i,j) + &
                           fp1(i,j) + &
                           fp2(i,j) + &
                           fp3(i,j) + &
                           fp4(i,j) + &
                           fp5(i,j) + &
                           fp6(i,j) + &
                           fp7(i,j) + &
                           fp8(i,j)
              rho1 = 1.d0 / rhod1(i,j)        
              u1(i,j) = ( fp1(i,j) - fp3(i,j) + fp5(i,j)            &
                        - fp6(i,j) - fp7(i,j) + fp8(i,j) ) * rho1 
              v1(i,j) = ( fp5(i,j) + fp2(i,j) + fp6(i,j)            &
                        - fp7(i,j) - fp4(i,j) - fp8(i,j) ) * rho1
!                
              u2(i,j) = u1(i,j)
              v2(i,j) = v1(i,j)
!
! (from previous force subroutine)
              psi(i,j) =rhopsi * (1.d0 -exp(- rhod1(i,j) / rhopsi))

           enddo
        enddo
!$acc end kernels
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "Completed subroutine hydrovarPRE"
#endif

        end subroutine hydrovarPRE
