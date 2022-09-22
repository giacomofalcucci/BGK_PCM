!----------------------------------------------------------
        subroutine laplace
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        p=0.d0

        do j = 1,ny
          do i = 1,nx
            pid  =rhod1(i,j)*cs2
            pnid =0.5d0*(gnn+gnnn)*cs2*psi(i,j)*psi(i,j)
!0.5d0*0.7365d0/3.d0*gnn*psi(i,j)*psi(i,j)+
!                                                          &
!                                                          0.2635d0/3.d0*gnnn*psi(i,j)*psi(i,j)
!                                                          !*psi(i,j)*psi(i,j)
            p(i,j)=pid+pnid
             if(p(i,j).gt.pmax)then
                pmax = p(i,j)
             endif
             if(p(i,j).lt.pmin)then
                pmin = p(i,j)
             endif

           write(65,*)i,j,p(i,j)

           enddo

           write(65,'(bn)')

         enddo

         do i=1,nx
            write(66,*)i,p(i,ny/2)
         enddo


            rhoaver=0.d0

        rhomax = -100000000.d0
        rhomin = 100000000.d0
        do i = 1, nx
          do j = 1, ny
            rhomax = max(rhomax,rhod1(i,j))
            rhomin = min(rhomin,rhod1(i,j))

         enddo
        enddo

         write(*,*) 'rho_min=',rhomin,'rho_max=',rhomax
         write(*,*) '......................................'

         write(*,*) '..computing mean rho...'

         do j=1,ny
          do i=1,nx
              rhoaver = rhoaver + rhod1(i,j)
          enddo
         enddo

         write(*,*)'A_1 =',gnn+gnnn
         write(*,*)'A_2 =',gnn+1.5d0*gnnn


              rhoaver = rhoaver/(nx*ny)

         write(*,*) 'rhoaver =', rhoaver

         radius = ((rhoaver*(nx*ny)-rhomin*(nx*ny))/(3.1415926536d0* &
     &             (rhomax-rhomin)))**(0.5d0)

         write(*,*) 'Bubble radius (Case 1)', radius

          surf_tens = (p(nx/2,ny/2)-p(nx/16,ny/16))*radius

          write(*,*) 'Laplace test: gamma = ', surf_tens

          return
          end subroutine laplace
 
