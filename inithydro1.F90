!--------------------------------------------------
        subroutine inithydro1
!---------------------------------------------------
! ------- modules
        use storage
        implicit double precision(a-h,o-z)

#ifdef NVFORTRAN        
        real*8 rand 
#else
        real*4 rand ! hack for intel compiler
#endif        

! Set Liquid Front at t = 0

      interf = int(nx/40)
      T_in   = T0
      u0     = Re * kvisc / length

! set everything 
        do j = 0, ny+1 !1, ny
           do i = 0, nx+1 !1, nx
              u1(i,j) = Re * kvisc / length
              v1(i,j) = 0
              u2(i,j) = Re * kvisc / length
              v2(i,j) = 0
           enddo
        enddo


! Case 1 (single bubble)
        if (icond == 1) then 
!                
! for a squared box
!           radius2 = (ny*0.16)**2
! for a rectangular box
           radius2 = (ny*0.25)**2
           do j = 0,ny+1
              do i =0,nx+1
!              
                 u1(i,j) = u0
                 v1(i,j) = v0
                 u2(i,j) = u0
                 v2(i,j) = v0
!            
                 if (((i-nx/2)**2+(j-ny/2)**2).lt.radius2) then
                    rhod1(i,j)= 1.939027909286642
!                    rhod1(i,j)=  &
!                    2.410d0*(1.d0 + 0.01d0* (rand(0) - 0.5d0) * 2.d0) 
                 else
                    rhod1(i,j)=0.1572145251524053
!                    rhod1(i,j)=0.1250d0
                 endif
              enddo
           enddo
        endif

! Case 2 (many bubbles)
        if (icond == 2) then 
!                
! modified to take into account the possibility of localized peturbation
           do j = 1, ny
              do i = 1, nx
                 u1(i,j) = u0
                 v1(i,j) = v0
                 u2(i,j) = u0
                 v2(i,j) = v0
              enddo
           enddo

           do j = 1,ny
              do i =1,nx
                 rhod1(i,j)=rhoin*(1.d0+0.01d0*(rand(0)-0.5d0)*2.d0)
              enddo
           enddo
        endif

! Case 3 (two bubbles colliding)
        if (icond == 3) then
           radius2 = (ny*0.15)**2
!
! left bubble                
           do j = 0,ny+1
              do i = 0,nx/2
                 if (((i-nx/4)**2+(j-ny/2-60)**2).lt.radius2) then
                    rhod1(i,j)= 1.939027909286642
                    u1(i,j) = u0
                    v1(i,j) = v0
                    u2(i,j) = u0
                    v2(i,j) = v0
!
                 else
                    rhod1(i,j)=0.1572145251524053
                 endif
              enddo
           enddo
!
! right bubble                
           do j = 0,ny+1
              do i = nx/2+1, nx+1
                 if (((i-3*nx/4)**2+(j-ny/2)**2).lt.radius2) then
                    rhod1(i,j)= 1.939027909286642
!
                    u1(i,j) = -u0
                    v1(i,j) = -v0
                    u2(i,j) = -u0
                    v2(i,j) = -v0
!
                 else
                    rhod1(i,j)=0.1572145251524053
                 endif
              enddo
           enddo
        endif

! Case 4 (PCM - MELTING and SOLIDIFICATION)
        if (icond == 4) then
           do j=0, ny+1
            do i =0, interf
               T(i,j)   = T0
               rho(i,j) = rho1
             enddo
             do i=interf+1, nx+1
               T(i,j)   = T_in_sol
               rho(i,j) = rho0
             enddo
           enddo
          
           write(*,*) '***********************************************'
           write(*,*) 'alpha= ' , alpha
           write(*,*) 'kvisc= ' , kvisc
           write(*,*) 'u0= ' , u0
           write(*,*) 'T0= ' , T0
           write(*,*) '***********************************************'
          
           if (u0 .ge. (sqrt(1./3.))/3.d0) then
              write(*,*) '=============================='
              write(*,*) 'WARNING! Mach number too high!'
              write(*,*) '=============================='
           end if

           do j = 0, ny+1
              T(0,j) = T0
           end do
          
           do i = 0, nx+1
              u1(i,Ny+1) = 0.d0
              v1(i,Ny+1) = 0.d0
              u2(i,0) = 0.d0
              v2(i,0) = 0.d0
           end do

           do j = 1, ny
              u1(0,j)    = u0
              v1(0,j)    = v0
              u2(nx+1,j) = 0.d0
              v2(nx+1,j) = 0.d0
           end do

           do j=0, ny+1
            do i=0, nx+1
              phi(i,j)      = -1.0d0    !"-1" := SOLID, "+1" := LIQUID
              phi_prec(i,j) = -1.0d0
            enddo
           enddo
          
           do j=0, ny+1
            do i=0,interf
               phi(i,j)      = 1.0d0    !"-1" := SOLID, "+1" := LIQUID
               phi_prec(i,j) = 1.0d0
            enddo
           enddo


        endif

        if ((icond .gt. 4).or.(icond.lt.1)) then 
           write(6,*) "ERROR: option non supported"
           stop
        endif

        return
        end subroutine inithydro1

