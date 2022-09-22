!----------------------------------------------------------
        subroutine out1d(frce)
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

! transverse profiles
        do j = 1,ny
           write(11,*) j,rhod1(nx/2,j),u1(nx/2,j),v1(nx/2,j)
        enddo
        write(11,'(bn)')

        end

!----------------------------------------------------------
        subroutine out2d(frce)
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)
        
        do j = 1, ny
           do i = 1, nx
              write(51,99) i,j,rhod1(i,j),u1(i,j),v1(i,j)
           enddo
           write(51,'(bn)')
        enddo

        write(51,'(bn)')
        write(51,'(bn)')

 99     format(2I6,3(1x,e13.6))

        end

!----------------------------------------------------------
        subroutine energy
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

! Energy

        g_eff=gnn+gnnn

        en_pot=0.d0
        en_cin=0.d0

        do j = 1, ny
           do i = 1, nx
              en_pot=en_pot+g_eff*psi(i,j)*psi(i,j)*cs2	!/2.d0
              en_cin=en_cin+rhod1(i,j)*cs2
           enddo
        enddo

        en_pot_1=0.d0
        en_pot_2=0.d0

!        write(6,*)'i valori dei c_2 per la subroutine energy sono'
!        write(6,*) c1_2,c2_2,c4_2,c5_2,c8_2


         do j = 1, ny
          do i = 1, nx

! computing forces to find total potential energy

           en_pot_1=en_pot_1+0.5*gnn*psi(i,j)* &
     & (w(1)*c1_2*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     & w(5)*c2_2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2=en_pot_2+0.5d0*gnnn*psi(i,j)* &
     &  (w1*c1_2*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     &   w2*c2_2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+&
     &   w4*c4_2*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+&
     &   w5*c5_2*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+&
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+&
     &  w8*c8_2*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))


         enddo
        enddo

         write(86,98) istep,en_pot_1,en_pot_2,en_pot_1+en_pot_2,en_cin,&
     &                en_pot,(en_pot_1+en_pot_2)/en_cin

98      format(1I6,1x,6(1x,e13.6)) 

! potential energy without "C_i^2" (Prof. Succi 11 - 03 - 2007)

        en_pot_1_bis=0.d0
        en_pot_2_bis=0.d0

         do j = 1, ny
          do i = 1, nx


           en_pot_1_bis=en_pot_1_bis+0.5*gnn*psi(i,j)*&
     &  (4./9.)*psi(i,j)+&
     & (w(1)*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+&
     &  w(5)*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_bis=en_pot_2_bis+0.5*gnnn*psi(i,j)*&
     &  (w0*psi(i,j)+&
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+&
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+&
     &   w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+&
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+&
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+&
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))


         enddo
        enddo

         write(85,96) istep,en_pot_1_bis,en_pot_2_bis, &
     &                en_pot_1_bis+en_pot_2_bis,en_pot

96      format(1I6,1x,4(1x,e13.6))


!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
! Diagnostica a meno di G1 e G2
!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]


        en_pot_ter=0.d0

!        do j = 1, ny
!         do i = 1, nx
!            en_pot_ter=en_pot_ter+psi(i,j)*psi(i,j)*cs2   !/2.d0
!         enddo
!        enddo



        en_pot_1_ter=0.d0
        en_pot_2_ter=0.d0

        en_pot_1_bulk=0.d0
        en_pot_2_bulk=0.d0


!         write(*,*)'i w valgono,',w(0),w(1),w(5)


         do j = 1, ny
          do i = 1, nx

           en_pot_1_ter=en_pot_1_ter+0.5d0*psi(i,j)* &
     & ((4.d0/9.d0)*psi(i,j)+ &
     & w(1)*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     & w(5)*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_ter=en_pot_2_ter+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1))+ &
     &   w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+ &
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+ &
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+ &
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))



           en_pot_1_bulk=en_pot_1_bulk+0.5d0*psi(i,j)* &
     & ((4.d0/9.d0)*psi(i,j)+ &
     &  w(1)*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &  w(5)*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))                        ! SHAN-CHEN !!!!!!!!!!!!

          en_pot_2_bulk=en_pot_2_bulk+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w2*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w4*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w5*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)+ &
     &            psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w8*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))


         enddo
        enddo

97       format(1I6,1x,5(1x,e13.6))

         en_2_1=0.d0
         en_2_2=0.d0
         en_2_1_bulk=0.d0
         en_2_2_bulk=0.d0

         do j = 1, ny
          do i = 1, nx

          en_2_1=en_2_1+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i+1,j)+psi(i,j+1)+psi(i-1,j)+psi(i,j-1))+ &
     &   w2*(psi(i+1,j+1)+psi(i-1,j+1)+psi(i-1,j-1)+psi(i+1,j-1)))


          en_2_2=en_2_2+0.5d0*psi(i,j)* &
     &   (w4*(psi(i+2,j)+psi(i-2,j)+psi(i,j+2)+psi(i,j-2))+ &
     &   w5*(psi(i+2,j+1)+psi(i+2,j-1)+psi(i+1,j+2)+psi(i-1,j+2)+ &
     &            psi(i-2,j+1)+psi(i-2,j-1)+psi(i-1,j-2)+psi(i+1,j-2))+ &
     &   w8*(psi(i+2,j+2)+psi(i-2,j+2)+psi(i-2,j-2)+psi(i+2,j-2)))



          en_2_1_bulk=en_2_1_bulk+0.5d0*psi(i,j)* &
     &  (w0*psi(i,j)+ &
     &   w1*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w2*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))

          en_2_2_bulk=en_2_2_bulk+0.5d0*psi(i,j)* &
     &  (w4*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w5*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)+ &
     &            psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j))+ &
     &   w8*(psi(i,j)+psi(i,j)+psi(i,j)+psi(i,j)))


         enddo
        enddo

         write(81,91) istep,en_pot_1_ter,en_pot_2_bulk, &
     &                en_2_1,en_2_2, &
     &                en_2_1_bulk,en_2_2_bulk,en_cin

         write(89,*) istep,gnn*en_pot_1_ter+gnnn*(en_2_1+en_2_2)


91         format(1I6,1x,7(1x,e13.6))

        

!]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]

        nliq=0

        do i=1,nx
          do j=1,ny

            if(rhod1(i,j).gt.(0.7d0)) then
              nliq=nliq+1
            endif

          enddo
        enddo 


!  number of bubbles...

 
        write(84,*)istep,nliq

        return
        end

!----------------------------------------------------------
        subroutine out0d
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        densit = 0.0d0

        do j= 1, ny
           do i = 1, nx
              densit = densit + f0(i,j) +                 &
     &                   f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+ &
     &                   f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)
           enddo
        enddo
        densit = densit / dfloat(nx*ny) 

        do i = 1, nx
           do j= 1, ny
              param(i,j) = f0(i,j) +                      &
     &                   f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+ &
     &                   f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)
           enddo
        enddo
!         
	do i = 1, nx
	   do j= 1, ny
	      param(i,j) = (param(i,j)-rhoaver)**2.d0
	   enddo
	enddo
        paramord=0.d0
	do i = 1, nx
	   do j= 1, ny
	      paramord = paramord+param(i,j)
	   enddo
	enddo
        paramord = paramord / dfloat(nx*ny)
        paramord = dsqrt(paramord) / rhoaver

	umoy = 0.0d0
	vmoy = 0.0d0

	do j = 1, ny
	   do i = 1, nx
	      umoy = umoy + u1(i,j)
	      vmoy = vmoy + v1(i,j)
	   enddo
	enddo
	
	umoy = umoy / dfloat(nx*ny)
	vmoy = vmoy / dfloat(nx*ny)

! Tu uncomment for debugging
!        write(*,*)'Debug'  
!        write(*,'(6e18.10)') dfloat(istep),densit,umoy,vmoy

        write(92,'(I6,3(1x,e13.6))')istep,densit,umoy,vmoy

!        rhomax = -100000000.d0
!        rhomin = 100000000.d0
!        do i = 1, nx
!           do j = 1, ny
!              rhomax = max(rhomax,rhod1(i,j))
!              rhomin = min(rhomin,rhod1(i,j)) 
!           enddo
!        enddo

        densit = densit / dfloat(nx*ny)  ! bulk 

!        print *,' rhomax rhomin'
!        write(*,'(6e18.10)')densit,rhomax,rhomin
!        print*,'============================================='  

!!       p=0.d0
!!
!!       do j = 1,ny
!!          do i = 1,nx
!!            pid  =rhod1(i,j)*cs2
!!            pnid =0.5d0*(gnn+gnnn)*cs2*psi(i,j)*psi(i,j)     !0.5d0*0.7365d0/3.d0*gnn*psi(i,j)*psi(i,j)+
!!                                                          &  0.2635d0/3.d0*gnnn*psi(i,j)*psi(i,j) !*psi(i,j)*psi(i,j)
!!            p(i,j)=pid+pnid
!!             if(p(i,j).gt.pmax)then
!!                pmax = p(i,j)
!!             endif
!!             if(p(i,j).lt.pmin)then
!!                pmin = p(i,j)
!!             endif
!!           enddo
!!       enddo

        return
        end subroutine out0d

