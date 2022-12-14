 subroutine regular

 use shared

 real(kind=8)   ::  pxx,pyy,pxy,fneq
 real(kind=8)   ::  Txx,Tyy,Txy,gneq

 call equilibrium

        do j=1,Ny
         do i=1,Nx
           pxx=0.0d0
           pyy=0.0d0
           pxy=0.0d0
!
           Txx=0.0d0
           Tyy=0.0d0
           Txy=0.0d0
           do k=0,npop-1
             fneq = f(k,i,j) - feq(k,i,j)
             gneq = g(k,i,j) - geq(k,i,j)
             !non equilibrium part of the momentum flux tensor
             pxx=pxx + (cx(k)*cx(k) - cs2) *  fneq
             pyy=pyy + (cy(k)*cy(k) - cs2) *  fneq
             pxy=pxy + (cx(k)*cy(k))       *  fneq
!
             Txx=Txx + (cx(k)*cx(k) - cs2) *  gneq
             Tyy=Tyy + (cy(k)*cy(k) - cs2) *  gneq
             Txy=Txy + (cx(k)*cy(k))       *  gneq
           enddo

           do k=0,npop-1
             f(k,i,j)= feq(k,i,j) + ((0.5d0 * w_eq(k)) / (cs2*cs2)) * &
                         ((cx(k)*cx(k) - cs2) * pxx + &
                         ( cy(k)*cy(k) - cs2) * pyy +  &
                           2.0d0 *cx(k)*cy(k)  * pxy )
!
             g(k,i,j)= geq(k,i,j) + ((0.5d0 * w_eq(k)) / (cs2*cs2)) * &
                         ((cx(k)*cx(k) - cs2) * Txx + &
                         ( cy(k)*cy(k) - cs2) * Tyy +  &
                           2.0d0 *cx(k)*cy(k)  * Txy )

           enddo


         enddo
        enddo

        return 





end subroutine regular
