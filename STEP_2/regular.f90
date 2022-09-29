 subroutine regular

 use shared

 real(kind=8)   ::  pxx,pyy,pxy,fneq
 real(kind=8)   ::  Txx,Tyy,Txy,gneq

        call equilibrium

        do j=1,Ny
         do i=1,Nx
!
             !non equilibrium part of the momentum flux tensor
           pxx=    + (            - cs2) *  (fp(0,i,j) - feq(0,i,j)) & 
                   + (1           - cs2) *  (fp(1,i,j) - feq(1,i,j)) & 
                   + (            - cs2) *  (fp(2,i,j) - feq(2,i,j)) & 
                   + (1           - cs2) *  (fp(3,i,j) - feq(3,i,j)) & 
                   + (            - cs2) *  (fp(4,i,j) - feq(4,i,j)) & 
                   + (1           - cs2) *  (fp(5,i,j) - feq(5,i,j)) & 
                   + (1           - cs2) *  (fp(6,i,j) - feq(6,i,j)) & 
                   + (1           - cs2) *  (fp(7,i,j) - feq(7,i,j)) & 
                   + (1           - cs2) *  (fp(8,i,j) - feq(8,i,j))   

           pyy=    + (            - cs2) *  (fp(0,i,j) - feq(0,i,j)) & 
                   + (            - cs2) *  (fp(1,i,j) - feq(1,i,j)) & 
                   + (1           - cs2) *  (fp(2,i,j) - feq(2,i,j)) & 
                   + (            - cs2) *  (fp(3,i,j) - feq(3,i,j)) & 
                   + (1           - cs2) *  (fp(4,i,j) - feq(4,i,j)) & 
                   + (1           - cs2) *  (fp(5,i,j) - feq(5,i,j)) & 
                   + (1           - cs2) *  (fp(6,i,j) - feq(6,i,j)) & 
                   + (1           - cs2) *  (fp(7,i,j) - feq(7,i,j)) & 
                   + (1           - cs2) *  (fp(8,i,j) - feq(8,i,j))   

           pxy=    + (1                ) *  (fp(5,i,j) - feq(5,i,j)) & 
                   - (1                ) *  (fp(6,i,j) - feq(6,i,j)) & 
                   + (1                ) *  (fp(7,i,j) - feq(7,i,j)) & 
                   - (1                ) *  (fp(8,i,j) - feq(8,i,j))   
!
           Txx=    + (            - cs2) *  (gp(0,i,j) - geq(0,i,j)) &
                   + (1           - cs2) *  (gp(1,i,j) - geq(1,i,j)) &
                   + (            - cs2) *  (gp(2,i,j) - geq(2,i,j)) &
                   + (1           - cs2) *  (gp(3,i,j) - geq(3,i,j)) &
                   + (            - cs2) *  (gp(4,i,j) - geq(4,i,j)) &
                   + (1           - cs2) *  (gp(5,i,j) - geq(5,i,j)) &
                   + (1           - cs2) *  (gp(6,i,j) - geq(6,i,j)) &
                   + (1           - cs2) *  (gp(7,i,j) - geq(7,i,j)) &
                   + (1           - cs2) *  (gp(8,i,j) - geq(8,i,j))  

           Tyy=    + (            - cs2) *  (gp(0,i,j) - geq(0,i,j)) & 
                   + (            - cs2) *  (gp(1,i,j) - geq(1,i,j)) &
                   + (1           - cs2) *  (gp(2,i,j) - geq(2,i,j)) &
                   + (            - cs2) *  (gp(3,i,j) - geq(3,i,j)) &
                   + (1           - cs2) *  (gp(4,i,j) - geq(4,i,j)) &
                   + (1           - cs2) *  (gp(5,i,j) - geq(5,i,j)) &
                   + (1           - cs2) *  (gp(6,i,j) - geq(6,i,j)) &
                   + (1           - cs2) *  (gp(7,i,j) - geq(7,i,j)) &
                   + (1           - cs2) *  (gp(8,i,j) - geq(8,i,j))  

           Txy=    + (1                ) *  (gp(5,i,j) - geq(5,i,j)) &
                   - (1                ) *  (gp(6,i,j) - geq(6,i,j)) &
                   + (1                ) *  (gp(7,i,j) - geq(7,i,j)) &
                   - (1                ) *  (gp(8,i,j) - geq(8,i,j))  
!
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


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine regular"
#endif

        return 

end subroutine regular
