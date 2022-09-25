 subroutine moments2

 use shared

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j)
!$OMP DO
 do j = 1,Ny
    do i = 1,Nx
             rho2(i,j) = rho(i,j) 
             T2(i,j)   = T(i,j) 
   
             u2(i,j)    = u(i,j)
             v2(i,j)    = v(i,j) 
    end do
 end do
!$OMP END DO
!$OMP END PARALLEL
!!

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine moments2"
#endif

 end subroutine moments2
