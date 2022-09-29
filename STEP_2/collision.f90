 subroutine collision
 
 use shared

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j,k)
!$OMP DO
 do j = 1,Ny
    do i = 1,Nx

       if (flag(i,j) .eq. 0) then

          do k = 0,npop-1
             f(k,i,j) = f(k,i,j)-omega_f*(f(k,i,j)-feq(k,i,j))
             g(k,i,j) = g(k,i,j)-omega_g*(g(k,i,j)-geq(k,i,j))
          end do

       end if

    end do
 end do
!$OMP END DO
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine collision"
#endif

 end subroutine collision
