 subroutine collision
 
 use shared

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j,k)
!$OMP DO
 do j = 1,Ny
    do i = 1,Nx

       if (flag(i,j) .eq. 0) then

          do k = 0,npop-1
             ff(i,j,k) = ff(i,j,k)-omega_f*(ff(i,j,k)-ffeq(i,j,k))
             gg(i,j,k) = gg(i,j,k)-omega_g*(gg(i,j,k)-ggeq(i,j,k))
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
