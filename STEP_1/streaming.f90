 subroutine streaming

 use shared

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j)
!!!!$OMP DO
!$OMP SECTIONS
!$OMP SECTION
 do j = 1,Ny
    do i = Nx,1,-1
       if (flag(i,j) .eq. 0) then
          f(4,i,j) = f(4,i,  j+1)
          f(8,i,j) = f(8,i-1,j+1)
          g(4,i,j) = g(4,i,  j+1)
          g(8,i,j) = g(8,i-1,j+1)
       end if
    end do
 end do
!!!!$OMP END DO
!!!!
!!!!$OMP DO
!$OMP SECTION
 do j = Ny,1,-1
    do i = 1,Nx
       if (flag(i,j) .eq. 0) then
          f(2,i,j) = f(2,i,  j-1)
          f(6,i,j) = f(6,i+1,j-1)
          g(2,i,j) = g(2,i,  j-1)
          g(6,i,j) = g(6,i+1,j-1)
       end if
    end do
 end do
!!!!$OMP END DO
!!!!
!!!!$OMP DO
!$OMP SECTION
 do j = 1,Ny
    do i = 1,Nx
       if (flag(i,j) .eq. 0) then
          f(3,i,j) = f(3,i+1,j  )
          f(7,i,j) = f(7,i+1,j+1)
          g(3,i,j) = g(3,i+1,j  )
          g(7,i,j) = g(7,i+1,j+1)
       end if
    end do
 end do
!!!!$OMP END DO
!!!!
!!!!$OMP DO
!$OMP SECTION
 do j = Ny,1,-1
    do i = Nx,1,-1
       if (flag(i,j) .eq. 0) then
          f(1,i,j) = f(1,i-1,j  )
          f(5,i,j) = f(5,i-1,j-1)
          g(1,i,j) = g(1,i-1,j  )
          g(5,i,j) = g(5,i-1,j-1)
       end if
    end do
 end do
!!!!!!$OMP END DO
!!!!!!$OMP END PARALLEL

!$OMP END SECTIONS
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine streaming"
#endif

 end subroutine streaming
