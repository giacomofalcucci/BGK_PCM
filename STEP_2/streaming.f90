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
          ff(i,j,4) = ff(i,  j+1,4)
          ff(i,j,8) = ff(i-1,j+1,8)
          gg(i,j,4) = gg(i,  j+1,4)
          gg(i,j,8) = gg(i-1,j+1,8)
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
          ff(i,j,2) = ff(i,  j-1,2)
          ff(i,j,6) = ff(i+1,j-1,6)
          gg(i,j,2) = gg(i,  j-1,2)
          gg(i,j,6) = gg(i+1,j-1,6)
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
          ff(i,j,3) = ff(i+1,j  ,3)
          ff(i,j,7) = ff(i+1,j+1,7)
          gg(i,j,3) = gg(i+1,j  ,3)
          gg(i,j,7) = gg(i+1,j+1,7)
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
          ff(i,j,1) = ff(i-1,j  ,1)
          ff(i,j,5) = ff(i-1,j-1,5)
          gg(i,j,1) = gg(i-1,j  ,1)
          gg(i,j,5) = gg(i-1,j-1,5)
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
