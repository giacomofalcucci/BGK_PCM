 subroutine streaming

 use shared

 do j=1,Ny
    do i=1,Nx
       if (flag(i,j) .eq. 0) then
!          
          fp(0,i,j) = f(0,i  ,j  )
          fp(1,i,j) = f(1,i-1,j  )
          fp(2,i,j) = f(2,i,  j-1)
          fp(3,i,j) = f(3,i+1,j  )
          fp(4,i,j) = f(4,i,  j+1)
          fp(5,i,j) = f(5,i-1,j-1)
          fp(6,i,j) = f(6,i+1,j-1)
          fp(7,i,j) = f(7,i+1,j+1)
          fp(8,i,j) = f(8,i-1,j+1)
!          
          gp(0,i,j) = g(0,i  ,j  )
          gp(1,i,j) = g(1,i-1,j  )
          gp(2,i,j) = g(2,i,  j-1)
          gp(3,i,j) = g(3,i+1,j  )
          gp(4,i,j) = g(4,i,  j+1)
          gp(5,i,j) = g(5,i-1,j-1)
          gp(6,i,j) = g(6,i+1,j-1)
          gp(7,i,j) = g(7,i+1,j+1)
          gp(8,i,j) = g(8,i-1,j+1)
!          
       end if
    end do
 end do

!$OMP END SECTIONS
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine streaming"
#endif

 end subroutine streaming
