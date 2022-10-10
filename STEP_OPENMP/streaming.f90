 subroutine streaming

 use shared


 !$OMP PARALLEL DEFAULT(NONE)  &
 !$OMP PRIVATE(i,j)  &
 !$OMP SHARED(Nx,Ny)  &
 !$OMP SHARED(fp0,fp1,fp2,fp3,fp4)  &
 !$OMP SHARED(fp5,fp6,fp7,fp8    )  &
 !$OMP SHARED(f0,f1,f2,f3,f4)  &
 !$OMP SHARED(f5,f6,f7,f8    )  &
 !$OMP SHARED(gp0,gp1,gp2,gp3,gp4)  &
 !$OMP SHARED(gp5,gp6,gp7,gp8    )  &
 !$OMP SHARED(g0,g1,g2,g3,g4)  &
 !$OMP SHARED(g5,g6,g7,g8    )  
 !$OMP DO
 do j=1,Ny
    do i=1,Nx
          fp0(i,j) = f0(i  ,j  )
          fp1(i,j) = f1(i-1,j  )
          fp2(i,j) = f2(i,  j-1)
          fp3(i,j) = f3(i+1,j  )
          fp4(i,j) = f4(i,  j+1)
          fp5(i,j) = f5(i-1,j-1)
          fp6(i,j) = f6(i+1,j-1)
          fp7(i,j) = f7(i+1,j+1)
          fp8(i,j) = f8(i-1,j+1)
!          
          gp0(i,j) = g0(i  ,j  )
          gp1(i,j) = g1(i-1,j  )
          gp2(i,j) = g2(i,  j-1)
          gp3(i,j) = g3(i+1,j  )
          gp4(i,j) = g4(i,  j+1)
          gp5(i,j) = g5(i-1,j-1)
          gp6(i,j) = g6(i+1,j-1)
          gp7(i,j) = g7(i+1,j+1)
          gp8(i,j) = g8(i-1,j+1)
    end do
 end do
!$OMP END PARALLEL

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine streaming"
#endif

 end subroutine streaming
