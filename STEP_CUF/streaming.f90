 subroutine streaming

 use shared

  !$cuf kernel  do(2) <<<*,(128,4) >>>
 do j=1,Ny
    do i=1,Nx
          fp0_gpu(i,j) = f0_gpu(i  ,j  )
          fp1_gpu(i,j) = f1_gpu(i-1,j  )
          fp2_gpu(i,j) = f2_gpu(i,  j-1)
          fp3_gpu(i,j) = f3_gpu(i+1,j  )
          fp4_gpu(i,j) = f4_gpu(i,  j+1)
          fp5_gpu(i,j) = f5_gpu(i-1,j-1)
          fp6_gpu(i,j) = f6_gpu(i+1,j-1)
          fp7_gpu(i,j) = f7_gpu(i+1,j+1)
          fp8_gpu(i,j) = f8_gpu(i-1,j+1)
!          
          gp0_gpu(i,j) = g0_gpu(i  ,j  )
          gp1_gpu(i,j) = g1_gpu(i-1,j  )
          gp2_gpu(i,j) = g2_gpu(i,  j-1)
          gp3_gpu(i,j) = g3_gpu(i+1,j  )
          gp4_gpu(i,j) = g4_gpu(i,  j+1)
          gp5_gpu(i,j) = g5_gpu(i-1,j-1)
          gp6_gpu(i,j) = g6_gpu(i+1,j-1)
          gp7_gpu(i,j) = g7_gpu(i+1,j+1)
          gp8_gpu(i,j) = g8_gpu(i-1,j+1)
    end do
 end do

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine streaming"
#endif

 end subroutine streaming
