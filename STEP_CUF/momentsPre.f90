subroutine momentspre

 use shared

 !$cuf kernel  do(2) <<<*, (128,4) >>>
 do j = 1,Ny
    do i = 1,Nx

          rho_gpu(i,j) = fp0_gpu(i,j)+fp1_gpu(i,j)+fp2_gpu(i,j)+fp3_gpu(i,j)+fp4_gpu(i,j) &
                    +fp5_gpu(i,j)+fp6_gpu(i,j)+fp7_gpu(i,j)+fp8_gpu(i,j)
          T_gpu(i,j)   = gp0_gpu(i,j)+gp1_gpu(i,j)+gp2_gpu(i,j)+gp3_gpu(i,j)+gp4_gpu(i,j) &
                    +gp5_gpu(i,j)+gp6_gpu(i,j)+gp7_gpu(i,j)+gp8_gpu(i,j)
  
          irho = 1.d0/rho_gpu(i,j)
   
          u_gpu(i,j)    = irho*(fp1_gpu(i,j)-fp3_gpu(i,j)+fp5_gpu(i,j)-fp6_gpu(i,j)-fp7_gpu(i,j)+fp8_gpu(i,j))
          v_gpu(i,j)    = irho*(fp2_gpu(i,j)-fp4_gpu(i,j)+fp5_gpu(i,j)+fp6_gpu(i,j)-fp7_gpu(i,j)-fp8_gpu(i,j))
!
          rho2_gpu(i,j) = rho_gpu(i,j)
          T2_gpu(i,j)   = T_gpu(i,j)
          u2_gpu(i,j)   = u_gpu(i,j)
          v2_gpu(i,j)   = v_gpu(i,j)
    end do
 end do

!write(*,*) 'CHECK1 : T(nx/2,ny/2) =', T(nx/2,ny/2)

#ifdef DEBUG_GA
k = 0
 do i = 1,samplex
    gradTwall(i) = (T(k,2)-T(k,0))/2.d0 
    gradTwall2(i) = (T(k,3)-T(k,1))/2.d0 
    sum1 = 0.d0
    sum2 = 0.d0
    sum3 = 0.d0
    sum4 = 0.d0
    do j = 1,sampley+1
       sum1 = T(k,(j-1)*b)*u(k,(j-1)*b) + sum1
       sum2 = u(k,(j-1)*b) + sum2
       sum3 = T(k,(j-1)*b)*u(k,(j-1)*b)*rho(k,(j-1)*b) + sum1
       sum4 = u(k,(j-1)*b)*rho(k,(j-1)*b) + sum2
    end do
    Tavg(i) = sum1/sum2
    Tavg2(i) = sum3/sum4
    Nu(i) = -2.d0*Ny/(T0-Tavg(i))*gradTwall(i)
    Nu2(i) = -2.d0*Ny/(T0-Tavg(i))*gradTwall2(i)
    Nu3(i) = -2.d0*Ny/(T0-Tavg2(i))*gradTwall(i)
    k = k + a
 end do
write(*,*) 'CHECK2 : T(nx/2,ny/2) =', T(nx/2,ny/2)
#endif


#ifdef DEBUG_GA
 temp_aver = 0
 do j=1,Ny
   do i=1,Nx
    temp_aver = temp_aver + T(i,j)
    enddo
 enddo 
 temp_aver = temp_aver/(Nx*Ny)

write(*,*) 'DEBUG: T(nx/2,ny/2) =', T(nx/2,ny/2)
write(*,*) 'DEBUG: temp_aver =', temp_aver
#endif

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine momentspre"
#endif

end subroutine momentspre
