subroutine moments

 use shared

 real(kind=mykind)  :: temp_T, temp_u, temp_v 


 !$cuf kernel  do(2) <<<*, (128,4) >>>
 do j = 1,Ny
    do i = 1,Nx

          rho_gpu(i,j) = f0_gpu(i,j)+f1_gpu(i,j)+f2_gpu(i,j) & 
                        +f3_gpu(i,j)+f4_gpu(i,j)+f5_gpu(i,j) & 
                        +f6_gpu(i,j)+f7_gpu(i,j)+f8_gpu(i,j)

          irho = 1.d0/rho_gpu(i,j)

          temp_T   = g0_gpu(i,j)+g1_gpu(i,j)+g2_gpu(i,j)+g3_gpu(i,j)+g4_gpu(i,j) & 
                    +g5_gpu(i,j)+g6_gpu(i,j)+g7_gpu(i,j)+g8_gpu(i,j)
  
   
          temp_u    = irho*(f1_gpu(i,j)-f3_gpu(i,j)+f5_gpu(i,j)-f6_gpu(i,j)-f7_gpu(i,j)+f8_gpu(i,j))
          temp_v    = irho*(f2_gpu(i,j)-f4_gpu(i,j)+f5_gpu(i,j)+f6_gpu(i,j)-f7_gpu(i,j)-f8_gpu(i,j))
 
          T_gpu(i,j)     = (temp_T+T2_gpu(i,j))*0.50d0  ! THIS, I have to check....
   
          u_gpu(i,j)     = (temp_u+u2_gpu(i,j))*0.50d0
          v_gpu(i,j)     = (temp_v+v2_gpu(i,j))*0.50d0
    end do
 end do
 !$acc end kernels

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
!
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
        write(6,*) "DEBUG: Completed subroutine moments"
#endif

end subroutine moments
