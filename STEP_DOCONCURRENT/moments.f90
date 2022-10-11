subroutine moments

 use shared

 do concurrent(j=1:Ny,i=1:Nx)
          rho(i,j) = f0(i,j)+f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j) & 
                    +f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)
          irho = 1.d0/rho(i,j)
   
          u(i,j)    = irho*(f1(i,j)-f3(i,j)+f5(i,j)-f6(i,j)-f7(i,j)+f8(i,j))
          v(i,j)    = irho*(f2(i,j)-f4(i,j)+f5(i,j)+f6(i,j)-f7(i,j)-f8(i,j))

          T(i,j)   = g0(i,j)+g1(i,j)+g2(i,j)+g3(i,j)+g4(i,j) & 
                    +g5(i,j)+g6(i,j)+g7(i,j)+g8(i,j)
 end do


#ifdef DEBUG_GA
!k = 0    ! Col OMP DO deve stare FUORI, con la SECTION DENTRO....!!!!! :O
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
