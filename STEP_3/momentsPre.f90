subroutine momentspre

 use shared

!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j,k) 
!$OMP SECTIONS
!$OMP SECTION
 do j = 1,Ny
    do i = 1,Nx
       if (flag(i,j) .eq. 0) then

          rho(i,j) = 0.d0
          T(i,j) = 0.d0
          do k = 0,npop-1
             rho(i,j) = rho(i,j) + fp(k,i,j)
             T(i,j)   = T(i,j)   + gp(k,i,j)
          end do
  
          irho = 1.d0/rho(i,j)
   
          u(i,j)    = irho*(fp(1,i,j)-fp(3,i,j)+fp(5,i,j)-fp(6,i,j)-fp(7,i,j)+fp(8,i,j))
          v(i,j)    = irho*(fp(2,i,j)-fp(4,i,j)+fp(5,i,j)+fp(6,i,j)-fp(7,i,j)-fp(8,i,j))
!
          rho2(i,j) = rho(i,j)
          T2(i,j)   = T(i,j)
          u2(i,j)    = u(i,j)
          v2(i,j)    = v(i,j)
       end if
    end do
 end do

!write(*,*) 'CHECK1 : T(nx/2,ny/2) =', T(nx/2,ny/2)

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
