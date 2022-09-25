subroutine moments

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
             rho(i,j) = rho(i,j) + f(k,i,j)
             T(i,j) = T(i,j) + g(k,i,j)
          end do
  
          irho = 1.d0/rho(i,j)
   
          u(i,j)    = irho*(f(1,i,j)-f(3,i,j)+f(5,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j))
          v(i,j)    = irho*(f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j))
 
       end if
    end do
 end do

!write(*,*) 'CHECK1 : T(nx/2,ny/2) =', T(nx/2,ny/2)


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

!write(*,*) 'CHECK2 : T(nx/2,ny/2) =', T(nx/2,ny/2)
!

!$OMP END SECTIONS
!$OMP END PARALLEL


 temp_aver = 0
 do j=1,Ny
   do i=1,Nx
    temp_aver = temp_aver + T(i,j)
    enddo
 enddo 

!write(*,*) 'CHECK3 : T(nx/2,ny/2) =', T(nx/2,ny/2)
!write(*,*) 'CHECK4: temp_aver =', temp_aver

 temp_aver = temp_aver/(Nx*Ny)

!write(*,*) 'CHECK5: temp_aver =', temp_aver

end subroutine moments
