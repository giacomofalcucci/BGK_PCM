 subroutine mass_check

 use shared
 implicit none

 real(kind=mykind)  :: mass_aver

 mass = 0.d0
 mass_aver = 0.0d0

 do j = 1,Ny
    do i = 1,Nx
       mass = mass + rho(i,j)
       mass_aver = mass_aver + rho(i,j)
    end do
 end do

 mass = (Nx*Ny*rho0-mass)*100.d0/(Nx*Ny*rho0)
 mass_aver = mass_aver / (Nx*Ny)

 temp_aver = 0
 do j=1,Ny
   do i=1,Nx
    temp_aver = temp_aver + T(i,j)
    enddo
 enddo

temp_aver = temp_aver/(Nx*Ny)

 write(64,*) it,mass_aver,temp_aver 

 write(6,1001) mass_aver, temp_aver

 ! formats...
1001    format("CHECK: Mass_aver, Temp_aver",2(e14.6,1x))


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine mass_check"
#endif

 end subroutine mass_check
