 subroutine field_check

 use shared
 implicit none

 real(kind=mykind)  :: mass_aver, phi_aver

 mass_aver = 0.0d0
 temp_aver = 0.0d0
 phi_aver  = 0.0d0

 do j = 1,Ny
    do i = 1,Nx
       mass_aver = mass_aver + rho(i,j)
       temp_aver = temp_aver + T(i,j)
       phi_aver  = phi_aver  + phi(i,j)
    end do
 end do

 mass_aver = mass_aver /(Nx*Ny)
 temp_aver = temp_aver /(Nx*Ny)
 phi_aver  = phi_aver /(Nx*Ny)

 write(64,1002) it,mass_aver,temp_aver,phi_aver
 flush(64)


! formats...
1001    format("CHECK: Mass, Temp, Phi",3(e14.6,1x))
1002    format(i8,3(e14.6,1x))



#ifdef DEBUG
        write(6,1001) mass_aver, temp_aver, phi_aver
        write(6,*) "DEBUG: Completed subroutine field_check"
#endif

 end subroutine field_check
