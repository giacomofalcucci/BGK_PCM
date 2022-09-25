 subroutine dump_out 

 use shared

 open(unit=90,file='dumpFile.dat')

 write(90,*) it

 do j = 0,Ny+1
    do i = 0,Nx+1
       write(90,*) u(i,j), v(i,j), rho(i,j), (f(k,i,j) , k=0,npop-1), (feq(k,i,j) , k=0,npop-1)
       write(90,*) T(i,j), (g(k,i,j) , k=0,npop-1), (geq(k,i,j) , k=0,npop-1)
    end do
 end do

 close(90)

 end subroutine dump_out
