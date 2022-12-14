 subroutine dump_in

 use shared

 open(unit=90,file='dumpFile.dat')

 read(90,*) it

 do j = 0,Ny+1
    do i = 0,Nx+1
       read(90,*) u(i,j), v(i,j), rho(i,j), (f(k,i,j) , k=0,npop-1), (feq(k,i,j) , k=0,npop-1)
       read(90,*) T(i,j), (g(k,i,j) , k=0,npop-1), (geq(k,i,j) , k=0,npop-1)
    end do
 end do

 close(90)

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine dump_in"
#endif

 end subroutine dump_in
