 subroutine dump_out 

 use shared

 open(unit=90,file='dumpFile.dat')

 write(90,*) it

 do j = 0,Ny+1
    do i = 0,Nx+1
       write(90,*) u(i,j), v(i,j), rho(i,j), (ff(i,j,k) , k=0,npop-1), (ffeq(i,j,k) , k=0,npop-1)
       write(90,*) T(i,j), (gg(i,j,k) , k=0,npop-1), (ggeq(i,j,k) , k=0,npop-1)
    end do
 end do

 close(90)

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine dump_out"
#endif

 end subroutine dump_out
