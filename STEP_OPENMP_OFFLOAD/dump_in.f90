 subroutine dump_in

 use shared

 open(unit=90,file='dumpFile.dat')

 read(90,*) it

 do j = 0,Ny+1
    do i = 0,Nx+1
       read(90,*) u(i,j), v(i,j), rho(i,j), & 
                  f0(i,j), f1(i,j), f2(i,j), f3(i,j), f4(i,j), & 
                  f5(i,j), f6(i,j), f7(i,j), f8(i,j)
       read(90,*) T(i,j), &  
                  g0(i,j), g1(i,j), g2(i,j), g3(i,j), g4(i,j), & 
                  g5(i,j), g6(i,j), g7(i,j), g8(i,j)
    end do
 end do

 close(90)

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine dump_in"
#endif

 end subroutine dump_in
