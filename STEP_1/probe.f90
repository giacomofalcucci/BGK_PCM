! ================================
  subroutine probe
! ================================

  use shared

     write (666,*) it,             u(Nx/2,Ny/2),   v(Nz/2,Ny/2), & 
                    rho(Nx/2,Ny/2), T(Nx/2,Ny/2), phi(Nz/2,Ny/2)

! formats...
1001    format(i8,5(e14.6,1x))

 
#ifdef DEBUG
     write(6,*) "DEBUG: Completed subroutine probe"
#endif

  end subroutine probe
