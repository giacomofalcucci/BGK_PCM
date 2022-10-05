! ================================
  subroutine prof_x
! ================================

  use shared

  write(667,1005) it
  do i = 1, Nx
     write (667,1002) i,           & 
                      u(i,Ny/2),   & 
                      v(i,Ny/2),   & 
                      rho(i,Ny/2), & 
                      T(i,Ny/2),   & 
                      phi(i,Ny/2)
  end do
  write(667,'(a1)')
  write(667,'(a1)')


! formats...
1002    format(i8,5(e14.6,1x))
1005    format("# t=",i7)

 
#ifdef DEBUG
     write(6,*) "DEBUG: Completed subroutine prof_x"
#endif

  end subroutine prof_x
