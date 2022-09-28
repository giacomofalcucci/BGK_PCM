! ================================
  subroutine prof_y
! ================================

  use shared

  write(668,1005) it
  do j = 1, Ny
     write (668,1002) j,           & 
                      u(Nx/2,j),   & 
                      v(Nx/2,j),   & 
                      jho(Nx/2,j), & 
                      T(Nx/2,j),   & 
                      phi(Nx/2,j)
  end do
  write(668,'(a1)')
  write(668,'(a1)')


! formats...
1002    format(i8,5(e14.6,1x))
1005    format("# t=",i7)

 
#ifdef DEBUG
     write(6,*) "DEBUG: Completed subroutine prof_y"
#endif

  end subroutine prof_y
