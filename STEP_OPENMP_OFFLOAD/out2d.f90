  subroutine out2d

  use shared
  character*25 file_name2, file_name3

4000  format(i8.8)

      file_name2 = 'Output_xxxxxxxx.vtk'
      write(file_name2(8:15),4000) it
      open(523,file=file_name2,status='unknown')

      write (523, '(A)') '# vtk DataFile Version 2.0'
      write (523, '(A)') 'Griglia strutturata'
      write (523, '(A)') 'ASCII'
      write (523, '(A)') 'DATASET STRUCTURED_GRID'
      write (523, '(A,3I8)') 'DIMENSIONS', Nx+2, Ny+2, 1
      write (523, '(A,I15,A)') 'POINTS', (Nx+2) * (Ny+2) * 1, ' float'

      do l = 0, Ny+1
         do i = 0, Nx+1
            write (523,*) i, l, 1
         end do
      end do
!
      write (523, '(A,I15)') 'POINT_DATA', (Nx+2) * (Ny+2) * 1
      write (523, '(A)') 'SCALARS DENSITY float 1'
      write (523, '(A)') 'LOOKUP_TABLE default'
!
!GA     !$omp target exit data map(from:rho)
     do l = 0, Ny+1
         do i = 0, Nx+1
            write (523,*) rho(i, l)
         end do
     enddo

!      write (523, '(A,I15)') 'POINT_DATA', (Nx+2) * (Ny+2) * 1
      write (523, '(A)') 'SCALARS Temperature float 1'
      write (523, '(A)') 'LOOKUP_TABLE default'
!
!GA     !$omp target exit data map(from:T)
     do l = 0, Ny+1
         do i = 0, Nx+1
            write (523,*) T(i, l)
         end do
     enddo
!!!
!!!!      write (523, '(A,I15)') 'POINT_DATA', (Nx) * (Ny) * 1
      write (523, '(A)') 'SCALARS Phase_Field float 1'
      write (523, '(A)') 'LOOKUP_TABLE default'

!GA     !$omp target exit data map(from:phi)
     do l = 0, Ny+1
         do i = 0, Nx+1
            write (523,*) phi(i, l)
         end do
     enddo


      write (523, '(A)') 'VECTORS VELOCITY float'
!
!GA     !$omp target exit data map(from:u,v)
     do l = 0, Ny+1
         do i = 0, Nx+1
            write (523,*) u(i, l), v(i, l), 0.0d0
         end do
     enddo


     close(523)
 
#ifdef DEBUG
     write(6,*) "DEBUG: Completed subroutine out2d"
#endif

  end subroutine out2d
