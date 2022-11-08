  subroutine out2d_Temp

  use shared
  implicit none

  character*25 file_nameP
  

      file_nameP = 'OutTem_xxxxxxxx.vtk'
      write(file_nameP(8:15),4000) it
      open(524,file=file_nameP,status='unknown')

! ASCII header
      write(524,'(A)') '# vtk DataFile Version 2.0'
      write(524,'(A)') 'Campo'
      write(524,'(A)') 'BINARY'
      write(524,'(A)') 'DATASET STRUCTURED_POINTS'
      write(524,'(A,3I8)') 'DIMENSIONS', Nx, Ny, 1
      write(524,'(A7,I10,A1,I10,A1,I10)')  'ORIGIN ',1,' ',1,' ',1
      write(524,'(A8,I10,A1,I10,A1,I10)') 'SPACING ',1,' ',1,' ',1
      write(524,'(A10,I10)')'POINT_DATA ',Nx*Ny
#ifdef SINGLEPRECISION
      write(524,'(A)')'SCALARS Temperature float'
#else
      write(524,'(A)')'SCALARS Temperature double'
#endif
      write(524,'(A)')'LOOKUP_TABLE default'
      close(524)

! BINARY section 
      open(524,file=file_nameP,status='old', position='append', &
           form='unformatted',access='STREAM',CONVERT="BIG_ENDIAN")

      do j = 1, Ny
         do i = 1, Nx
            write (524) T(i,j)
         end do
      enddo

     close(524)

! format
4000  format(i8.8)
 
!#ifdef DEBUG
     write(6,*) "DEBUG: Completed subroutine out2d_Temp"
!#endif

  end subroutine out2d_Temp
