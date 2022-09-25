  subroutine output_paraview_structured

  use shared

  write(fn,*), it
  filename_1 = 'out_'//trim(adjustl(fn))//'.vts'
  filename = trim(adjustl(filename_1))

  write(Nxtot,*), Nx+2
  write(Nytot,*), Ny+2

  !========================== offset variables ==========================!
  ss = (Nx+2)*(Ny+2)*sizeof(sizereal)

  offset_velocity     = 0
  offset_density      = offset_velocity + sizeof(sizeint) + 2*ss
  offset_temp         = offset_density  + sizeof(sizeint) +   ss
  offset_phasef       = offset_temp     + sizeof(sizeint) +   ss
  offset_coordinates  = offset_phasef   + sizeof(sizeint) +   ss

  write(velocity,     *), offset_velocity
  write(density,      *), offset_density
  write(temperature,  *), offset_temp
  write(ph_field,     *), offset_phasef
  write(coordinates,  *), offset_coordinates
  !======================================================================!

  open(unit=96, file='part_1')
     string = '<?xml version="1.0"?>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
     write(96,"(A)"),trim(adjustl(string))
     string = '<StructuredGrid WholeExtent="1 '//trim(adjustl(Nxtot))//' 1 '//trim(adjustl(Nytot))//' 1 1">'
     write(96,"(A)"),trim(adjustl(string))
     string = '<Piece Extent="1 '//trim(adjustl(Nxtot))//' 1 '//trim(adjustl(Nytot))//' 1 1">'
     write(96,"(A)"),trim(adjustl(string))
     string = '<PointData>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<DataArray type="Float64" Name="velocity" NumberOfComponents="2" format="appended" &
	      offset="'//trim(adjustl(velocity))//'"/>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<DataArray type="Float64" Name="density" NumberOfComponents="1" format="appended" &
	      offset="'//trim(adjustl(density))//'"/>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<DataArray type="Float64" Name="temperature" NumberOfComponents="1" format="appended" &
              offset="'//trim(adjustl(temperature))//'"/>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<DataArray type="Float64" Name="ph_field" NumberOfComponents="1" format="appended" &
              offset="'//trim(adjustl(ph_field))//'"/>'
     write(96,"(A)"),trim(adjustl(string))
     string = '</PointData>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<Points>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<DataArray type="Float64" NumberOfComponents="3" format="appended" &
              offset="'//trim(adjustl(coordinates))//'"/>'
     write(96,"(A)"),trim(adjustl(string))
     string = '</Points>'
     write(96,"(A)"),trim(adjustl(string))
     string = '</Piece>'
     write(96,"(A)"),trim(adjustl(string))
     string ='</StructuredGrid>'
     write(96,"(A)"),trim(adjustl(string))
     string = '<AppendedData encoding="raw">'
     write(96,"(A)"),trim(adjustl(string))
  close(96)
  open(unit=96,file='part_2',form='unformatted',access='stream')
     write(96),'_',2*ss
        do j = 0,Ny+1
           do i = 0,Nx+1
              write(96) u(i,j), v(i,j)
           end do
        end do
     write(96),2*ss
        do j = 0,Ny+1
           do i = 0,Nx+1
              write(96) rho(i,j)
           end do
        end do
     write(96),ss
        do j = 0,Ny+1
           do i = 0,Nx+1
              write(96) T(i,j)
           end do
        end do
     write(96),ss
        do j = 0,Ny+1
           do i = 0,Nx+1
              write(96) phi(i,j)
           end do
        end do
     write(96),3*ss
        do j = 0,Ny+1
           do i = 0,Nx+1
              write(96) x(i),y(j),0.d0
           end do
        end do
  close(96)
  open(unit=96,file='part_3')
     write(96,*)
     string = '</AppendedData>'
     write(96,"(A)"),trim(adjustl(string))
     string = '</VTKFile>'
     write(96,"(A)"),trim(adjustl(string))
  close(96)

  command = 'cat '//trim(adjustl('part_1'))//' '//trim(adjustl('part_2'))//' '//trim(adjustl('part_3'))//' > '//filename
  call system(command)
  command = 'rm -f part_1 part_2 part_3'
  call system(command)
 
  end subroutine output_paraview_structured
