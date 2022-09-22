!----------------------------------------------------------
        subroutine calcsurf
!----------------------------------------------------------
! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        do j=1,ny
           do i=1,nx
              iflag(i,j)=0
              if (rhod1(i,j).gt.0.865d0) then 
                 iflag(i,j)=1
              endif
           enddo
        enddo

        isurf=0

        do j=1,ny
           do i=1,nx
              if(iflag(i,j).ne.iflag(i,j+1).or. &
                 iflag(i,j).ne.iflag(i,j-1).or. &
                 iflag(i,j).ne.iflag(i+1,j).or. &
                 iflag(i,j).ne.iflag(i-1,j)) then
                 isurf=isurf+1
              endif
           enddo
        enddo    

        isurf=isurf/2.d0

        write(88,*)istep,isurf

        return
        end subroutine calcsurf

