!----------------------------------------------------------
        subroutine resume 
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        write(6,*)'################################'
        write(6,*)'reading  popolations ....'
        write(6,*)'################################'

        do j=0,ny+1
           do i=0,nx+1
              read(111)f0(i,j)
              read(111)f1(i,j)
              read(111)f2(i,j)
              read(111)f3(i,j)
              read(111)f4(i,j)
              read(111)f5(i,j)
              read(111)f6(i,j)
              read(111)f7(i,j)
              read(111)f8(i,j)
           enddo
        enddo

        do j=1,ny
           do i=1,nx
              read(113)u1(i,j),v1(i,j)
           enddo
        enddo

        do j=-1,ny+2
           do i=-1,nx+2
              read(112)rhod1(i,j)
           enddo
        enddo

#ifdef DEBUG
        write(6,*) "Completed subroutine resume"
#endif

        end subroutine resume

