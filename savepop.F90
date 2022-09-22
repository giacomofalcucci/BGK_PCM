!----------------------------------------------------------
        subroutine savepop
!----------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        write(6,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
        write(6,*)'savin populations  ....'
        write(6,*)'@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'  

        rewind(111)
        rewind(113)
        rewind(112)

        do j=0,ny+1
           do i=0,nx+1
              write(111)f0(i,j)      
              write(111)f1(i,j)      
              write(111)f2(i,j)      
              write(111)f3(i,j)      
              write(111)f4(i,j)      
              write(111)f5(i,j)      
              write(111)f6(i,j)      
              write(111)f7(i,j)      
              write(111)f8(i,j)      
           enddo
        enddo 
     
        do j=1,ny
           do i=1,nx
              write(113)u1(i,j),v1(i,j)
           enddo
        enddo

        do j=-1,ny+2
           do i=-1,nx+2
              write(112)rhod1(i,j)
           enddo
        enddo

        write(6,*)'.....done!!! :^)'
   
        return
        end subroutine savepop
       
