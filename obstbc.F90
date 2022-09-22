!--------------------------------------------------------
        subroutine obstbc
!--------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        k = nx / 4

        do j = ny/2-nobst/2+1,ny/2+nobst/2
           f1(k+1,j) = f3(k+1,j)
           f3(k  ,j) = f1(k  ,j)
        enddo

        do j = ny/2-nobst/2,ny/2+nobst/2+1
           f5(k+1,j) = f7(k+1,j)
           f8(k+1,j) = f6(k+1,j)
           f7(k,  j) = f5(k,  j)
           f6(k,  j) = f8(k,  j)
        enddo

        return
        end subroutine obstbc
