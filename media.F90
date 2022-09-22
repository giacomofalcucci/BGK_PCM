!----------------------------------------------------------
        subroutine media
!--------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)

        do i=1,nx
           do j=1,ny
              u1(i,j)=(u1(i,j)+u2(i,j))*0.5d0
              v1(i,j)=(v1(i,j)+v2(i,j))*0.5d0
           end do
        end do

        return
        end subroutine media
