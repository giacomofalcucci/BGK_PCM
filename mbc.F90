!-------------------------------------------------------------
        subroutine mbc
!-------------------------------------------------------------

! ------- modules 
        use storage
        implicit double precision(a-h,o-z)
! EAST case
        do j = 1,ny
           f1(0,j) = f1(nx,j)
           f5(0,j) = f5(nx,j)
           f8(0,j) = f8(nx,j)
        enddo

! WEST case
        do j = 1,ny
           f3(nx+1,j) = f3(1,j)
           f6(nx+1,j) = f6(1,j)
           f7(nx+1,j) = f7(1,j)
        enddo

! NORTH case
        do i = 1,nx
           f4(i,ny+1) = f2(i,ny)
           f8(i,ny+1) = f6(i,ny)
           f7(i,ny+1) = f5(i,ny)
        enddo

! SOUTH case
        do i = 1,nx
           f2(i,0) = f4(i,1)
           f6(i,0) = f8(i,1)
           f5(i,0) = f7(i,1)
        enddo

        return
        end subroutine mbc
