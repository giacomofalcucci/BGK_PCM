!-----------------------------------------------------------
        subroutine pbc
!-----------------------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)


!$acc kernels
!$acc loop independent
        do j = 1, ny
            f0(0   ,j) = f0(nx,j)
            f0(nx+1,j) = f0( 1,j)
            f1(0   ,j) = f1(nx,j)
            f1(nx+1,j) = f1( 1,j)
            f2(0   ,j) = f2(nx,j)
            f2(nx+1,j) = f2( 1,j)
            f3(0   ,j) = f3(nx,j)
            f3(nx+1,j) = f3( 1,j)
            f4(0   ,j) = f4(nx,j)
            f4(nx+1,j) = f4( 1,j)
            f5(0   ,j) = f5(nx,j)
            f5(nx+1,j) = f5( 1,j)
            f6(0   ,j) = f6(nx,j)
            f6(nx+1,j) = f6( 1,j)
            f7(0   ,j) = f7(nx,j)
            f7(nx+1,j) = f7( 1,j)
            f8(0   ,j) = f8(nx,j)
            f8(nx+1,j) = f8( 1,j)
        enddo
!$acc end kernels

!$acc kernels
!$acc loop independent
        do i = 1, nx
            f0(i,0   ) = f0(i,ny)
            f0(i,ny+1) = f0(i, 1)
            f1(i,0   ) = f1(i,ny)
            f1(i,ny+1) = f1(i, 1)
            f2(i,0   ) = f2(i,ny)
            f2(i,ny+1) = f2(i, 1)
            f3(i,0   ) = f3(i,ny)
            f3(i,ny+1) = f3(i, 1)
            f4(i,0   ) = f4(i,ny)
            f4(i,ny+1) = f4(i, 1)
            f5(i,0   ) = f5(i,ny)
            f5(i,ny+1) = f5(i, 1)
            f6(i,0   ) = f6(i,ny)
            f6(i,ny+1) = f6(i, 1)
            f7(i,0   ) = f7(i,ny)
            f7(i,ny+1) = f7(i, 1)
            f8(i,0   ) = f8(i,ny)
            f8(i,ny+1) = f8(i, 1)
        enddo
!$acc end kernels

        f0(0   ,   0) = f0(nx  ,   0)
        f0(nx+1,   0) = f0(nx+1,  ny)
        f0(nx+1,ny+1) = f0(   1,ny+1)
        f0(0   ,ny+1) = f0(   0,   1)

        f1(0   ,   0) = f1(nx  ,   0)
        f1(nx+1,   0) = f1(nx+1,  ny)
        f1(nx+1,ny+1) = f1(   1,ny+1)
        f1(0   ,ny+1) = f1(   0,   1)

        f2(0   ,   0) = f2(nx  ,   0)
        f2(nx+1,   0) = f2(nx+1,  ny)
        f2(nx+1,ny+1) = f2(   1,ny+1)
        f2(0   ,ny+1) = f2(   0,   1)

        f3(0   ,   0) = f3(nx  ,   0)
        f3(nx+1,   0) = f3(nx+1,  ny)
        f3(nx+1,ny+1) = f3(   1,ny+1)
        f3(0   ,ny+1) = f3(   0,   1)

        f4(0   ,   0) = f4(nx  ,   0)
        f4(nx+1,   0) = f4(nx+1,  ny)
        f4(nx+1,ny+1) = f4(   1,ny+1)
        f4(0   ,ny+1) = f4(   0,   1)

        f5(0   ,   0) = f5(nx  ,   0)
        f5(nx+1,   0) = f5(nx+1,  ny)
        f5(nx+1,ny+1) = f5(   1,ny+1)
        f5(0   ,ny+1) = f5(   0,   1)

        f6(0   ,   0) = f6(nx  ,   0)
        f6(nx+1,   0) = f6(nx+1,  ny)
        f6(nx+1,ny+1) = f6(   1,ny+1)
        f6(0   ,ny+1) = f6(   0,   1)

        f7(0   ,   0) = f7(nx  ,   0)
        f7(nx+1,   0) = f7(nx+1,  ny)
        f7(nx+1,ny+1) = f7(   1,ny+1)
        f7(0   ,ny+1) = f7(   0,   1)

        f8(0   ,   0) = f8(nx  ,   0)
        f8(nx+1,   0) = f8(nx+1,  ny)
        f8(nx+1,ny+1) = f8(   1,ny+1)
        f8(0   ,ny+1) = f8(   0,   1)

        
#ifdef DEBUG
        write(6,*) "Completed subroutine pbc"
#endif

        end subroutine pbc
