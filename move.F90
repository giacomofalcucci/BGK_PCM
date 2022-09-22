! --------------------------------------------------
        subroutine move
!---------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)

!$acc kernels
!$acc loop independent collapse(2)
        do j = 1,ny
           do i = 1, nx
!              fp0(i,j) = f0(i,j)
              fp1(i,j) = f1(i-1,j)
              fp2(i,j) = f2(i,j-1)
              fp3(i,j) = f3(i+1,j)
              fp4(i,j) = f4(i,j+1)
              fp5(i,j) = f5(i-1,j-1)
              fp6(i,j) = f6(i+1,j-1)
              fp7(i,j) = f7(i+1,j+1)
              fp8(i,j) = f8(i-1,j+1)
!
              gp1(i,j) = g1(i-1,j)
              gp2(i,j) = g2(i,j-1)
              gp3(i,j) = g3(i+1,j)
              gp4(i,j) = g4(i,j+1)
              gp5(i,j) = g5(i-1,j-1)
              gp6(i,j) = g6(i+1,j-1)
              gp7(i,j) = g7(i+1,j+1)
              gp8(i,j) = g8(i-1,j+1)
           enddo
        enddo
!$acc end kernels

!!!###GF##!!! border fix
!!!###GF##!!
!!!###GF##!!        j = 0
!!!###GF##!!!$acc kernels
!!!###GF##!!!$acc loop independent 
!!!###GF##!!        do i = 0, nx+1 
!!!###GF##!!!           fp0(i,j)=f0(i,j)
!!!###GF##!!           fp1(i,j)=f1(i,j)
!!!###GF##!!           fp2(i,j)=f2(i,j)
!!!###GF##!!           fp3(i,j)=f3(i,j)
!!!###GF##!!           fp4(i,j)=f4(i,j)
!!!###GF##!!           fp5(i,j)=f5(i,j)
!!!###GF##!!           fp6(i,j)=f6(i,j)
!!!###GF##!!           fp7(i,j)=f7(i,j)
!!!###GF##!!           fp8(i,j)=f8(i,j)
!!!###GF##!!!
!!!###GF##!!           gp1(i,j)=g1(i,j)
!!!###GF##!!           gp2(i,j)=g2(i,j)
!!!###GF##!!           gp3(i,j)=g3(i,j)
!!!###GF##!!           gp4(i,j)=g4(i,j)
!!!###GF##!!           gp5(i,j)=g5(i,j)
!!!###GF##!!           gp6(i,j)=g6(i,j)
!!!###GF##!!           gp7(i,j)=g7(i,j)
!!!###GF##!!           gp8(i,j)=g8(i,j)
!!!###GF##!!        enddo
!!!###GF##!!!$acc end kernels
!!!###GF##!!!
!!!###GF##!!        j = ny+1
!!!###GF##!!!$acc kernels
!!!###GF##!!!$acc loop independent 
!!!###GF##!!        do i = 0, nx+1 
!!!###GF##!!!           fp0(i,j)=f0(i,j)
!!!###GF##!!           fp1(i,j)=f1(i,j)
!!!###GF##!!           fp2(i,j)=f2(i,j)
!!!###GF##!!           fp3(i,j)=f3(i,j)
!!!###GF##!!           fp4(i,j)=f4(i,j)
!!!###GF##!!           fp5(i,j)=f5(i,j)
!!!###GF##!!           fp6(i,j)=f6(i,j)
!!!###GF##!!           fp7(i,j)=f7(i,j)
!!!###GF##!!           fp8(i,j)=f8(i,j)
!!!###GF##!!!
!!!###GF##!!           gp1(i,j)=g1(i,j)
!!!###GF##!!           gp2(i,j)=g2(i,j)
!!!###GF##!!           gp3(i,j)=g3(i,j)
!!!###GF##!!           gp4(i,j)=g4(i,j)
!!!###GF##!!           gp5(i,j)=g5(i,j)
!!!###GF##!!           gp6(i,j)=g6(i,j)
!!!###GF##!!           gp7(i,j)=g7(i,j)
!!!###GF##!!           gp8(i,j)=g8(i,j)
!!!###GF##!!        enddo
!!!###GF##!!!$acc end kernels
!!!###GF##!!!
!!!###GF##!!        i = 0
!!!###GF##!!!$acc kernels
!!!###GF##!!!$acc loop independent 
!!!###GF##!!        do j = 0, ny+1
!!!###GF##!!!           fp0(i,j)=f0(i,j)
!!!###GF##!!           fp1(i,j)=f1(i,j)
!!!###GF##!!           fp2(i,j)=f2(i,j)
!!!###GF##!!           fp3(i,j)=f3(i,j)
!!!###GF##!!           fp4(i,j)=f4(i,j)
!!!###GF##!!           fp5(i,j)=f5(i,j)
!!!###GF##!!           fp6(i,j)=f6(i,j)
!!!###GF##!!           fp7(i,j)=f7(i,j)
!!!###GF##!!           fp8(i,j)=f8(i,j)
!!!###GF##!!!
!!!###GF##!!           gp1(i,j)=g1(i,j)
!!!###GF##!!           gp2(i,j)=g2(i,j)
!!!###GF##!!           gp3(i,j)=g3(i,j)
!!!###GF##!!           gp4(i,j)=g4(i,j)
!!!###GF##!!           gp5(i,j)=g5(i,j)
!!!###GF##!!           gp6(i,j)=g6(i,j)
!!!###GF##!!           gp7(i,j)=g7(i,j)
!!!###GF##!!           gp8(i,j)=g8(i,j)
!!!###GF##!!        enddo
!!!###GF##!!!$acc end kernels
!!!###GF##!!!
!!!###GF##!!        i = nx+1
!!!###GF##!!!$acc kernels
!!!###GF##!!!$acc loop independent 
!!!###GF##!!        do j = 0, ny+1
!!!###GF##!!!           fp0(i,j)=f0(i,j)
!!!###GF##!!           fp1(i,j)=f1(i,j)
!!!###GF##!!           fp2(i,j)=f2(i,j)
!!!###GF##!!           fp3(i,j)=f3(i,j)
!!!###GF##!!           fp4(i,j)=f4(i,j)
!!!###GF##!!           fp5(i,j)=f5(i,j)
!!!###GF##!!           fp6(i,j)=f6(i,j)
!!!###GF##!!           fp7(i,j)=f7(i,j)
!!!###GF##!!           fp8(i,j)=f8(i,j)
!!!###GF##!!!
!!!###GF##!!           gp1(i,j)=g1(i,j)
!!!###GF##!!           gp2(i,j)=g2(i,j)
!!!###GF##!!           gp3(i,j)=g3(i,j)
!!!###GF##!!           gp4(i,j)=g4(i,j)
!!!###GF##!!           gp5(i,j)=g5(i,j)
!!!###GF##!!           gp6(i,j)=g6(i,j)
!!!###GF##!!           gp7(i,j)=g7(i,j)
!!!###GF##!!           gp8(i,j)=g8(i,j)
!!!###GF##!!        enddo
!!!###GF##!!!$acc end kernels
!!!###GF##!!!
#ifdef DEBUG
        write(6,*) "Completed subroutine move"
#endif
        end subroutine move
