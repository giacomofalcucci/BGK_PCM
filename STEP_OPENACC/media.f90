 subroutine media

 use shared

!$acc kernels
 do j = 1,Ny
    do i = 1,Nx
             T(i,j)     = (T(i,j)+T2(i,j))*0.50d0  ! THIS, I have to check....
   
             u(i,j)     = (u(i,j)+u2(i,j))*0.50d0
             v(i,j)     = (v(i,j)+v2(i,j))*0.50d0
    end do
 end do
!$acc end kernels


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine media"
#endif

 end subroutine media
