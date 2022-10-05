 subroutine media

 use shared

 do concurrent(j=1:Ny, i=1:Nx)
             T(i,j)     = (T(i,j)+T2(i,j))*0.50d0  ! THIS, I have to check....
   
             u(i,j)     = (u(i,j)+u2(i,j))*0.50d0
             v(i,j)     = (v(i,j)+v2(i,j))*0.50d0
 end do

#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine media"
#endif

 end subroutine media
