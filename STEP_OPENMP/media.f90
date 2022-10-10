 subroutine media

 use shared

!$OMP PARALLEL DEFAULT(NONE)  &
!$OMP PRIVATE(i,j)            &
!$OMP SHARED(Nx,Ny)           &
!$OMP SHARED(u,v,T)  &
!$OMP SHARED(u2,v2,T2)  
!$OMP DO
 do j = 1,Ny
    do i = 1,Nx
             T(i,j)     = (T(i,j)+T2(i,j))*0.50d0  ! THIS, I have to check....
   
             u(i,j)     = (u(i,j)+u2(i,j))*0.50d0
             v(i,j)     = (v(i,j)+v2(i,j))*0.50d0
    end do
 end do
 !$OMP END PARALLEL


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine media"
#endif

 end subroutine media
