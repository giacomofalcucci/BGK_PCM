 subroutine media

 use shared

!GA1 !$OMP PARALLEL DEFAULT(NONE)  &
!GA1 !$OMP PRIVATE(i,j)            &
!GA1 !$OMP SHARED(Nx,Ny)           &
!GA1 !$OMP SHARED(u,v,T)  &
!GA1 !$OMP SHARED(u2,v2,T2)  
!GA1 !$OMP DO
!$OMP target teams distribute parallel do collapse(2)

 do j = 1,Ny
    do i = 1,Nx
             T(i,j)     = (T(i,j)+T2(i,j))*0.50d0  ! THIS, I have to check....
   
             u(i,j)     = (u(i,j)+u2(i,j))*0.50d0
             v(i,j)     = (v(i,j)+v2(i,j))*0.50d0
    end do
 end do
!GA1  !$OMP END PARALLEL


#ifdef DEBUG
        write(6,*) "DEBUG: Completed subroutine media"
#endif

 end subroutine media
