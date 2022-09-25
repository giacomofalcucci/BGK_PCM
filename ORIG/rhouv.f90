 subroutine rhouv

 use shared
 use omp_lib

 call omp_set_num_threads(n_th)


!$OMP PARALLEL DEFAULT(shared) &
!$OMP PRIVATE(i,j)
!$OMP DO
 do j = 1,Ny
    do i = 1,Nx
       rho(i,j) = f(0,i,j)+f(1,i,j)+f(2,i,j)+f(3,i,j)+f(4,i,j)+f(5,i,j)+f(6,i,j)+f(7,i,j)+f(8,i,j)
       u(i,j)   = (1.d0/rho(i,j))*(f(1,i,j)-f(3,i,j)+f(6,i,j)-f(6,i,j)-f(7,i,j)+f(8,i,j))
       v(i,j)   = (1.d0/rho(i,j))*(f(2,i,j)-f(4,i,j)+f(5,i,j)+f(6,i,j)-f(7,i,j)-f(8,i,j))
    end do
 end do
!$OMP END DO
!$OMP END PARALLEL

 end subroutine rhouv
