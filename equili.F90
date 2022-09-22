!-------------------------------------------------
        subroutine equili
!-------------------------------------------------
! ------- modules
        use storage
        implicit double precision(a-h,o-z)

        do j = 1, ny
           do i = 1, nx
              usq = u1(i,j) * u1(i,j) 
              vsq = v1(i,j) * v1(i,j)
              sumsq = (usq + vsq) / cs22
              sumsq2 = sumsq * (1.0d0 - cs2) / cs2
              u22 = usq / cssq 
              v22 = vsq / cssq
              ui = u1(i,j) / cs2
              vi = v1(i,j) / cs2
              uv = ui * vi
              rhoij = rhod1(i,j)

              feq(0,i,j) = (4.d0/9.d0)*rhoij*(1.0d0 - sumsq)

              feq(1,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + u22 + ui)
              feq(2,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + v22 + vi)
              feq(3,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + u22 - ui)
              feq(4,i,j) = (1.d0/9.d0)*rhoij*(1.0d0 - sumsq + v22 - vi)

              feq(5,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 +ui+vi+uv)
              feq(6,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 -ui+vi-uv)
              feq(7,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 -ui-vi+uv)
              feq(8,i,j) = (1.d0/36.d0)*rhoij*(1.0d0 + sumsq2 +ui-vi-uv)
           enddo
        enddo

        return
        end subroutine equili
