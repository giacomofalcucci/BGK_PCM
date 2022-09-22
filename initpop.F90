! --------------------------------------------------
        subroutine initpop
!---------------------------------------------------

! ------- modules
        use storage
        implicit double precision(a-h,o-z)
           
        if (dump.eq.1)then
           call resume
        else
           do j = 1, ny
              do i = 1, nx
!
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
!
                 f0(i,j) = cte04*rhoij*(1.0d0 - sumsq)

                 f1(i,j) = cte09*rhoij*(1.0d0 - sumsq + u22 + ui)
                 f2(i,j) = cte09*rhoij*(1.0d0 - sumsq + v22 + vi)
                 f3(i,j) = cte09*rhoij*(1.0d0 - sumsq + u22 - ui)
                 f4(i,j) = cte09*rhoij*(1.0d0 - sumsq + v22 - vi)

                 f5(i,j) = cte36*rhoij*(1.0d0 + sumsq2 +ui+vi+uv)
                 f6(i,j) = cte36*rhoij*(1.0d0 + sumsq2 -ui+vi-uv)
                 f7(i,j) = cte36*rhoij*(1.0d0 + sumsq2 -ui-vi+uv)
                 f8(i,j) = cte36*rhoij*(1.0d0 + sumsq2 +ui-vi-uv)


                 g0(i,j) = cte04*T(i,j)*(1.0d0 - sumsq)

                 g1(i,j) = cte09*T(i,j)*(1.0d0 - sumsq + u22 + ui)
                 g2(i,j) = cte09*T(i,j)*(1.0d0 - sumsq + v22 + vi)
                 g3(i,j) = cte09*T(i,j)*(1.0d0 - sumsq + u22 - ui)
                 g4(i,j) = cte09*T(i,j)*(1.0d0 - sumsq + v22 - vi)

                 g5(i,j) = cte36*T(i,j)*(1.0d0 + sumsq2 +ui+vi+uv)
                 g6(i,j) = cte36*T(i,j)*(1.0d0 + sumsq2 -ui+vi-uv)
                 g7(i,j) = cte36*T(i,j)*(1.0d0 + sumsq2 -ui-vi+uv)
                 g8(i,j) = cte36*T(i,j)*(1.0d0 + sumsq2 +ui-vi-uv)

              enddo
           enddo
        endif

        rhoaver = 0.d0
        do j = 1, ny
           do i = 1, nx
               rhoaver = rhoaver + f0(i,j) +  & 
     &                   f1(i,j)+f2(i,j)+f3(i,j)+f4(i,j)+ &
     &                   f5(i,j)+f6(i,j)+f7(i,j)+f8(i,j)    
           enddo
        enddo

        rhoaver = rhoaver / dfloat(nx*ny)
        dinvrho = 1.d0 / rhoaver
        print*,'average density at t=0:', rhoaver
        write(*,*) rhod1(25,25), rhod1(32,32), rhod1(60,60)
        write(*,*) 'rhoaver =', rhoaver
        write(*,*) 'pop -->', f5(32,32)
        return

#ifdef DEBUG
        write(6,*) "Completed subroutine initpop"
#endif


        end subroutine initpop
