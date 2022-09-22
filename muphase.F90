! ================================
        program bgk2D
! ================================
! ------- modules
        use storage
        implicit double precision(a-h,o-z)

! ------- input parameters

        omega=1.d0 

        call SYSTEM_CLOCK(icountI0, icount_rate, icount_max)
        call input
        call alloca

! -------initialisation

        call inithydro1
        call initpop
        call SYSTEM_CLOCK(icountI1, icount_rate, icount_max)
        time_init = real(icountI1-icountI0)/(icount_rate)
        write(6,1001) time_init

! ------- MAIN LOOP

        call SYSTEM_CLOCK(icountL0, icount_rate, icount_max)
        icountT0=icountL0
        do 10 istep = 1,nsteps
 
          if(icond.eq.4) then
            call BC_f
            call BC_g
          else
           call pbcdens
           call pbc
          endif
           call move
           call hydrovarPRE
           call collis
!           call hydrovarPOST
!
           if (mod(istep,ndiag).eq.0) then
              call out0d
              call energy
              call calcsurf
           endif
!
           if (mod(istep,nout).eq.0) then
              call out1d(frce)
              call out2d(frce)
           endif
!
           if (mod(istep,100).eq.0) then
              call SYSTEM_CLOCK(icountL1, icount_rate, icount_max)
              time_inn_loop = real(icountL1-icountL0)/(icount_rate)
              write(6,1003)(time_inn_loop)/100,istep,nsteps
              call SYSTEM_CLOCK(icountL0, icount_rate, icount_max)
           endif
!
           if (mod(istep,200000).eq.0) then
              call savepop
           endif
!
           if(icond == 1) then 
              if (mod(istep,nsteps).eq.0) then
                 call laplace    !  it needs onle for Case 1
              endif
           endif
!
!-------- end of main loop
            
10      continue

        call SYSTEM_CLOCK(icountT1, icount_rate, icount_max)
        time_inn_loop = real(icountT1-icountT0)/(icount_rate)
        write(6,*) "=========================================="
        write(6,*) "INFO: Run completed" 
        write(6,1004)(time_inn_loop)
        write(6,1005)(time_inn_loop)/nsteps
        write(6,1006) (float(nx)/1000.0*float(ny)/1000.0)             &
                     *(float(nsteps)/time_inn_loop)
        write(6,*) "=========================================="

! formats...
1001    format(" INFO: Init   time           ",1(e14.6,1x))
1003    format(" INFO: Mean   time           ",1(e14.6,1x),i8,"/",i8)
1004    format(" INFO: Main loop time        ",1(e14.6,1x))
1005    format(" INFO: Mean time for timestep",1(e14.6,1x))
1006    format(" INFO: Mlups                 ",1(f14.6,1x))

        end program bgk2D

