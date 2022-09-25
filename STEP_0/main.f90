! ================================
 program main
! ================================

 use shared
 implicit none

 integer:: icountT0, icountT1, icount_rate, icount_max
 integer:: icountL0, icountL1
 real(kind=8)    :: time_inn_loop

 call setup

! init section
 if (dump .eq. 0) then
    it  = 0
    it0 = 1
    call initialization
    open(unit=150,file='log.dat')
    write(150,*) 'iteration','    ','mass'
    write(150,*) '===================='
    close(150)
 else
    call dump_in
    it0 = it+1
 end if

 deltaOut = int(iter/noutput)
 if (it .eq. 0) then
    itOut = deltaOut
    call out2d
 else
    itOut = it + deltaOut
 end if

 open(unit=164,file='non-dim_bulk_numbers.dat',status='unknown')
 open(unit=165,file='Temp_trend.dat',status='unknown')

! time loop section
 call SYSTEM_CLOCK(icountL0, icount_rate, icount_max)
 icountT0=icountL0
 do it = it0,iter

#ifdef DEBUG
    write(6,*) "DEBUG: starting iteration ",it
#endif

    call BC_f
    call BC_g
    call streaming
    call moments
    call moments2
!
    call phase_field
    call forces
!
    call regular
!
    call equilibrium
    call collision
    call moments
    call media
    call mass_check


    if(mod(it,1000).eq.0) then
      call diag    
    endif

    open(unit=150,file='log.dat',position='append')
    write(150,"(I8,A,F10.4)") it,' ',mass
    close(150)

!!!!    do i=1,Nx
!!!!      if(phi(i,int(Ny/2)).eq.0) then
!!!!        write(128,*) it, i, phi(i,int(Ny/2))
!!!!      endif
!!!!    enddo

    if (it == it0) then
       write(*,*) it
    end if 
 
    if (it == itOut .or. it == iter) then

       call out2d

       write(*,"(I8,4F10.4)") it, mass, Nu(samplex), Nu2(samplex), Nu3(samplex)
       write(*,"(F10.4)") gradTwall2(samplex), Tavg(samplex)

       open(unit=160,file='Nusselt.dat')
       do i = 1,samplex
          write(160,"(I8,A,F10.4,A,F10.4,A,F10.4)") i,' ',Nu(i),' ',Nu2(i),' ',Nu3(i)
       end do
       close(160)

       write(*,*) 'saving results...'
       !call dump_out
       itOut = itOut + deltaOut
       !write(*,*) 'dumping completed'
    end if

    if(mod(it,10000).eq.0) then
       call dump_out
       !itOut = itOut + deltaOut
       write(*,*) 'dumping completed'
    endif

! Here I'm signal....
    if (mod(it,100).eq.0) then
!        write(6,*) "timestep ", it, "done"
       call SYSTEM_CLOCK(icountL1, icount_rate, icount_max)
       time_inn_loop = real(icountL1-icountL0)/(icount_rate)
       write(6,1003)(time_inn_loop)/100,it,iter
       call SYSTEM_CLOCK(icountL0, icount_rate, icount_max)
    endif

 end do

 call SYSTEM_CLOCK(icountT1, icount_rate, icount_max)
 time_inn_loop = real(icountT1-icountT0)/(icount_rate)
 write(6,*) "=========================================="
 write(6,*) "INFO: Run completed"
 write(6,1004)(time_inn_loop)
 write(6,1005)(time_inn_loop)/iter
 write(6,1006) (float(nx)/1000.0*float(ny)/1000.0)             &
               *(float(iter)/time_inn_loop)
 write(6,*) "=========================================="

! formats...
1001    format(" INFO: Init   time           ",1(e14.6,1x))
1003    format(" INFO: Mean   time           ",1(e14.6,1x),i8,"/",i8)
1004    format(" INFO: Main loop time        ",1(e14.6,1x))
1005    format(" INFO: Mean time for timestep",1(e14.6,1x))
1006    format(" INFO: Mlups                 ",1(f14.6,1x))

 end program main
