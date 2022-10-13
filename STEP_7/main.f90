! ================================
 program main
! ================================
!
 use shared
 use check_mem
! 
 implicit none

 integer:: icountT0, icountT1, icount_rate, icount_max
 integer:: icountL0, icountL1
 integer:: isignal
 real(kind=dp)    :: time_inn_loop
 real(kind=dp)    :: mem_start, mem_stop

 write(6,*) "===================================="
 write(6,*) " Serial Optimized version (Step_7) "
 write(6,*) "===================================="

 call setup
!
 isignal=500

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
 mem_start = get_mem();
 write(6,1007) mem_start
 call SYSTEM_CLOCK(icountL0, icount_rate, icount_max)
 icountT0=icountL0
 do it = it0,iter

#ifdef DEBUG
    write(6,*) "DEBUG: starting iteration ",it
#endif

    call BC_f
    call BC_g
    call streaming
    call momentsPre
!
    call phase_field
    call forces
!
    call collision
    call moments


    if(mod(it,1000).eq.0) then
      call diag    
      call probe
      call prof_x
      call prof_y
      call field_check

      open(unit=150,file='log.dat',position='append')
      write(150,"(I8,A,F10.4)") it,' ',mass
      close(150)
    endif
  
    if (it == it0) then
       write(*,*) "starting at ", it
    end if 
 
    if (it == itOut .or. it == iter) then

!      vtk dump
       call out2d
!
!GA      write(6,1010) it, mass 
!GA      write(6,1011) Nu(samplex), Nu2(samplex), Nu3(samplex)
!GA      write(6,1012) gradTwall2(samplex), Tavg(samplex)
!GA
       open(unit=160,file='Nusselt.dat')
       do i = 1,samplex
          write(160,"(I8,A,F10.4,A,F10.4,A,F10.4)")  & 
                                       i,' ',Nu(i),' ',Nu2(i),' ',Nu3(i)
       end do
       close(160)

       itOut = itOut + deltaOut
    end if

!GA    if(mod(it,10000).eq.0) then
!GA       call dump_out
!GA       !itOut = itOut + deltaOut
!GA       write(*,*) 'dumping completed'
!GA       Dump to fix/check
!GA    endif

! Here I'm signal....
    if (mod(it,isignal).eq.0) then
       call SYSTEM_CLOCK(icountL1, icount_rate, icount_max)
       time_inn_loop = real(icountL1-icountL0)/(icount_rate)
       write(6,1003)(time_inn_loop)/isignal,it,iter
       call SYSTEM_CLOCK(icountL0, icount_rate, icount_max)
    endif

 end do

 call SYSTEM_CLOCK(icountT1, icount_rate, icount_max)
 time_inn_loop = real(icountT1-icountT0)/(icount_rate)
 mem_stop = get_mem();
 write(6,*) "=========================================="
 write(6,*) "INFO: Run completed"
 write(6,1004)(time_inn_loop)
 write(6,1005)(time_inn_loop)/iter
 write(6,1006) (float(nx)/1000.0*float(ny)/1000.0)             &
               *(float(iter)/time_inn_loop)
 write(6,1007) mem_stop
 write(6,*) "=========================================="

! formats...
1001    format("INFO: Init   time           ",1(e14.6,1x))
1003    format("INFO: Mean   time           ",1(e14.6,1x),i8,"/",i8)
1004    format("INFO: Main loop time        ",1(e14.6,1x))
1005    format("INFO: Mean time for timestep",1(e14.6,1x))
1006    format("INFO: Mlups                 ",1(f14.6,1x))
1007    format("INFO: Memory used (MB)      ",1(f14.6,1x))
1010    format("CHECK: it, mass             ",1(i8,f10.4))
1011    format("CHECK: Nu, Nu2, Nu3         ",3(f10.4,1x))
1012    format("CHECK: gradTwall2, Tavg     ",2(f10.4,1x))


 end program main
