 program main

 use shared

 call setup

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
!!    call output_paraview_structured
 else
    itOut = it + deltaOut
 end if

 open(unit=164,file='non-dim_bulk_numbers.dat',status='unknown')
 open(unit=165,file='Temp_trend.dat',status='unknown')

 do it = it0,iter

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


    if(mod(it,200).eq.0) then
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

!!       call output_paraview_structured
 
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

 end do

 end program main
