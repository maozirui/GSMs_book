    program SPH3D
    use config_parameter
    implicit none

    real(8) ntotal, itype(maxn), maxtimestep, d, m, i, yesorno
    real(8) x(dim,maxn), vx(dim,maxn), mass(maxn), rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn), hsml(maxn), dt,dx
    real(8) s1,s2,rho0

    call time_elapsed(s1)

    if (taylor)  then
        dt = 1.e-9
    else if (hvisphere) then
        dt = .5e-9
    else if (tungstencube) then
        dt = 1.e-4
    else if (TNT) then
        dt =5e-8
    else if (landslide)then
        dt =2.0e-6
        !
    end if

    call input(x,vx,mass,rho,p,u,itype,hsml,ntotal,dx,rho0)
1   write(*,*)'  ***************************************************'
    write(*,*)'          Please input the maximal time period '
    write(*,*)'  ***************************************************'
    read(*,*) maxtimestep
    !maxtimestep = 50000
    call time_integration(x, dx,vx, mass, rho, p, itype, hsml, ntotal, maxtimestep, dt,rho0)
    write(*,*)'  ***************************************************'
    write(*,*) 'Are you going to run more time steps ? (0=No, 1=yes)'
    write(*,*)'  ***************************************************'
    read (*,*) yesorno
    !yesorno = 0
    if(yesorno.ne.0) then
        go to 1
        !   call time_print
        call time_elapsed(s2)
        write (*,*)'        Elapsed CPU time = ', s2-s1
    end if

    pause
    end program