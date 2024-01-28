    subroutine input(x, vx, mass, rho, p, u, itype, hsml, ntotal,dx,rho0)
    use config_parameter
    !----------------------------------------------------------------------
    !     Subroutine for loading or generating initial particle information

    !     x-- coordinates of particles                                 [out]
    !     vx-- velocities of particles                                 [out]
    !     mass-- mass of particles                                     [out]
    !     rho-- dnesities of particles                                 [out]
    !     p-- pressure  of particles                                   [out]
    !     u-- internal energy of particles                             [out]
    !     itype-- types of particles                                   [out]
    !     hsml-- smoothing lengths of particles                        [out]
    !     ntotal-- total particle number                               [out]
    implicit none
    !
    real(8) itype(maxn), ntotal
    real(8) x(dim, maxn), vx(dim, maxn), mass(maxn),p(maxn), u(maxn), hsml(maxn), rho(maxn),x0,y0

    real(8) i, d, im,di,dj,jx,jy,j
    real(8) tempim,tempitype
    real(8) Lx,Ly,Nx,Ny,dx,rho0

    !     load initial particle information from external disk file
    if(config_input) then
        write(*,*)'  **************************************************'
        write(*,*)'       Reading initial particle configuration...   '
        !
        open(1,file='../data/ini_xv.dat',status='old')
        open(2,file='../data/ini_state.dat',status='old')
        open(3,file='../data/ini_other.dat',status='old')
        !
        ntotal = 0
        do while (.not.eof(1))
            ntotal = ntotal + 1
            i = ntotal
            read(1,*) tempim, (x(d,i),d = 1,dim), (vx(d,i),d = 1,dim)
            read(2,*) tempim, mass(i), rho(i), p(i), u(i)
            read(3,*) tempim, tempitype, hsml(i)
            im = NINT(tempim)
            itype(i)=NINT(tempitype)
            !          if(vx(1,i).eq.2500) then
            !             vx(1,i)=4000
            !          endif
            if(itype(i).eq.7) then
                hsml(i)=4*hsml(i) ! 1.5 particle per h in this problem feng
            else if(itype(i).eq.101) then
                hsml(i)=1.2*hsml(i)
            else if(itype(i).eq.200)then
                hsml(i)=1.2*hsml(i)
            endif
        enddo
        !
        write(*,*)'      Total number of particles   ', ntotal
        write(*,*)'  **************************************************'

        !
        close(1)
        close(2)
        close(3)

    else
        if (flowtype<3) then
            Ny=50
            Nx=Ny
            Lx=25
            Ly=25
            dx=Ly/Ny
            ntotal=Nx*Ny
            do i=1,ntotal
                jx=int(i/Ny)
                jy=mod(i,Ny)
                if (jy.eq.0.0) then
                    jy=Ny
                    jx=jx-1
                endif
                x(1,i)=(jx+0.5)*dx
                x(2,i)=(jy-0.5)*dx
                vx(1,i)=0.0
                vx(2,i)=0.0
                mass(i)=250
                rho(i)=1000
                p(i)=0.0
                u(i)=0.0
                itype(i)=200;
                hsml(i)=1.2*dx
                if (jx.eq.0.or.jx.eq.9.or.jx.eq.19.or.jx.eq.29.or.jx.eq.39.or.jx.eq.49.or.jy.eq.1.or.jy.eq.10.or.jy.eq.20.or.jy.eq.30.or.jy.eq.40.or.jy.eq.50) mass(i)=200
            enddo
            write(*,*) 'Current total number of particles = ', ntotal
        elseif(flowtype==3) then
            Ny=50
            Nx=25
            Lx=12.5
            Ly=25
            dx=Ly/Ny
            ntotal=Nx*Ny
            do i=1,ntotal
                jx=int(i/Ny)
                jy=mod(i,Ny)
                if (jy.eq.0.0) then
                    jy=Ny
                    jx=jx-1
                endif
                x(1,i)=(jx+0.5)*dx
                x(2,i)=(jy-0.5)*dx
                vx(1,i)=0.0
                vx(2,i)=0.0
                mass(i)=250
                rho(i)=1000
                p(i)=0.0
                u(i)=0.0
                itype(i)=200
                hsml(i)=1.2*dx
                if (jx.eq.0.or.jx.eq.9.or.jx.eq.19.or.jx.eq.29.or.jx.eq.39.or.jx.eq.49.or.jy.eq.1.or.jy.eq.10.or.jy.eq.20.or.jy.eq.30.or.jy.eq.40.or.jy.eq.50) mass(i)=200
            enddo
            write(*,*) 'Current total number of particles = ', ntotal
        elseif(flowtype==4) then ! Splash of water
            Ny=100
            Nx=250
            Lx=25
            Ly=10
            dx=Ly/Ny
            ntotal=Nx*Ny
            do i=1,ntotal
                jx=int(i/Ny)
                jy=mod(i,Ny)
                if (jy.eq.0.0) then
                    jy=Ny
                    jx=jx-1
                endif
                x(1,i)=(jx+0.5)*dx
                x(2,i)=(jy-0.5)*dx
                vx(1,i)=0.0
                vx(2,i)=0.0
                mass(i)=10
                rho(i)=1000
                p(i)=rho(i)*9.8*(Ly-x(2,i))
                u(i)=0.0
                itype(i)=200
                hsml(i)=1.2*dx
               ! if (jx.eq.0.or.jx.eq.9.or.jx.eq.19.or.jx.eq.29.or.jx.eq.39.or.jx.eq.49.or.jy.eq.1.or.jy.eq.10.or.jy.eq.20.or.jy.eq.30.or.jy.eq.40.or.jy.eq.50) mass(i)=200
            enddo
            
            ! generate the water drop in circle
            x0=12.5
            y0=12.5
            
            ntotal=ntotal+1
            x(1,ntotal)=x0
            x(2,ntotal)=y0
            vx(1,ntotal)=0
            vx(2,ntotal)=-2
            mass(ntotal)=10
            rho(ntotal)=1000
            p(ntotal)=0.0
            u(ntotal)=0.0
            itype(ntotal)=200
            hsml(ntotal)=1.2*dx
            
            do i=1,20
                do j=1,6*i
                    ntotal=ntotal+1
                    x(1,ntotal)=cosd(60*(j-1)/i)*i*dx+x0
                    x(2,ntotal)=sind(60*(j-1)/i)*i*dx+y0
                    vx(1,ntotal)=0
                    vx(2,ntotal)=-2
                    mass(ntotal)=10
                    rho(ntotal)=1000
                    p(ntotal)=0.0
                    u(ntotal)=0.0
                    itype(ntotal)=200
                    hsml(ntotal)=1.2*dx
                enddo
            enddo          
            write(*,*) 'Current total number of particles = ', ntotal
        endif
    endif

    rho0=rho(1)


    end