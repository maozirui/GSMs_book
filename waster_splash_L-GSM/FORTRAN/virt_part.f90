    subroutine virt_part1(itimestep, dx, ntotal,nvirt,hsml,mass,x,vx,rho,p,itype,stress)
    !----------------------------------------------------------------------
    !   Subroutine to determine the information of virtual particles
    !     itimestep : Current time step                                 [in]
    !     ntotal : Number of particles                                  [in]
    !     nvirt  : Number of virtual particles                         [out]
    !     hsml   : Smoothing Length                                 [in|out]
    !     mass   : Particle masses                                  [in|out]
    !     x      : Coordinates of all particles                     [in|out]
    !     vx     : Velocities of all particles                      [in|out]
    !     rho    : Density                                          [in|out]
    !     u      : internal energy                                  [in|out]
    !     itype   : type of particles                               [in|out]

    use config_parameter
    implicit none
    !
    real(8)  ntotal, nvirt1, itype(maxn),i,im,d, j,dd,nvirt2,nvirt
    real(8) hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),rho(maxn), u(maxn), p(maxn),stress(maxn, dim, dim)
    real(8) tempim,tempitype,dx,Nx,Ny,Lx,Ly,layer,L,nb
    integer(4) stat,itimestep

    if (vp_input) then

        open(1,file='../data/xv_vp.dat',status='old',iostat=stat)
        open(2,file='../data/state_vp.dat',status='old')
        open(3,file='../data/other_vp.dat',status='old')
        read(1,*) nvirt1
        !      nvirt1=0
        !      i = ntotal
        do j = 1, nvirt1
            i = ntotal + j

            !      do while (.not.is_iostat_end(1))
            !      do while (nvirt1.lt.1)
            !         nvirt1 = nvirt1 + 1
            !         i = i + 1
            read(1,*)tempim, (x(d, i),d = 1, dim), (vx(d, i),d = 1, dim)
            read(2,*)tempim, mass(i), rho(i), p(i)
            read(3,*)tempim, tempitype, hsml(i)
            !         im = NINT(tempim)
            itype(i)=NINT(tempitype)
            hsml(i)=hsml(i)*1.2 ! 1.5 particles per h in this problem
        enddo

        call no_slipbc(ntotal,nvirt1,x,hsml,vx,p,stress)

        !     Built-in Monaghan type virtual particle for particular probs
        !     Set the 'itype' to zero to be identified.
        close(1)
        close(2)
        close(3)

    else
        !!!!!!!!input virtual particles with a small obstacle.
        !!
        if (flowtype==1) then
            layer=1.0 ! layers of virtual particles in the boundary
            Lx=80
            Ly=30

            ! set left boundary
            nvirt1=0.0
            L=0.0
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(0.5-i)*dx
                    x(2,ntotal+nvirt1)=(int((nvirt1-1)/layer+1)-0.5)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo

            nb=nvirt1

            ! set bottom bondary
            L=0.0
            do while (L<50)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-1.5)*dx
                    x(2,ntotal+nvirt1)=(0.5-i)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-2.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(1,ntotal+nvirt1)
            enddo

            L=54.75
            do while (L<Lx)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(int((nvirt1-nb+5)/layer+1)-0.5)*dx
                    x(2,ntotal+nvirt1)=(0.5-i)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-2.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(1,ntotal+nvirt1)
            enddo


            !
            nvirt2=0
            do i=1,10
                nvirt2=nvirt2+1
                x(1,nvirt1+nvirt2+ntotal)=50.25+i*0.4
                x(2,nvirt1+nvirt2+ntotal)=i*0.3-0.25
                vx(1,nvirt1+nvirt2+ntotal)=0.0
                vx(2,nvirt1+nvirt2+ntotal)=0.0
                mass(ntotal+nvirt1+nvirt2)=250
                rho(ntotal+nvirt1+nvirt2)=1000
                p(ntotal+nvirt1+nvirt2)=0.0
                u(ntotal+nvirt1+nvirt2)=0.0
                itype(ntotal+nvirt1+nvirt2)=-3.0
                hsml(ntotal+nvirt1+nvirt2)=1.2*dx
                do d=i,dim
                    do dd=1,dim
                        stress(ntotal+nvirt1+nvirt2,dd,d)=0.0
                    enddo
                enddo
            enddo

            do i=1,5
                nvirt2=nvirt2+1
                x(1,nvirt1+nvirt2+ntotal)=54.25
                x(2,nvirt1+nvirt2+ntotal)=i*0.5-0.25
                vx(1,nvirt1+nvirt2+ntotal)=0.0
                vx(2,nvirt1+nvirt2+ntotal)=0.0
                mass(ntotal+nvirt1+nvirt2)=250
                rho(ntotal+nvirt1+nvirt2)=1000
                p(ntotal+nvirt1+nvirt2)=0.0
                u(ntotal+nvirt1+nvirt2)=0.0
                itype(ntotal+nvirt1+nvirt2)=-3.0
                hsml(ntotal+nvirt1+nvirt2)=1.2*dx
                do d=i,dim
                    do dd=1,dim
                        stress(ntotal+nvirt1+nvirt2,dd,d)=0.0
                    enddo
                enddo
            enddo



            nvirt=nvirt1+nvirt2
        elseif (flowtype==2) then
            !!!!! containning right boundary
            layer=1.0 ! layers of virtual particles in the boundary
            Lx=80
            Ly=50

            ! set left boundary
            nvirt1=0.0
            L=0.0
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(0.5-i)*dx
                    x(2,ntotal+nvirt1)=(int((nvirt1-1)/layer+1)-0.5)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo
            L=0.0
            nb=nvirt1
            ! set right wall
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=80.25
                    x(2,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-0.5)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo

            nb=nvirt1
            ! set bottom bondary
            L=0.0
            do while (L<Lx)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-1.5)*dx
                    x(2,ntotal+nvirt1)=(0.5-i)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-2.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(1,ntotal+nvirt1)
            enddo
            nvirt=nvirt1
        elseif (flowtype==3) then
            layer=1.0 ! layers of virtual particles in the boundary
            Lx=12.5
            Ly=30
            ! set left boundary
            nvirt1=0.0
            L=0.0
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(0.5-i)*dx
                    x(2,ntotal+nvirt1)=(int((nvirt1-1)/layer+1)-0.5)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo
            L=0.0
            nb=nvirt1
            ! set right wall
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=Lx+0.5*dx
                    x(2,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-0.5)*dx+2.0
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo
            
            nb=nvirt1
            ! set bottom bondary
            L=0.0
            do while (L<Lx+2)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-0.5)*dx
                    x(2,ntotal+nvirt1)=(0.5-i)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=250
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=0.0
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-2.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(1,ntotal+nvirt1)
            enddo
            
            do i=1,6
                nvirt1=nvirt1+1
                if (i<6) then
                    x(1,ntotal+nvirt1)=12.75+(i-1)*dx
                    x(2,ntotal+nvirt1)=1.75
                elseif (i==6) then
                    x(1,ntotal+nvirt1)=-0.25
                    x(2,ntotal+nvirt1)=-0.25
                endif
                vx(1,ntotal+nvirt1)=0.0
                vx(2,ntotal+nvirt1)=0.0
                mass(ntotal+nvirt1)=250
                rho(ntotal+nvirt1)=1000
                p(ntotal+nvirt1)=0.0
                u(ntotal+nvirt1)=0.0
                if (i<4) then
                    itype(ntotal+nvirt1)=-2.0
                else
                    itype(ntotal+nvirt1)=-3.0
                endif 
                hsml(ntotal+nvirt1)=1.2*dx
                do d=i,dim
                    do dd=1,dim
                        stress(ntotal+nvirt1,dd,d)=0.0
                    enddo
                enddo
            enddo
        elseif (flowtype==4) then   ! Splash of water
            layer=1.0 ! layers of virtual particles in the boundary
            Lx=25
            Ly=15
            ! set left boundary
            nvirt1=0.0
            L=0.0
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(0.5-i)*dx
                    x(2,ntotal+nvirt1)=(int((nvirt1-1)/layer+1)-0.5)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=10
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=rho(ntotal+nvirt1)*9.8*(10.0-x(2,ntotal+nvirt1))
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo
            L=0.0
            nb=nvirt1
            ! set right wall
            do while (L<Ly)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=Lx+0.5*dx
                    x(2,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-0.5)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=10
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=rho(ntotal+nvirt1)*9.8*(10.0-x(2,ntotal+nvirt1))
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-1.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(2,ntotal+nvirt1)+dx
            enddo
            
            nb=nvirt1
            ! set bottom bondary
            L=0.0
            do while (L<Lx-0.5*dx)
                do i=1,layer
                    nvirt1=nvirt1+1
                    x(1,ntotal+nvirt1)=(int((nvirt1-nb-1)/layer+1)-0.5)*dx
                    x(2,ntotal+nvirt1)=(0.5-i)*dx
                    vx(1,ntotal+nvirt1)=0.0
                    vx(2,ntotal+nvirt1)=0.0
                    mass(ntotal+nvirt1)=10
                    rho(ntotal+nvirt1)=1000
                    p(ntotal+nvirt1)=rho(ntotal+nvirt1)*9.8*(10.0-x(2,ntotal+nvirt1))
                    u(ntotal+nvirt1)=0.0
                    itype(ntotal+nvirt1)=-2.0
                    hsml(ntotal+nvirt1)=1.2*dx
                    do d=i,dim
                        do dd=1,dim
                            stress(ntotal+nvirt1,dd,d)=0.0
                        enddo
                    enddo
                enddo
                L=x(1,ntotal+nvirt1)
            enddo
            
            do i=1,2
                nvirt1=nvirt1+1
                x(1,ntotal+nvirt1)=(i-1)*(Lx+dx)-0.5*dx
                x(2,ntotal+nvirt1)=-dx/2
                vx(1,ntotal+nvirt1)=0.0
                vx(2,ntotal+nvirt1)=0.0
                mass(ntotal+nvirt1)=10
                rho(ntotal+nvirt1)=1000
                p(ntotal+nvirt1)=rho(ntotal+nvirt1)*9.8*(10.0-x(2,ntotal+nvirt1))
                u(ntotal+nvirt1)=0.0
                itype(ntotal+nvirt1)=-3.0
                hsml(ntotal+nvirt1)=1.2*dx
                do d=i,dim
                    do dd=1,dim
                        stress(ntotal+nvirt1,dd,d)=0.0
                    enddo
                enddo
            enddo    
                
        endif
    endif
    nvirt=nvirt1+nvirt2
        call no_slipbc(ntotal,nvirt,x,hsml,vx,p,stress,itype,rho)
    end
    !
    !      subroutine virt_part2(itimestep, ntotal,nvirt2,hsml,mass,x,vx,rho,u,p,itype)
    !!----------------------------------------------------------------------
    !!   Subroutine to determine the information of virtual particles
    !!     itimestep : Current time step                                 [in]
    !!     ntotal : Number of particles                                  [in]
    !!    `
    !!     mass   : Particle masses                                  [in|out]
    !!     x      : Coordinates of all particles                     [in|out]
    !!     vx     : Velocities of all particles                      [in|out]
    !!     rho    : Density                                          [in|out]
    !!     u      : internal energy                                  [in|out]
    !!     itype   : type of particles                               [in|out]
    !      use config_parameter
    !      implicit none
    !!
    !      real(8)  ntotal, nvirt2, itype(maxn),i,im,d
    !      real(8) hsml(maxn),mass(maxn),x(dim,maxn),vx(dim,maxn),rho(maxn), u(maxn), p(maxn)
    !      integer(4) itimestep
    !      integer(4) stat
    !!
    !
    !   if (vp_input) then
    !
    !      open(1,file='../data/xv_vp.dat',status='old',iostat=stat)
    !      open(2,file='../data/state_vp.dat',status='old')
    !      open(3,file='../data/other_vp.dat',status='old')
    !      read(1,*) nvirt2
    !!      nvirt2=ntotal
    !	do i = ntotal+1, ntotal+nvirt2
    !	  x(1, i) = x(1, i-ntotal )
    !        x(2, i) = -x(2, i-ntotal )
    !        vx(1, i) = -vx(1, i-ntotal )
    !	  vx(2, i) = -vx(2, i-ntotal )
    !!       z-directional coordinate/velocity are mirrored in this case:
    !!        x(3, i)  = -x(3, i-ntotal )
    !!        vx(3, i) = -vx(3, i-ntotal )
    !
    !	  rho (i)  =  rho( i-ntotal )
    !	  mass(i)  = mass( i-ntotal )
    !	  p(i)     = p( i-ntotal )
    !
    !	  itype(i) = itype( i-ntotal )
    !	  hsml(i)  =   hsml( i-ntotal )
    !      enddo
    !      endif
    !!
    !      if (mod(itimestep,print_step).eq.0) then
    !        if (int_stat) then
    !         print *,' >> Statistics: Virtual boundary particles:'
    !         print *,'          Number of Type2 virtual particles:',nvirt2
    !        endif
    !      endif
    !!
    !      end
