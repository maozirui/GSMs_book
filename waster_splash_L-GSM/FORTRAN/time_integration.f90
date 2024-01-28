    subroutine time_integration(x,ddx,vx, mass, rho, p, itype,hsml, ntotal, maxtime, dt,rho0 )
    !----------------------------------------------------------------------
    !      x-- coordinates of particles                       [input/output]
    !      vx-- velocities of particles                       [input/output]
    !      mass-- mass of particles                                  [input]
    !      rho-- densities of particles                       [input/output]
    !      p-- pressure  of particles                         [input/output]
    !      u-- internal energy of particles                   [input/output]
    !      c-- sound velocity of particles                          [output]
    !      s-- entropy of particles, not used here                  [output]
    !      e-- total energy of particles                            [output]
    !      itype-- types of particles                               [input]
    !           =1   ideal gas
    !           =2   water
    !           =3   tnt
    !      hsml-- smoothing lengths of particles              [input/output]
    !      ntotal-- total particle number                            [input]
    !      maxtimestep-- maximum timesteps                           [input]
    !      dt-- timestep                                             [input]
    !      f_R : reaction progress (0 means no reaction, 1 means total reaction)

    use config_parameter
    implicit none
    !
    real(8) itype(maxn), ntotal, maxtime
    real(8) x(dim, maxn), vx(dim, maxn), mass(maxn)
    real(8) rho(maxn), p(maxn), u(maxn), c(maxn), s(maxn), e(maxn),stress(maxn, dim, dim),stress_prerate(maxn, dim, dim),p_prerate(maxn),stress_mid(maxn, dim, dim),p_mid(maxn)
    real(8) hsml(maxn), dt, f_R(maxn)
    !
    real(8)  dx(dim,maxn), dvx(dim, maxn), du(maxn),e_traceless(maxn,dim,dim)
    real(8)   drho(maxn),  av(dim, maxn), ds(maxn)
    real(8)         t(maxn), tdsdt(maxn)
    real(8)  x_min(dim, maxn), v_min(dim, maxn)
    real(8)            rho_min(maxn), u_min(maxn)
    real(8)  time, temp_rho, temp_u,internal_plstic,internal_artvis,internal_p,ahdudt(maxn),AArea(maxn)
    real(8) i, j, k,d,d1,d2,gg(maxn),ddx,nvirt
    real(8),save :: current_ts=0, nstart=0, period=0, max_time=0
    real(8) d_i_p, d_i_h, d_i_v,d_i_plstic,d_i_pl,accumulated_strain(maxn)
    real(8) internal_heat,nvirt1,nvirt2,err,rho0
    integer(4)  itimestep,no
    !
    do i = 1, ntotal
        do d = 1, dim
            av(d, i) = 0.
            x_min(d,i)=x(d,i)
            AArea(i)=0.0
            mass(i)=250
        enddo
        rho_min(i) = rho(i) ! new to solids
    enddo
    !
    if (period.eq.0) then 
        no=0
    else
        no=int((period)/save_period)
    endif
    
    
    do itimestep = nstart+1, nstart+100000000  !call single_step in every timestep
        !
        current_ts=current_ts+1


        !     If not first time step, then update  density and
        !     velocity half a time step
        if (itimestep .ne. 1) then
            do i=1,ntotal
                rho_min(i) = rho(i)
                temp_rho=0.
                    if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
                    rho(i) = rho(i) +(dt/2.)*( drho(i)+ temp_rho)
                !
                do d = 1, dim
                    v_min(d, i) = vx(d, i)
                    vx(d, i) = vx(d, i) + (dt/2.)*dvx(d, i)
                enddo
            enddo

        endif


        !---  Definition of variables out of the function vector:
        ! Define the intial stress tensor

        call single_step(itimestep, dt,ddx, ntotal, hsml, mass, x, x_min,vx, v_min, rho, rho_min, p,  t, dx, dvx, drho,itype, av, nvirt1,nvirt2,rho0,AArea)
        
        period=period+dt  

        if (mod(itimestep,print_step).eq.0) then
            call CPU_TIME(time)
            !         write(*,*)'______________________________________________'
            write(*,*)'  Progress    percentage =',     period/(maxtime+max_time)*100
            write(*,*)'Time remaining = ', time/period*((maxtime+max_time)-period)
        endif

        if (itimestep .eq. 1) then
            do i=1,ntotal
                
                if (.not.summation_density ) then
                    temp_rho=0.
                    if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
                    rho(i) = rho(i) + (dt/2.)* (drho(i)+temp_rho)
                endif
                
                do d = 1, dim
                    vx(d, i) = vx(d, i) + (dt/2.) * dvx(d, i) + av(d, i)

                    x(d, i) = x(d, i) + dt * vx(d, i)
                enddo
                

            enddo
            !
        else   ! NOT the first step
            do i=1,ntotal
                if (.not.summation_density ) then
                    temp_rho=0.
                    if (dim.eq.1) temp_rho=-nsym*rho(i)*vx(1,i)/x(1,i)
                    rho(i) = rho_min(i) + dt*(drho(i)+temp_rho)
                endif
                !!!!!!!  integration

                do d = 1, dim
                    vx(d, i) = v_min(d, i) + dt * dvx(d, i) + av(d, i)
                    x_min(d, i) = x(d, i)
                    x(d, i) = x(d, i) + dt * vx(d, i)
                enddo
            enddo
        endif

        
        time = time + dt
        
       !do i=ntotal+1,ntotal+nvirt1
       !     vx(1,i)=0.0
       !     vx(2,i)=0.0
       !     p(i)=0.0
       !     rho(i)=1000.0
       ! enddo
        
        !Save step result
        if (int(period/save_period)>no) then
            no=no+1
            call output(x, vx, mass, rho, p, c,itype, hsml, ntotal+nvirt1, itimestep, period)
        endif
        !        call CPU_TIME(time)
        !         write(*,*)'current final time=', time
        
        if (period>maxtime+max_time) goto 145

    enddo
    !
145    nstart=current_ts
    max_time=maxtime+max_time

    !
      
    end