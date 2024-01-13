    subroutine single_step(itimestep, dt, ddx,ntotal, hsml, mass, x, x_min,vx, vx_mid,  rho, rho_mid, p, t, dx, dvx, drho,itype, av, nvirt1,nvirt2,rho0,AArea)
    !----------------------------------------------------------------------
    !   Subroutine to determine the right hand side of a differential
    !   equation in a single step for performing time integration

    !   In this routine and its subroutines the SPH algorithms are performed.
    !     itimestep: Current timestep number                            [in]
    !     dt       : Timestep                                           [in]
    !     ntotal   :  Number of particles                               [in]
    !     hsml     :  Smoothing Length                                  [in]
    !     mass     :  Particle masses                                   [in]
    !     x        :  Particle position                                 [in]
    !     vx       :  Particle velocity                                 [in]
    !     u        :  Particle internal energy (specific)               [in]
    !     s        :  Particle entropy (not used here)                  [in]
    !     rho      :  Density                                       [in/out]
    !     p        :  Pressure                                         [out]
    !     t        :  Temperature                                   [in/out]
    !     tdsdt    :  Production of viscous entropy t*ds/dt            [out]
    !     dx       :  dx = vx = dx/dt                                  [out]
    !     dvx      :  dvx = dvx/dt, force per unit mass                [out]
    !     du       :  du  = du/dt                                      [out]
    !     ds       :  ds  = ds/dt                                      [out]
    !     drho     :  drho =  drho/dt                                  [out]
    !     itype    :  Type of particle                                 [in]
    !     av       :  Monaghan average velocity                        [out]
    use config_parameter
    implicit none
    !
    real(8) ntotal, itype(maxn)
    real(8) dt, hsml(maxn), mass(maxn), x(dim,maxn),x_min(dim, maxn)
    real(8) vx(dim,maxn), u(maxn), s(maxn), rho(maxn), p(maxn)
    real(8) t(maxn), tdsdt(maxn), dx(dim,maxn), dvx(dim,maxn)
    real(8) du(maxn), ds(maxn), drho(maxn), av(dim, maxn),B
    real(8) rho_mid(maxn), vx_mid(dim,maxn), du_pl, du_vis, du_heat, f_R(maxn),stress(maxn, dim, dim),stress_prerate(maxn, dim, dim),p_prerate(maxn),plstic_du(maxn),du_p,p_du(maxn),du_plstic

    real(8) i, d, nvirt1, nvirt2, niac, pair_i(max_interaction),xxx,totalstrain(maxn),e_traceless(maxn,dim,dim),alpha,AArea(maxn)
    real(8) pair_j(max_interaction), ns(maxn), j, k, nvirt,gg(maxn)
    real(8) w(max_interaction), dwdx(dim,max_interaction)
    real(8) indvxdt(dim,maxn),exdvxdt(dim,maxn),ardvxdt(dim,maxn),arstrdvxdt(dim,maxn),err,ddx,d_x,d_y,v_x,v_y,a_x,a_y
    real(8) avdudt(maxn), ahdudt(maxn), c(maxn),rho0,acc,min_dis,x1,v1,a1,near(maxn)
    integer(4) itimestep,numb(maxn,15),countiac(maxn)

    if (itimestep.eq.1) write(*,'(A)') 'Variables in single_step has been defined!'
    
    
    !     Initialization
    do i=1,ntotal

        do d=1,dim
            indvxdt(d,i) = 0.
            ardvxdt(d,i) = 0.
            exdvxdt(d,i) = 0.
        enddo
    enddo
    
    ! update the pressure
    if (itimestep.ne.1) then
        do i=1,ntotal
            B=500*9.8*250*1000/7.0
            p(i)=B*((rho(i)/rho0)**7.0-1.0)
        enddo
    endif
    
    !---  Generation of Monaghan virtual (boundary) particles:

    nvirt1 = 0
    if (virtual_part1) then
        call virt_part1(itimestep,ddx, ntotal,nvirt1,hsml,mass,x,vx,rho,p,itype,stress)
        !     Mid-step (time) density defined for solids due to strain rate model
        do i = ntotal+1, ntotal+nvirt1
            do d = 1, dim
                vx_mid(d, i) = 0.
            enddo
            rho_mid(i) = rho(i)
        enddo
    endif
    !      endif

    !    xxx=itype(2500)
    !---  Generation of Libersky virtual (reflected) particles:
    nvirt2 = 0
    if (virtual_part2) then
        !        call virt_part2(itimestep, ntotal, nvirt2,hsml,mass,x,vx,rho,p,itype)
        !     Mid-step (time) density defined for solids due to strain rate model
        do i = ntotal+nvirt1+1, ntotal+nvirt1+nvirt2
            vx_mid(1, i) = vx_mid(1, (i-(ntotal+nvirt1)))
            vx_mid(2, i) = -vx_mid(2, (i-(ntotal+nvirt1)))
            rho_mid(i) = rho_mid(i-(ntotal+nvirt1))
        enddo
    endif
    !---  Interaction parameters, calculating neighboring particles
    !     and optimzing smoothing length
    !     direct_find和tree search未改动,link list已改动
    niac=0
    !         call CPU_TIME(time)
    !         write(*,*)'time before searching=', time
    ! Define the pressure variable



    if (nnps.eq.1) then
        call direct_find(itimestep, ntotal+nvirt1,hsml,x,niac,pair_i,pair_j,w,dwdx,ns,min_dis)
        if (itimestep.eq.1)   write(*,'(A)') 'direct_find has been called!'
    else if (nnps.eq.2) then
        call link_list(itimestep,ntotal+nvirt,hsml,x,niac,pair_i,pair_j,w,dwdx,ns,hsml,itype)
        if (itimestep.eq.1)   write(*,'(A)') 'link_list has been called!'
    else if (nnps.eq.3) then
        call tree_search(itimestep, ntotal+nvirt,hsml,x,niac,pair_i,pair_j,w,dwdx,ns)
        if (itimestep.eq.1) write(*,'(A)') 'tree_search has been called!'
    else if (nnps.eq.4) then
        call direct_GSM(itimestep,hsml,ntotal,nvirt1,nvirt2,x,vx,p,mass,rho,itype,niac,pair_i,pair_j,w,dwdx,AArea,numb,countiac)
        if (itimestep.eq.1) write(*,'(A)') 'GSM scheme has been called!'
    endif
    
    
        
    !         call CPU_TIME(time)
    !         write(*,*)'time after searching=', time


    !     kernel gradient correct
    if (kernel_correct) then
        call  kernelgrad_correct(ntotal,niac,x,hsml,pair_i,pair_j,itype,rho,mass,dwdx)
        write(*,'(A)') 'kernelgrad_correct has been called!'
    end if

    !---  Density approximation or change rate
    !     未改动
    if (summation_density) then
        call sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,itype,rho)
        write(*,'(A)') 'sum_density has been called!'
    else
        call con_density(ntotal,mass,niac,pair_i,pair_j,dwdx,vx, itype,x,rho, drho,AArea)
        if (itimestep.eq.1)  write(*,'(A)') 'con_density has been called!'
    endif
    !---  Internal forces:

    !        call CPU_TIME(time)
    !         write(*,*)'time before internal force=', time

    call int_force(itimestep,dt,ntotal,hsml,mass,vx, vx_mid,niac,rho, rho_mid,rho0, pair_i,pair_j,dwdx,itype,x,x_min,t,c,p,indvxdt,AArea,near)

    !         call CPU_TIME(time)
    !         write(*,*)'time after internal force=', time

    if (itimestep.eq.1) write(*,'(A)') 'int_force has been called!'


    !---  Artificial viscosity:
    if (visc_artificial) then
        call art_visc(ntotal,itype,hsml,mass,x,vx,niac,rho,pair_i,pair_j,w,dwdx,ardvxdt,avdudt,alpha,AArea)
        if (itimestep.eq.1)   write(*,'(A)') 'art_visc has been called!'
    end if

    if(stress_artificial)then
        call art_stress(ntotal+nvirt,hsml,mass,x,vx,niac,rho,c, pair_i,pair_j,itype,w,dwdx,stress,p,arstrdvxdt)
        if (itimestep.eq.1) write(*,'(A)') 'art_stress has been called!'
    endif
    !---  External forces:
    !     未改动
    if (ex_force) then
        call ext_force(itimestep,ntotal,nvirt1,nvirt2,mass,x,niac,pair_i,pair_j,itype, hsml, exdvxdt,vx)
        !
        if (itimestep.eq.1)   write(*,'(A)') 'ext_force has been called!'
    endif
    !---  Artificial heat:
    if (heat_artificial) then
        call art_heat(ntotal,hsml, mass,x,vx,niac,rho,u,c,pair_i,pair_j,w,dwdx,ahdudt)
        if (itimestep.eq.1)    write(*,'(A)') 'art_heat has been called!'
    end if
    !---  Artificial compressibility:
    !     XSPH: average velocity of each partile for avoiding penetration
    if (average_velocity) then
        call av_vel(ntotal,mass,itype,niac,pair_i,pair_j, w,dwdx, x,vx, rho, av,countiac,near)
        !      write(*,'(A)') 'av_vel has been called!'
    end if
    !     Calculating the neighboring particles and undating HSML
    if (sle.ne.0) then
        call h_upgrade(dt,ntotal,ns, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)
        write(*,'(A)') 'h_upgrade has been called!'
    end if
    !---  ÃÜ¶ÈÐÞÕý£¬Ã¿¸ô30²½
    if (nor_density) then
        if (mod(itimestep,50)==0) then
            call sum_density(ntotal+nvirt,hsml,mass,niac,pair_i,pair_j,w,itype,rho)
            write(*,'(A)') 'sum_density has been called!'
        end if
    end if

    acc=0
    do i=1,ntotal
        do d=1,dim
            dvx(d,i) = indvxdt(d,i) + exdvxdt(d,i) + ardvxdt(d,i)+arstrdvxdt(d,i)
            if(indvxdt(d,i).gt.1e30) then
                !             write(*,'(A)') 'single step:indvxdt is Infinity!!'
                !             pause
                !          else if(ardvxdt(d,i).gt.1e30) then
                write(*,'(A)') 'single step:indvxdt is Infinity!!'
                pause
            endif
        enddo
        
    enddo
    
    ! update the critical time step
    dt=1e-5
    !do i=1,ntotal
    !    do j=1,countiac(i)
    !        k=numb(i,j)
    !        x1=sqrt((x(1,i)-x(1,k))**2+(x(2,i)-x(2,k))**2)
    !        v1=sqrt((vx_mid(1,i)-vx_mid(1,k))**2+(vx_mid(2,i)-vx_mid(2,k))**2)
    !        a1=sqrt((dvx(1,i)-dvx(1,k))**2+(dvx(2,i)-dvx(2,k))**2)
    !        if (0.05*x1/v1<dt) dt=0.05*x1/v1
    !        if (sqrt(0.05*x1/a1)<dt) dt=sqrt(0.05*x1/a1)
    !    enddo
    !enddo
    !
    !do i=ntotal+1,ntotal+nvirt1
    !    do j=1,ntotal
    !        if (itype(j)<200) then
    !            x1=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)
    !            v1=sqrt((vx_mid(1,j))**2+(vx_mid(2,j))**2)
    !            a1=sqrt((dvx(1,j))**2+(dvx(2,j))**2)
    !            if (0.05*x1/v1<dt) dt=0.05*x1/v1
    !            if (sqrt(0.05*x1/a1)<dt) dt=sqrt(0.05*x1/a1)
    !        endif 
    !    enddo
    !enddo
    !
    !do i=1,ntotal
    !    if (itype(i)<200) then
    !        do j=1,ntotal
    !            if (itype(j)<200) then
    !               x1=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)
    !               v1=sqrt((vx_mid(1,i)-vx_mid(1,j))**2+(vx_mid(2,i)-vx_mid(2,j))**2)
    !               a1=sqrt((dvx(1,i)-dvx(1,j))**2+(dvx(2,i)-dvx(2,j))**2)
    !               if (0.05*x1/v1<dt) dt=0.05*x1/v1
    !               if (sqrt(0.05*x1/a1)<dt) dt=sqrt(0.05*x1/a1)
    !            endif
    !        enddo
    !        
    !    endif
    !enddo
    !
            

            
    if (nnps==1) dt=5.0e-5

    !      call CPU_TIME(time)
    !         write(*,*)'time after forces calculating =', time

    end