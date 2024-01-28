    subroutine int_force(itimestep,dt,ntotal,hsml,mass,vx,vx_mid,niac,rho,rho_mid,rho0,pair_i,pair_j,dwdx,itype,x,x_mid,t,c,p,dvxdt,AArea,near)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the internal forces on the right hand side
    !   of the Navier-Stokes equations, i.e. the pressure gradient and the
    !   gradient of the viscous stress tensor, used by the time integration.
    !   Moreover the entropy production due to viscous dissipation, tds/dt,
    !   and the change of internal energy per mass, de/dt, are calculated.

    !     itimestep: Current timestep number                            [in]
    !     dt     :   Time step                                          [in]
    !     ntotal : Number of particles                                  [in]
    !     hsml   : Smoothing Length                                     [in]
    !     mass   : Particle masses                                      [in]
    !     vx     : Velocities of all particles                          [in]
    !     niac   : Number of interaction pairs                          [in]
    !     rho    : Density                                              [in]
    !     eta    : Dynamic viscosity                                    [in]
    !     pair_i : List of first partner of interaction pair            [in]
    !     pair_j : List of second partner of interaction pair           [in]
    !     dwdx   : Derivative of kernel with respect to x, y and z      [in]
    !     itype  : Type of particle (material types)                    [in]
    !     u      : Particle internal energy                             [in]
    !     x      : Particle coordinates                                 [in]
    !     itype  : Particle type                                        [in]
    !     t      : Particle temperature                             [in/out]
    !     c      : Particle sound speed                                [out]
    !     p      : Particle pressure                                   [out]
    !     dvxdt  : Acceleration with respect to x, y and z             [out]
    !     dedt   : Change of specific internal energy                  [out]
    use config_parameter
    implicit none
    !
    real(8)  ntotal,niac,pair_i(max_interaction),pair_j(max_interaction), itype(maxn),AArea(maxn)
    real(8) dt
    real(8) hsml(maxn), mass(maxn), vx(dim,maxn),x_mid(dim,maxn)
    real(8)  vx_mid(dim, maxn), rho_mid(maxn),rho(maxn)
    real(8)  dwdx(dim,max_interaction),u(maxn),x(dim,maxn), t(maxn)
    real(8)  c(maxn), p(maxn),dvxdt(dim,maxn),dedt(maxn),f_R(maxn)
    real(8) i, j, k, d, err,gg(maxn),I3
    real(8) e(ntotal, dim, dim), r(ntotal, dim, dim),e_traceless(ntotal, dim, dim),stress(maxn, dim, dim),totalstrain(maxn),stress_prerate(maxn, dim, dim),p_prerate(maxn)
    real(8) p_dvxdt(dim,maxn)
    real(8) v_dvxdt(dim,maxn),rho0,B,near(maxn)
    integer(4) itimestep
    !----     Option of consideration of plastic work effects

    if (itimestep.eq.1) write(*,'(A)') 'Variables in int_force have been defined!'
    dvxdt(1:dim, 1:maxn)=0.
    !	call pressure (ntotal, itype, rho, u, p, c, x, ttt, dt)  !feng


    !   calculating stress
    !      call  shear_stress(itimestep, ntotal, itype,mass, rho, rho_mid,x_mid,vx, vx_mid, dt,u, niac, pair_i, pair_j, dwdx, e_traceless,stress, t, p,stress_prerate,p_prerate,plstic_dedt,gg)
    !
    !
    !
    !    if (itimestep.eq.1)  write(*,'(A)') 'shear_stress has been called!'
    !
    !!     SPH summation of viscous force and viscous energy
    !!     visc is true, so the related programs should be rewriten
    !
    !      err=stress_prerate(1, 1, 2)
    !      err=stress_prerate(1, 2, 1)
    !      err=stress(1, 1, 2)
    !      err=stress(1, 2, 1)
    !
    !      if (visc) then
    !         call viscous_forces(ntotal, itype, niac, pair_i, pair_j,dwdx, mass, rho, vx, stress, v_dvxdt) !v_dvxdt is beyond limit of number
    !      if (itimestep.eq.1)   write(*,'(A)') 'viscous_forces has been called!'
    !
    !      endif
    !
    !
    !!
    !
    !
    
   
    
    !     SPH summation of Pressure force and pressure work
    call SPH_pressure(ntotal, itype, niac, pair_i, pair_j,dwdx, p, mass, rho, vx, p_dvxdt,AArea,x,near) !p_dvxdt is beyond limit of number
    if (itimestep.eq.1) write(*,'(A)') 'SPH_pressure has been called!'
    !!     Combination of pressure and viscous effects:
    do i=1,ntotal
        do d = 1, dim
            dvxdt(d, i) = p_dvxdt(d, i)
        enddo
    enddo


    
    end