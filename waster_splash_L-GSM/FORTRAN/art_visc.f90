    subroutine art_visc(ntotal,itype,hsml,mass,x,vx,niac,rho, pair_i,pair_j,w,dwdx,dvxdt,dedt,para,AArea)
    !----------------------------------------------------------------------
    !     Subroutine to calculate the artificial viscosity (Monaghan, 1992)

    !     ntotal : Number of particles (including virtual particles)    [in]
    !     hsml   : Smoothing Length                                     [in]
    !     mass   : Particle masses                                      [in]
    !     x      : Coordinates of all particles                         [in]
    !     vx     : Velocities of all particles                          [in]
    !     niac   : Number of interaction pairs                          [in]
    !     rho    : Density                                              [in]
    !     c      : sound speed                                          [in]
    !     pair_i : List of first partner of interaction pair            [in]
    !     pair_j : List of second partner of interaction pair           [in]
    !     w      : Kernel for all interaction pairs                     [in]
    !     dwdx   : Derivative of kernel with respect to x, y and z      [in]
    !     dvxdt  : Acceleration with respect to x, y and z             [out]
    !     dedt   : Change of specific internal energy                  [out]
    use config_parameter
    implicit none
    !
    real(8) ntotal, niac, pair_i(max_interaction),pair_j(max_interaction)
    real(8) hsml(maxn), mass(maxn), x(dim,maxn),vx(dim,maxn),rho(maxn), c(maxn), w(max_interaction),dwdx(dim,max_interaction), dvxdt(dim,maxn), dedt(maxn),AArea(maxn)
    real(8) i,j,k,d,alph
    real(8) dx, dvx(dim), alpha, beta, etq, piv,muv, vr, rr, h, mc, mrho, mhsml,itype(maxn),para
    !     Parameters for the artificial viscosity:
    !     Bulk viscosity
    parameter( alpha = 0.25)
    !     Shear viscosity
    parameter( beta  = 0.0)
    !     Parameter to avoid singularities
    parameter( etq   = 0.1)
    !
    para=alpha

    dvxdt(1:dim,1:maxn) = 0.

    !     Calculate SPH sum for artificial viscosity
    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        if (itype(j)<0.and.nnps.ne.4) goto 11
        mhsml= (hsml(i)+hsml(j))/2.
        vr = 0.e0
        rr = 0.e0
        do d=1,dim
            dvx(d) = vx(d,i) - vx(d,j)
            dx     =  x(d,i) -  x(d,j)
            vr     = vr + dvx(d)*dx
            rr     = rr + dx*dx
        enddo
        !     Artificial viscous force only if v_ij * r_ij < 0
        if (vr.lt.0.e0) then
            !     Calculate muv_ij = hsml v_ij * r_ij / ( r_ij^2 + hsml^2 etq^2 )
            muv = 0.5*vr/(rr + 0.25*etq*etq)
            !     Calculate PIv_ij = (-alpha muv_ij c_ij + beta muv_ij^2) / rho_ij
            !          £¡Must Be Assigned

            mrho = 0.5e0*(rho(i) + rho(j))
            mc=sqrt(200*9.8*25.0)
            piv  = (beta*muv - alpha*mc)*muv/mrho 
            !     Calculate SPH sum for artificial viscous force
            do d=1,dim
                if (nnps.eq.4) then ! L-GSM
                    if (flowtype==3) then ! discharge flow
                        if (itype(j)<0) then
                            if ((-0.25<=x(1,i)<=12.75.and.-0.25<=x(2,i)<=30).or.(12.75<=x(1,i)<=14.75.and.-0.25<=x(2,i)<=1.75)) then 
                                goto 27
                            else
                                goto 11
                            endif
                        endif
                    endif
27                  h=-piv*dwdx(d,k)*AArea(i)/mass(i)*rho(i)
                    if (d==1.and.itype(j)==-2) h=0.0
                    if (d==2.and.itype(j)==-1) h=0.0
                    dvxdt(d,i) = dvxdt(d,i) + mrho*h
                    
                else ! SPH
                    h = -piv*dwdx(d,k)
                    dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
                    dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
                endif


                if(isNaN(dwdx(d,k))) then
                    write(*,'(A)') 'art_visc:dwdx is NaN!'
                    pause
                endif
            enddo
            if(isNaN(muv)) then
                write(*,'(A)') 'art_visc:muv is NaN!'
                pause
            endif
            if(isNaN(mc)) then
                write(*,'(A)') 'art_visc:mc is NaN!'
                pause
            endif
            if(isNaN(c(i)).or.isNaN(c(j))) then
                write(*,'(A)') 'art_visc:c(i) or c(j) is NaN!'
                pause
            endif
        endif
11  enddo

    !    Change of specific internal energy:


    end