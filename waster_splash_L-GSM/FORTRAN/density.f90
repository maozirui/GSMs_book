    subroutine sum_density(ntotal,hsml,mass,niac,pair_i,pair_j,w,itype,rho)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the density with SPH summation algorithm.

    !     ntotal : Number of particles                                  [in]
    !     hsml   : Smoothing Length                                     [in]
    !     mass   : Particle masses                                      [in]
    !     niac   : Number of interaction pairs                          [in]
    !     pair_i : List of first partner of interaction pair            [in]
    !     pair_j : List of second partner of interaction pair           [in]
    !     w      : Kernel for all interaction pairs                     [in]
    !     itype   : type of particles                                   [in]
    !     x       : Coordinates of all particles                        [in]
    !     rho    : Density                                             [out]
    use config_parameter
    implicit none
    !
    real(8) ntotal, niac, pair_i(max_interaction),pair_j(max_interaction), itype(maxn)
    real(8) hsml(maxn),mass(maxn), w(max_interaction),rho(maxn)
    real(8) i, j, k, d
    real(8) selfdens, hv(dim), r, wi(maxn)
    !     wi(maxn)---integration of the kernel itself
    do d=1,dim
        hv(d) = 0.e0
    enddo
    !     Self density of each particle: Wi(i) (Kernel for distance 0)
    !     and take contribution of particle itself:
    r=0.
    !     Firstly calculate the integration of the kernel over the space
    do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv) !selfdens->w: Kernel for all interaction pairs
        wi(i)=selfdens*mass(i)/rho(i)
    enddo
    !
    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        wi(i) = wi(i) + mass(j)/rho(j)*w(k)
        wi(j) = wi(j) + mass(i)/rho(i)*w(k)
    enddo
    !     Secondly calculate the rho integration over the space
    do i=1,ntotal
        call kernel(r,hv,hsml(i),selfdens,hv)
        rho(i) = selfdens*mass(i)
    enddo
    !     Calculate SPH sum for rho:
    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        rho(i) = rho(i) + mass(j)*w(k)
        rho(j) = rho(j) + mass(i)*w(k)
    enddo
    !     Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
    if (nor_density) then
        do i=1, ntotal
            rho(i)=rho(i)/wi(i)
        enddo
    endif
    !
    end
    !
    subroutine con_density(ntotal,mass,niac,pair_i,pair_j,dwdx,vx,itype,x,rho, drhodt,AArea)
    !----------------------------------------------------------------------
    !     Subroutine to calculate the density with SPH continuiity approach.

    !     ntotal : Number of particles                                  [in]
    !     mass   : Particle masses                                      [in]
    !     niac   : Number of interaction pairs                          [in]
    !     pair_i : List of first partner of interaction pair            [in]
    !     pair_j : List of second partner of interaction pair           [in]
    !     dwdx   : derivation of Kernel for all interaction pairs       [in]
    !     vx     : Velocities of all particles                          [in]
    !     itype   : type of particles                                   [in]
    !     x      : Coordinates of all particles                         [in]
    !     rho    : Density                                              [in]
    !     drhodt : Density change rate of each particle                [out]
    use config_parameter
    implicit none
    !
    real(8) ntotal,niac,pair_i(max_interaction),pair_j(max_interaction), itype(maxn),AArea(maxn)
    real(8) mass(maxn), dwdx(dim, max_interaction),vx(dim,maxn), x(dim,maxn), rho(maxn), drhodt(maxn)
    real(8) i,j,k,d, ii
    real(8) vcc, dvx(dim),beta,da,db


    !
    do i = 1, ntotal
        drhodt(i) = 0.
    enddo
    !
    if (nnps.eq.4) then   ! GSM algorithm
        do k=1,niac
            i = pair_i(k)
            j = pair_j(k)
                vcc=0
                if (flowtype==3) then ! discharge flow
                        if (itype(j)<0) then
                            if ((-0.25<=x(1,i)<=12.75.and.-0.25<=x(2,i)<=30).or.(12.75<=x(1,i)<=14.75.and.-0.25<=x(2,i)<=1.75)) then 
                                goto 27
                            else
                                vx(1:dim,j)=vx(d:dim,i)
                            endif
                        endif
                endif
                
27              do d=1,dim
                    !if ((d==1.and.itype(j)==-2).or.(d==2.and.itype(j)==-1)) goto 100  
                    vcc = vcc + rho(i)*(vx(d,i) - vx(d,j))*dwdx(d,k)*AArea(i)/mass(i)*rho(i)
100             enddo
                
                drhodt(i) = drhodt(i) + vcc
101     enddo
        



    else ! SPH algorithm

        do k=1,niac
            i = pair_i(k)
            j = pair_j(k)

            !! Define the BC velocity-----no slip at bottom while free-slip at left
            !if (itype(i)>0.and.itype(j)<0) then
            !    if (itype(j).eq.(-2.0)) then
            !      beta=1.5
            !      db=-x(2,j)
            !      da=x(2,i)
            !      if (1+db/da<beta) beta=1+db/da
            !      vx(1:dim,j)=(1-beta)*vx(1:dim,i) ! no_slip bc
            !    endif
            !    if (itype(j).eq.(-1.0)) then
            !      vx(1,j)=-vx(1,i)
            !      vx(2,j)=vx(2,i)
            !    endif
            !endif
            if (itype(j)>0) then
                do d=1,dim
                    dvx(d) = vx(d,i) - vx(d,j)
                enddo
                vcc=0
                do d=1,dim
                    vcc = vcc + dvx(d)*dwdx(d,k)
                enddo
                drhodt(i) = drhodt(i) + mass(j)*vcc
                drhodt(j) = drhodt(j) + mass(i)*vcc
            endif
        enddo
    endif

    end