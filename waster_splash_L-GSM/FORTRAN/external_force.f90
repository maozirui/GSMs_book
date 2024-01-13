    subroutine ext_force(itimestep,ntotal,nvirt1,nvirt2,mass,x,niac,pair_i,pair_j,itype,hsml,dvxdt,v)
    !--------------------------------------------------------------------------
    !     Subroutine to calculate external forces, e.g. GRAVITY or BODY force
    !     The forces from the interactions with BOUNDARY VIRTUAL PARTICLES
    !     are also calculated here as external forces.

    !     ntotal  : Number of particles                                 [in]
    !     mass    : Particle masses                                     [in]
    !     x       : Coordinates of all particles                        [in]
    !     pair_i : List of first partner of interaction pair            [in]
    !     pair_j : List of second partner of interaction pair           [in]
    !     itype   : type of particles                                   [in]
    !     hsml   : Smoothing Length                                     [in]
    !     Pr: Poissoin's ratio                                          [in]
    !     E: Young's modulus
    !     dvxdt   : Acceleration with respect to x, y and z            [out]
    use config_parameter
    implicit none
    !
    real(8) ntotal, nvirt1,nvirt2,itype(maxn), niac,pair_i(max_interaction), pair_j(max_interaction)
    real(8) mass(maxn), x(dim,maxn), hsml(maxn),dvxdt(dim,maxn),v(dim,maxn)
    real(8) i, j, k, d, kn, Pr, E, Fn, dn, mu, Ve, delta1, rr1,delta2,Xe,deltax,deltay
    real(8) dx(dim), rr, f, rr0, dd, n1, n2,p1,p2,rx,ry,r,c,cx,cy
    integer(4) itimestep
    !
    do i = 1, ntotal+nvirt1+nvirt2
        do d = 1, dim
            dvxdt(d, i) = 0.
        enddo
    enddo
    !     Consider self-gravity
    if (self_gravity) then
        do i = 1, ntotal
            dvxdt(2, i) =-9.8
        enddo
    endif
    !!     Boundary particle force and penalty anti-penetration force.
    !if (hvisphere) then
    !    rr0 = 6.e-004
    !    dd = 4.e+007
    !    n1 = 12
    !    n2 = 4
    !else if (tungstencube) then
    !    rr0 = 0.5
    !    dd = 0.5
    !    n1 = 6
    !    n2 = 4
    !else if (taylor) then
    !    rr0 = 5.e-4
    !    dd = 1.
    !    n1 = 6
    !    n2 = 4
    !else if (TNT) then
    !    rr0 = 1.e-3
    !    dd = 1.e5
    !    n1 = 6
    !    n2 = 4
    !endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!calculate the electric field!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !     if (itimestep.le.50000) then
    !     goto 1000
    !
    !     else
    !       c=2.5e-5
    !       cx=0.3*c   !0.25
    !       cy=1.5*c   !1.0
    !       do i=1, ntotal
    !         do j=1, ntotal
    !           if (i.ne.j) then
    !             r=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)
    !             rx=x(1,i)-x(1,j)
    !             ry=x(2,i)-x(2,j)
    !             if (r.le.0.02) then
    !               r=0.02
    !               rx=0.02*rx/sqrt(rx*rx+ry*ry)
    !               ry=0.02*ry/sqrt(rx*rx+ry*ry)
    !             endif
    !
    !             if (x(1,i).ge.0.1) then
    !               dvxdt(1,i)=dvxdt(1,i)-4.*(10.0*x(1,i)-2.0)*cx*rx/r**3
    !             else if (x(1,i).ge.0.03) then
    !                dvxdt(1,i)=dvxdt(1,i)+2.25*(5.62-54.0*x(1,i))*cx*rx/r**3
    !             else if (x(1,i).ge.-0.03) then
    !                dvxdt(1,i)=dvxdt(1,i)+3.0*cx*rx/r**3
    !             else if (x(1,i).ge.-0.08) then
    !                dvxdt(1,i)=dvxdt(1,i)+0.5*(80.0*x(1,i)+6.4)*cx*rx/r**3
    !             else
    !                dvxdt(1,i)=dvxdt(1,i)
    !             endif
    !
    !             if (x(1,i).ge.0.07) then
    !                dvxdt(2,i)=dvxdt(2,i)
    !             else if (x(1,i).ge.0) then
    !                dvxdt(2,i)=dvxdt(2,i)+11.5/(80.0*x(1,i)+1.0)*cy*ry/r**3
    !             else if (x(1,i).ge.-0.025) then
    !                dvxdt(2,i)=dvxdt(2,i)+0.30*(40.0*x(1,i)+1.0)*cy*ry/r**3
    !             else
    !                dvxdt(2,i)=dvxdt(2,i)
    !             endif
    !
    !           else
    !             dvxdt(1,i)=dvxdt(1,i)
    !             dvxdt(2,i)=dvxdt(2,i)
    !           endif
    !         enddo
    !
    !       enddo
    !     endif

    !    Boundary particle force and penalty anti-penetration force.
1000 rr0 = 0.1
    dd = 5*9.8*25
    n1 = 4
    n2 = 2
    !
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Boundary without friction .see Liu's book Page 141
    !

    if (nnps.eq.4) then !GSM case
        do i=1,ntotal
            do j=ntotal+1,ntotal+nvirt1
                dx(1)=x(1,i)-x(1,j)
                dx(2)=x(2,i)-x(2,j)
                rr=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)
                if (rr.lt.rr0) then
                    f = ((rr0/rr)**n1-(rr0/rr)**n2)/rr**2
                    if (itype(j).eq.-1.0) then
                        dvxdt(1,i)=dvxdt(1,i)+dd*dx(1)*f
                    elseif (itype(j).eq.-2.0) then
                        dvxdt(2,i)=dvxdt(2,i)+dd*dx(2)*f
                    elseif (itype(j).eq.-3.0) then
                        do d = 1, dim
                            dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
                        enddo
                    endif
                endif
            enddo
        enddo                    

    else ! SPH case
        do  k=1,niac
            i = pair_i(k)
            j = pair_j(k)

            if(itype(i)>0.and.itype(j)<0) then
                rr = 0.
                do d=1,dim
                    dx(d) =  x(d,i) -  x(d,j)
                    rr = rr + dx(d)*dx(d)
                enddo
                rr = sqrt(rr)
                rr0= 0.5
                if(rr.lt.rr0) then
                    f = ((rr0/rr)**n1-(rr0/rr)**n2)/rr**2
                    if (itype(j).eq.-1.0) then
                        dvxdt(1,i)=dvxdt(1,i)+dd*dx(1)*f
                    elseif (itype(j).eq.-2.0) then
                        dvxdt(2,i)=dvxdt(2,i)+dd*dx(2)*f
                    elseif (itype(j).eq.-3.0) then
                        do d = 1, dim
                            dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
                        enddo
                    endif

                endif
            endif
        enddo
    endif
    !
    !!!!!!!!!!!!!!!!!By HuMan Considering the boundary friction
    !    Pr=0.3    !Poissoin's ratio
    !    E=0.84e+006  !Young's modulus
    !    mu=0.3
    !       do  k=1,niac
    !        i = pair_i(k)
    !        j = pair_j(k)
    !
    !         if(itype(i).ne.itype(j)) then
    !          rr = 0.
    !
    !          do d=1,dim
    !            dx(d) =  x(d,i) -  x(d,j)
    !            rr = rr + dx(d)*dx(d)
    !          enddo
    !          rr = sqrt(rr)
    !
    !          if(rr.eq.0)then
    !          deltax=0.
    !          deltay=0.
    !          else
    !          deltax=dx(1)/rr
    !          deltay=dx(2)/rr
    !          endif
    !
    !
    !!          rr0= 0.5*(hsml(i)+hsml(j))
    !!
    !!          if(rr.lt.rr0) then
    !!          kn=2*sqrt(hsml(i)*hsml(j)/(hsml(i)+hsml(j)))/((1-Pr**2)/E)/3
    !!          dn=hsml(i)+hsml(j)-rr
    !!          Fn= kn*(dn**1.5)
    !!
    !!           f = ((rr0/rr)**n1-(rr0/rr)**n2)/rr**2
    !!            do d = 1, dim
    !!            dvxdt(d, i) = dvxdt(d, i) + dd*dx(d)*f
    !!             dvxdt(1, i) = 0.
    !!             dvxdt(2, i) = 9.8
    !!             v(1, i) = 0.
    !!             v(2, i) = 0.
    !!            enddo
    !!            endif
    !
    !            rr1= 0.5*(hsml(i)+hsml(j))
    !
    !           if (rr.lt.rr1)then
    !          kn=2*sqrt(hsml(i)*hsml(j)/(hsml(i)+hsml(j)))/((1-Pr**2)/E)/3
    !          dn=0.5*(hsml(i)+hsml(j))-rr
    !          Fn= kn*(dn**1.5)
    !          delta1=v(1,i)*v(1,i)+v(2,i)*v(2,i)
    !
    !          !!!!!!!!!!!!!1delta2=x(1,i)*x(1,i)+x(2,i)*x(2,i)
    !
    !          if(delta1.gt.0)then
    !          Ve=1/sqrt(delta1)
    !              do d = 1, dim
    !             dvxdt(1, i) = dvxdt(1, i) - mu*Ve*v(1,i)*Fn/mass(i)
    !              dvxdt(2, i) = dvxdt(2, i) - 0.0005*mu*Ve*v(2,i)*Fn/mass(i)
    !             enddo
    !          endif
    !          endif
    !        endif
    !        enddo


    !
    end