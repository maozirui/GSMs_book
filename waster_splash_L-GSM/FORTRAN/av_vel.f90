    subroutine av_vel(ntotal,mass,itype,niac,pair_i,pair_j,w,dwdx, x,vx, rho, av,countiac,near)
    !----------------------------------------------------------------------
    !     Subroutine to calculate the average velocity to correct velocity
    !     for preventing penetration (monaghan, 1992)
    !     This is XSPH tech. as Eq. (4.17) in SPH Book (Liu G.R. and Liu M.B.)
    !------------------------------------------------------------------------
    !     ntotal : Number of particles                                  [in]
    !     mass   : Particle masses                                      [in]
    !     niac   : Number of interaction pairs                          [in]
    !     pair_i : List of first partner of interaction pair            [in]
    !     pair_j : List of second partner of interaction pair           [in]
    !     w      : Kernel for all interaction pairs                     [in]
    !     vx     : Velocity of each particle                            [in]
    !     rho    : Density of each particle                             [in]
    !     av     : Average velocity of each particle                    [out]
    use config_parameter
    implicit none
    !
    real(8) ntotal, niac, pair_i(max_interaction),av_rho(maxn),drhodx(dim,maxn)
    real(8) pair_j(max_interaction)
    real(8) mass(maxn),w(max_interaction),dwdx(dim,max_interaction)
    real(8) vx(dim,maxn), rho(maxn), av(dim, maxn),itype(maxn),x(dim,maxn)
    real(8) i,j,k,d
    real(8)   vcc, dvx(dim), epsilon,dvdx(dim,maxn),near(maxn)
    integer(4) countiac(maxn),kn(maxn)
    !     epsilon --- a small constants chosen by experence, may lead to instability.
    !     for example, for the 1 dimensional shock tube problem, the E <= 0.3
    epsilon = 0.02
    !
    do i = 1, ntotal
        do d = 1, dim
            av(d,i) = 0.
            dvdx(d,i)=0.
            drhodx(d,i)=0.
        enddo
        kn(i)=0
        av_rho(i)=0.
    enddo
    !
    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)

        do d=1,dim
            if ((d==1.and.itype(j)==-2).or.(d==2.and.itype(j)==-1)) goto 100
            dvdx(d,i) = dvdx(d,i) - (vx(d,i) - vx(d,j))*dwdx(d,k)
  !          drhodx(d,i) = drhodx(d,i) + (rho(j)-rho(i))*dwdx(d,k)
100     enddo
        
    enddo
    !
    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        if (near(i)==0.and.itype(j)>0) then
            do d=1,dim
                av(d, i) = av(d,i) + (vx(d,i)+vx(d,j))/2 + (x(d,i)-x(d,j))/2*dvdx(d,i)
   !           av_rho(i) = av_rho(i) + (rho(i)+rho(j))/2 + (x(d,i)-x(d,j))/2*drhodx(d,i)
            enddo
            kn(i) = kn(i) + 1
        endif
    enddo
    !
    
    do i = 1, ntotal
        if (kn (i) .ne. 0) then
            do d = 1, dim
                av(d,i) = epsilon * (av(d,i)/kn(i)-vx(d,i))
            enddo
  !          av_rho(i) = epsilon*(av_rho(i)/kn(i)/dim-rho(i))
        endif
    enddo
 
    !
    end