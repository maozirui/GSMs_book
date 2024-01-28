    subroutine direct_find(itimestep, ntotal,hsml,x,niac,pair_i,pair_j,w,dwdx,countiac,min_dis)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the smoothing funciton for each particle and
    !   the interaction parameters used by the SPH algorithm. Interaction
    !   pairs are determined by directly comparing the particle distance
    !   with the corresponding smoothing length.

    !     itimestep : Current time step                                 [in]
    !     ntotal    : Number of particles                               [in]
    !     hsml      : Smoothing Length                                  [in]
    !     x         : Coordinates of all particles                      [in]
    !     niac      : Number of interaction pairs                      [out]
    !     pair_i    : List of first partner of interaction pair        [out]
    !     pair_j    : List of second partner of interaction pair       [out]
    !     w         : Kernel for all interaction pairs                 [out]
    !     dwdx      : Derivative of kernel with respect to x, y and z  [out]
    !     countiac  : Number of neighboring particles                  [out]
    use config_parameter
    implicit none
    !
    real(8)  ntotal,niac,pair_i(max_interaction),pair_j(max_interaction), countiac(maxn)
    real(8) hsml(maxn), x(dim,maxn), w(max_interaction),dwdx(dim,max_interaction)
    real(8) i, j, d,  sumiac, maxiac, miniac, noiac,maxp, minp, scale_k
    real(8) dxiac(dim), driac, r, mhsml, tdwdx(dim),min_dis
    integer(4) itimestep
    !
    if (skf.eq.1) then
        scale_k = 2
    else if (skf.eq.2) then
        scale_k = 3
    else if (skf.eq.3) then
        scale_k = 3
    else if (skf.eq.4) then
        scale_k = 2
    endif
    !
    do i=1,ntotal
        countiac(i) = 0
    enddo
    !
    min_dis=100000
    do i=1,ntotal-1
        do j = i+1, ntotal
            mhsml = (hsml(i)+hsml(j))/2.
            dxiac(1) = x(1,i) - x(1,j)
            if (abs(dxiac(1)).ge.(scale_k*mhsml))  go to 1
            driac    = dxiac(1)*dxiac(1)
            do d=2,dim
                dxiac(d) = x(d,i) - x(d,j)
                if (abs(dxiac(d)).ge.(scale_k*mhsml))  go to 1
                driac    = driac + dxiac(d)*dxiac(d)
            enddo
            if (sqrt(driac)<min_dis) min_dis=sqrt(driac)
            
            if (sqrt(driac).lt.scale_k*mhsml) then
                if (niac.lt.max_interaction) then
                    !     Neighboring pair list, and totalinteraction number and
                    !     the interaction number for each particle
                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(driac)
                    !              if(r.le.5.) then
                    !                 write(*,'(A,2x,f20.12)') 'The two particles are too close!! r=',r
                    !!                 pause
                    !              endif
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1
                    !     Kernel and derivations of kernel
                    call kernel(r,dxiac,mhsml,w(niac),tdwdx)
                    do d=1,dim
                        dwdx(d,niac) = tdwdx(d)
                    enddo
                else
                    write(*,*) ' >>> ERROR <<< : Too many interactions'
                    pause
                endif
            endif
1       enddo
    enddo
    !     Statistics for the interaction
    sumiac = 0
    maxiac = 0
    miniac = 1000
    noiac  = 0
    do i=1,ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i).gt.maxiac) then
            maxiac = countiac(i)
            maxp = i
        endif
        if (countiac(i).lt.miniac) then
            miniac = countiac(i)
            minp = i
        endif
        if (countiac(i).eq.0)      noiac  = noiac + 1
    enddo
    !
    if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
            !          print *,' >> Statistics: interactions per particle:'
            !          print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
            !          print *,'**** Particle:',minp, ' minimal interactions:',miniac
            print *,'Average :',real(sumiac)/real(ntotal)
            print *,'Total pairs : ',niac
            !          print *,'**** Particles with no interactions:',noiac
        endif
    endif
    !
    end