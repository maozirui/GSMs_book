    subroutine direct_GSM(itimestep,hsml,ntotal,nvirt,nvirt2,x,vx,p,mass,rho,itype,niac,pair_i,pair_j,w,dwdx,AArea,numb,countiac)
    !----------------------------------------------------------------------
    !   Subroutine to calculate the smoothing funciton for each particle and
    !   the interaction parameters used by the GSM algorithm. Interaction
    !   pairs are determined by directly comparing the particle distance.

    !     itimestep : Current time step                                 [in]
    !     ntotal    : Number of particles                               [in]
    !     hsml      : Smoothing Length                                  [in]
    !     x         : Coordinates of all particles                      [in]
    !     mass      : mass of all particles                             [in]
    !     rho       : density of all particles                          [in]
    !     niac      : Number of interaction pairs                      [out]
    !     pair_i    : List of first partner of interaction pair        [out]
    !     pair_j    : List of second partner of interaction pair       [out]
    !     w         : Kernel for all interaction pairs                 [out]
    !     dwdx      : Derivative of kernel with respect to x, y and z  [out]
    use config_parameter
    implicit none
    character (40) :: xname
    character (40) :: yname
    real(8) ntotal,niac,pair_i(max_interaction),pair_j(max_interaction),nvirt,itype(maxn)
    real(8) mass(maxn), rho(maxn), x(dim,maxn), dwdx(dim,max_interaction),w(max_interaction),XX(dim,maxn)
    real(8) distj,anglej,dx,dy,C1,C2,C3,C4,C,distt,a,nc,sig(maxn,15)
    real(8) uangle,dangle,sigma1,sigma2,sigma3,sigma4,nvirt2
    real(8) dxiac(dim), driac, r, mhsml, tdwdx(dim),sumiac, maxiac, miniac, noiac,maxp, minp,origin_L,origin_N
    real(8) angle(maxn,15),upper(dim),lower(dim),temp_dist(15),dis,dis_cr,hsml(maxn)
    real(8) dsx(15),dsy(15),dxiac1,min_dis,vx(dim,maxn),p(maxn),area(dim,maxn),AArea(maxn),x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6,Last,Next
    integer(4) jp1,jm1,jj,jjj,numbj,k,l,ll,i, j, d,itimestep,nabes_dim,nabes_max
    integer(4) countiac(maxn),nabes(max_interaction),nabes_first(maxn),ntotall,numb(maxn,15),temp_numb(15)
    integer(4),save:: countiac_save(maxn),numb_save(maxn,15)
    real(8),save:: itype_save(maxn),sig_save(maxn,15)
    !
    !      do i=1,max_interaction
    !         do d=1,dim
    !            dwdx(d,i)=0.0
    
    !         enddo
    !         w(i)=0.0
    !      enddo
    !
    
    !if (itimestep==1) then
    !    goto 44 ! remeshing
    !else
    !    goto 55 ! recheck
    !endif
    !!!!!!!!!!!!!!!!!!!!!!!!!!!! option 1: Triangulation scheme!!!!!!!!!!!!!!!!!!!!!!!!!!
44    ntotall=ntotal+nvirt
    if (mod(itimestep,10).eq.1) then
      call triangulation_order3 (x, itimestep,ntotall, itype,hsml(1), countiac, numb)
      countiac_save=countiac
      numb_save=numb
    else
        countiac=countiac_save
         numb=numb_save
      endif 
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! option 2: direct find scheme !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !      do i=1,ntotal
    !         countiac(i) = 0
    !         do j=1,15
    !            numb(i,j)=0.0
    !         enddo
    !      enddo
    !
    !      do i=1,ntotal-1
    !        do j = i+1, ntotal
    !          mhsml = (hsml(i)+hsml(j))/2.
    !          dxiac(1) = x(1,i) - x(1,j)
    !
    !          if (abs(dxiac(1)).ge.(1.4*mhsml))  go to 1       !!!!!!!!!!!! 1.4*1.1=1.54 dx
    !          driac    = dxiac(1)*dxiac(1)
    !          do d=2,dim
    !            dxiac(d) = x(d,i) - x(d,j)
    !            if (abs(dxiac(d)).ge.(1.4*mhsml))  go to 1  !!!!!!!!!!!! 1.4*1.1=1.54 dx
    !            driac    = driac + dxiac(d)*dxiac(d)
    !          enddo
    !          if (sqrt(driac).lt.1.4*mhsml) then           !!!!!!!!!!!! 1.4*1.1=1.54 dx
    !               countiac(i) = countiac(i) + 1
    !               countiac(j) = countiac(j) + 1
    !               numb(i,countiac(i))=j
    !               numb(j,countiac(j))=i
    !          endif
    !1       enddo
    !      enddo
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !! determine the order of the 6 nearest particles based on their position
    !! relationship according to the interested particle i; and then set the
    !! interaction pairs based on the nearest particles and their anlges according
    !! to each interested particle i.
    min_dis=1000000
    do i=1,ntotal
        do j=1,countiac(i)
            k=numb(i,j)
            dx=x(1,k)-x(1,i)
            dy=x(2,k)-x(2,i)
            angle(i,j)=atan2(dy,dx)
            if (angle(i,j).lt.0) angle(i,j)=angle(i,j)+pi*2.0
        enddo
    enddo
    !

    !
    do i=1,ntotal
100     k=0
        do j=1,countiac(i)-1
            if (angle(i,j).gt.angle(i,j+1)) then
                anglej=angle(i,j+1)
                angle(i,j+1)=angle(i,j)
                angle(i,j)=anglej
                numbj=numb(i,j+1)
                numb(i,j+1)=numb(i,j)
                numb(i,j)=numbj
                k=k+1
            endif
        enddo
        if (k.gt.0) go to 100
    enddo

!    ! generate ghost particles based on the maximum angle of each real particle
    !!
!    nvirt2=0.0
!    do i=1,ntotal
!        if (countiac(i)>1) then
!        do j=1,countiac(i)
!            k=numb(i,j)
!            jj=j-1
!            if (jj.lt.1) jj=countiac(i)
!            jjj=j+1
!            if (jjj.gt.countiac(i)) jjj=1
!            jp1=numb(i,jjj)
!            jm1=numb(i,jj)
!            uangle=angle(i,jjj)-angle(i,j)
!            if (uangle.lt.0) uangle=uangle+2*pi
!
!            if (uangle.ge.pi) then
!                nc=2
!                do l=1,nc
!                    ll=numb(i,l)
!                    nvirt2=nvirt2+1
!                    countiac(i)=countiac(i)+1
!                    numb(i,countiac(i))=ntotal+nvirt+nvirt2
!                    angle(i,countiac(i))=angle(i,jjj)-uangle/(nc+1)*(nc-l+1)
!                    if (angle(i,countiac(i))<0.0) angle(i,countiac(i))=angle(i,countiac(i))+2*pi
!                    x(1,ntotal+nvirt+nvirt2)=x(1,i)+0.1*cos(angle(i,countiac(i)))
!                    x(2,ntotal+nvirt+nvirt2)=x(2,i)+0.1*sin(angle(i,countiac(i)))
!                    vx(1,ntotal+nvirt+nvirt2)=vx(1,i)
!                    vx(2,ntotal+nvirt+nvirt2)=vx(2,i)
!                    mass(ntotal+nvirt+nvirt2)=mass(i)
!                    rho(ntotal+nvirt+nvirt2)=rho(i)
!                    p(ntotal+nvirt+nvirt2)=0
!                    itype(ntotal+nvirt+nvirt2)=-200.0
!                    itype(i)=100
!146                enddo
!                goto 301
!            endif
!        enddo
!        itype(i)=200
!        endif
!301 enddo
!
!    do i=1,ntotal
!401     k=0
!        do j=1,countiac(i)-1
!            if (angle(i,j).gt.angle(i,j+1)) then
!                anglej=angle(i,j+1)
!                angle(i,j+1)=angle(i,j)
!                angle(i,j)=anglej
!                numbj=numb(i,j+1)
!                numb(i,j+1)=numb(i,j)
!                numb(i,j)=numbj
!                k=k+1
!            endif
!        enddo
!        if (k.gt.0) go to 401
!    enddo
      !write(yname, '(A)') 'result_after.dat'
      !open (2, file = yname)
      !do i=1,ntotal
      !      write(2,'(11(I4,2x))') i, numb(i,1), numb(i,2), numb(i,3), numb(i,4), numb(i,5),  numb(i,6),  numb(i,7),  numb(i,8),  numb(i,9),  numb(i,10)
      !enddo
      !close(2)

!!     generate ghost particles based on the maximum angle of each real particle
!    !!
!    nvirt2=0.0
!    do i=1,ntotal
!        do j=1,countiac(i)
!            k=numb(i,j)
!            jj=j-1
!            if (jj.lt.1) jj=countiac(i)
!            jjj=j+1
!            if (jjj.gt.countiac(i)) jjj=1
!            jp1=numb(i,jjj)
!            jm1=numb(i,jj)
!            uangle=angle(i,jjj)-angle(i,j)
!            if (uangle.lt.0) uangle=uangle+2*pi
!
!            if (uangle.ge.pi) then
!                itype(i)=100
!                goto 301
!            endif
!        enddo
!        itype(i)=200
!301 enddo
    
    ! search for the free surface particles in a constitutive way
    itype(1:ntotal)=200
    sig(1:ntotal,1:15)=0.0
    do i=1,ntotal
        if (itype(i).eq.200) then
            do j=1,countiac(i)
                k=numb(i,j)
                jjj=j+1
                if (jjj.gt.countiac(i)) jjj=1
                jp1=numb(i,jjj)
                uangle=angle(i,jjj)-angle(i,j)
                if (uangle.lt.0) uangle=uangle+2*pi
                
                if (uangle.ge.0.99*pi) then
                    itype(i)=100 ! surface particle
                    sig(i,j)=1 ! =1 surface edge; =0 not surface edge  
                    sig(i,jjj)=1 ! the opposite direction should also be denoted as surface edge
                    if (itype(k)>0) then
                        itype(k)=100
                        goto 221  ! searching right chain
                    else
                        goto 321 ! searching left chain
                    endif
                    
                endif
            enddo
            goto 242
            
221         Next=k
            Last=i
223         do j=1,countiac(Next)
                k=numb(Next,j)
                if (k.eq.Last) then
                    sig(Next,j)=1
                    jj=j-1
                    if (jj.lt.1) jj=countiac(Next)
                    jm1=numb(Next,jj)
                    sig(Next,jj)=1 ! surface edge
                    if (itype(jm1)>0) then
                        itype(jm1)=100
                    else
                        goto 321
                    endif
                    Last=Next
                    Next=jm1
                    if (Next.eq.i) then
                        goto 242
                    else
                        goto 223
                    endif
                endif
            enddo
            write (*,*) 'BBBBBBBBBBBBBBBBBBBBBBB',Next
            pause 
            
321         Next=jp1
            Last=i
            if (itype(Next)>0) then
                itype(Next)=100
            else
                goto 242
            endif
323         do j=1,countiac(Next)
                k=numb(Next,j)
                if (k.eq.Last) then
                    sig(Next,j)=1
                    jjj=j+1
                    if (jjj>countiac(Next)) jjj=1
                    jp1=numb(Next,jjj)
                    sig(Next,jjj)=1
                    if (itype(jp1)>0) then
                        itype(jp1)=100
                    else
                        goto 242
                    endif
                    Last=Next
                    Next=jp1
                    if (Next==i) then
                        goto 242
                    else 
                        goto 323
                    endif
                endif
            enddo
            write (*,*) 'BBBBBBBBBBBBBBBBBBBBBBB',Next
            pause
                    
        endif
242    enddo
        itype_save=itype
        sig_save=sig
        goto 334            
                    
!!! Recheck the mesh!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!55  min_dis=1000000
!    countiac=countiac_save
!    numb=numb_save
!    itype=itype_save
!    sig=sig_save
!    do i=1,ntotal
!        do j=1,countiac(i)
!            k=numb(i,j)
!            dx=x(1,k)-x(1,i)
!            dy=x(2,k)-x(2,i)
!            angle(i,j)=atan2(dy,dx)
!            if (angle(i,j).lt.0) angle(i,j)=angle(i,j)+pi*2.0
!        enddo
!    enddo
!    !
!    do i=1,ntotal
!        if (countiac(i)>2) then
!            do j=2,countiac(i)
!                dangle=(angle(i,j+1)-angle(i,1))*(angle(i,j)-angle(i,1))
!                if (dangle<0) then
!                    goto 44 ! remeshing
!                else
!                    uangle=(angle(i,j+1)-angle(i,1))
!                    dangle=(angle(i,j)-angle(i,1))
!                    if (uangle<dangle) goto 44
!                endif
!            enddo
!        endif
!    enddo
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                

334 niac=0.0
    
    do i=1,ntotal
        AArea(i)=0.0
        sigma1=0.0
        sigma2=0.0
        sigma3=0.0
        sigma4=0.0
        if (countiac(i)<2) then
            itype(i)=50
            goto 60
        endif
        
        do j=1,countiac(i)
            k=numb(i,j)
            jj=j-1
            if (jj.lt.1) jj=countiac(i)
            jjj=j+1
            if (jjj.gt.countiac(i)) jjj=1
            jp1=numb(i,jjj)
            jm1=numb(i,jj)
            uangle=angle(i,jjj)-angle(i,j)
            if (uangle.lt.0) uangle=uangle+2*pi
            dangle=angle(i,j)-angle(i,jj)
            if (dangle.lt.0) dangle=dangle+2*pi
                
            !if (uangle.lt.pi) then
            !    do d=1,dim
            !        upper(d)=(x(d,i)+x(d,k)+x(d,jp1))/3.0
            !    enddo
            !else
            !    do d=1,dim
            !        upper(d)= (x(d,i)+x(d,k))/2.0
            !    enddo
            !endif
            !
            !if (dangle.lt.pi) then
            !    do d=1,dim
            !        lower(d)=(x(d,i)+x(d,k)+x(d,jm1))/3.0
            !    enddo
            !else
            !    do d=1,dim
            !        lower(d)=(x(d,i)+x(d,k))/2.0
            !    enddo
            !endif
            if ((sig(i,j).eq.1.and.sig(i,jjj).eq.1.and.sig(i,jj).ne.1.0).or.(sig(i,j).eq.1.and.sig(i,jjj).eq.1.and.sig(i,jj).eq.1.0.and.uangle.gt.pi)) then
                do d=1,dim
                    upper(d)= (x(d,i)+x(d,k))/2.0
                enddo
            else
                do d=1,dim
                    upper(d)=(x(d,i)+x(d,k)+x(d,jp1))/3.0
                enddo
            endif
            
            if ((sig(i,j).eq.1.and.sig(i,jj).eq.1.and.sig(i,jjj).ne.1.0).or.(sig(i,j).eq.1.and.sig(i,jj).eq.1.and.sig(i,jjj).eq.1.0.and.dangle.gt.pi)) then
                do d=1,dim
                    lower(d)=(x(d,i)+x(d,k))/2.0
                enddo
            else
                do d=1,dim
                    lower(d)=(x(d,i)+x(d,k)+x(d,jm1))/3.0
                enddo
            endif
            dsx(j)=upper(2)-lower(2)
            dsy(j)=upper(1)-lower(1)
            
            x1=x(1,i)
            y1=x(2,i)
            x2=(x(1,i)+x(1,k))/2.0
            y2=(x(2,i)+x(2,k))/2.0
            x3=upper(1)
            y3=upper(2)
            x4=x(1,i)
            y4=x(2,i)
            x5=(x(1,i)+x(1,k))/2.0
            y5=(x(2,i)+x(2,k))/2.0
            x6=lower(1)
            y6=lower(2)
            AArea(i)=AArea(i)+abs(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)/2.0+abs(x4*y5-x5*y4+x5*y6-x6*y5+x6*y4-x4*y6)/2.0
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!option 1: area!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !  x1=x(1,i)
          !  y1=x(2,i)
          !  x2=(x(1,i)+x(1,k))/2.0
          !  y2=(x(2,i)+x(2,k))/2.0
          !  x3=upper(1)
          !  y3=upper(2)
          !  x4=x(1,i)
          !  y4=x(2,i)
          !  x5=(x(1,i)+x(1,k))/2.0
          !  y5=(x(2,i)+x(2,k))/2.0
          !  x6=lower(1)
          !  y6=lower(2)
          !  AArea(i)=AArea(i)+abs(x1*y2-x2*y1+x2*y3-x3*y2+x3*y1-x1*y3)/2.0+abs(x4*y5-x5*y4+x5*y6-x6*y5+x6*y4-x4*y6)/2.0
          !  
          !enddo
          ! if (AArea(i).eq.0) write(*,*) '*************',i
          ! if(isNaN(area(1,i))) then
          !       write(*,*) 'direct_GSM:area(1,i) is NaN!!',i
          !       pause
          !   endif
          !
          !   if(isNaN(area(2,i))) then
          !       write(*,*) 'direct_GSM:area(2,i) is NaN!!',i
          !       pause
          !   endif
          !
          !do j=1,countiac(i)
          !   niac=niac+1
          !   k=numb(i,j)
          !   pair_i(niac)=i
          !   pair_j(niac)=k
          !   dwdx(1,niac)=dsx(j)/AArea(i)/2.0
          !   dwdx(2,niac)=dsy(j)/AArea(i)/2.0
          !   if(isNaN(dwdx(1,niac))) then
          !       write(*,*) 'direct_GSM:dwdx1 is NaN!!',i
          !       pause
          !   endif
          !
          !   if(isNaN(dwdx(2,niac))) then
          !       write(*,*) 'direct_GSM:dwdx2 is NaN!!',i
          !       pause
          !   endif
          !
          !enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!option 2: advanced !!!!!!!!!!!!!!!!!!!!!!!!!!
            sigma1=sigma1+dsx(j)*(x(1,k)-x(1,i))/2.0
            sigma2=sigma2+dsx(j)*(x(2,k)-x(2,i))/2.0
            sigma3=sigma3+dsy(j)*(x(1,k)-x(1,i))/2.0
            sigma4=sigma4+dsy(j)*(x(2,k)-x(2,i))/2.0
        enddo
        
        C=sigma1*sigma4-sigma2*sigma3
        if(C.eq.0.0) then
            write(*,*) 'direct_GSM:1/Cccc is NaN!!',i
            goto 60
        endif
        
        C1=sigma4/C
        C2=(-1.0)*sigma2/C
        C3=(-1.0)*sigma3/C
        C4=sigma1/C
        !          write (*,*) 'countiac(i)=',countiac(i)
        do j=1,countiac(i)
            niac=niac+1
            k=numb(i,j)
            pair_i(niac)=i
            pair_j(niac)=k
            dwdx(1,niac)=(dsx(j)*C1+dsy(j)*C2)/2.0
            dwdx(2,niac)=(dsx(j)*C3+dsy(j)*C4)/2.0
            !if (itype(i).eq.100) then
            !    dwdx(1,niac)=dwdx(1,niac)/4.0
            !    dwdx(2,niac)=dwdx(2,niac)/4.0
            !endif
            
            if(isNaN(dwdx(1,niac))) then
                write(*,*) 'direct_GSM:dwdx1 is NaN!!',i
                pause
            endif
        
            if(isNaN(dwdx(2,niac))) then
                write(*,*) 'direct_GSM:dwdx2 is NaN!!',i
                pause
            endif
        
        enddo
60  enddo
     



    !
    !     if (mod(itimestep,print_step).eq.0) then
    !     !     Statistics for the interaction
    !     sumiac = 0
    !     maxiac = 0
    !     miniac = 1000
    !     noiac  = 0
    !     do i=1,ntotal
    !       sumiac = sumiac + countiac(i)
    !       if (countiac(i).gt.maxiac) then
    !  maxiac = countiac(i)
    !  maxp = i
    !endif
    !if (countiac(i).lt.miniac) then
    !  miniac = countiac(i)
    !         minp = i
    !endif
    !       if (countiac(i).eq.0)      noiac  = noiac + 1
    !     enddo
    !       if (int_stat) then
    !         print *,' >> Statistics: interactions per particle:'
    !         print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
    !         print *,'**** Particle:',minp, ' minimal interactions:',miniac
    !         print *,'Average :',real(sumiac)/real(ntotal)
    !         print *,'Total pairs : ',niac
    !         print *,'**** Particles with no interactions:',noiac
    !       endif
    !     endif

    

    end

