    subroutine SPH_pressure(ntotal, itype, niac, pair_i, pair_j,dwdx, p, mass, rho, vx, dvxdt,AArea,x,near)

    !      Calculate SPH sum for pressure force -p,a/rho
    !      and the internal energy change de/dt due to -p/rho vc,c
    use config_parameter
    implicit none
    !
    real(8) ntotal, itype(maxn),niac,pair_i(max_interaction),pair_j(max_interaction),x(dim,maxn)
    real(8) p(maxn), mass(maxn), rho(maxn),dwdx(dim, max_interaction), vx(dim, maxn),dvxdt(dim, maxn),AArea(maxn)
    real(8) i, j, k, d,near(maxn)
    real(8) h, rhoij
    !     Initialization

    dvxdt(1:dim, 1:maxn) = 0.e0
    near(1:ntotal)=0.

    !     SPH summation
    do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        if (nnps.eq.4) then
                do d=1,dim
                    if (flowtype==3) then ! discharge flow
                        if (itype(j)<0) then
                            if ((-0.25<=x(1,i)<=12.75.and.-0.25<=x(2,i)<=30).or.(12.75<=x(1,i)<=14.75.and.-0.25<=x(2,i)<=1.75)) then 
                                goto 27
                            else
                                p(j)=0
                                goto 27
                            endif
                        endif
                    endif
27                  h = (-1)*(p(j)+p(i))*dwdx(d,k)*AArea(i) ! minus because -p=stress
                    if (flowtype==4.and.itype(i)==100) h=h*2 ! for area consideration
                    if (0.06.ge.sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2))  then
                            if (x(d,j)>x(d,i)) then
                                h=abs(h*(0.1/(sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)+0.04))**0.1)*(-1)
                            else
                                h=abs(h*(0.1/(sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)+0.04))**0.1)
                            endif
                            near(i)=1
                            near(j)=1
                    endif                    
                    
                    !if (d==1.and.itype(j)==-2) h=0.0
                    !if (d==2.and.itype(j)==-1) h=0.0
                    dvxdt(d,i) = dvxdt(d,i) + h/mass(i)
                    !_____________________________________________________________________________________________________
                    if(isNaN(dwdx(d,k))) then
                        write(*,'(A)') 'SPH pressure:dwdx is NaN!'
                        pause
                    endif
                    if(isNaN(p(i))) then
                        write(*,*) 'SPH pressure:p is NaN!'
                        write(*,*) 'i ',i
                    elseif (isNaN(p(j))) then
                        write(*,*) 'itype(j):  ',p(j)
                        write(*,*) 'j  ',j

                        pause
                    endif
                    if(isNaN(rho(i))) then
                        write(*,'(A)') 'SPH pressure:rho is NaN!'
                        pause
                    endif
                    if(isNaN(rhoij)) then
                        write(*,'(A)') 'SPH pressure:rhoij is NaN!'
                        pause
                    endif
                    if(isNaN(h)) then
                        write(*,'(A)') 'SPH pressure:h is NaN!'
                        pause
                    endif
                enddo

        else
            if (itype(j)>0) then
                !---- For SPH algorithm 1
                if(pa_sph.eq.1) then
                    rhoij = 1.e0/(rho(i)*rho(j))
                    do d=1,dim
                        h = (p(i) + p(j))*dwdx(d,k)
                        h = h*rhoij
                        dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
                        dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
                    enddo

                    !---- For SPH algorithm 2
                else if (pa_sph.eq.2) then
                    !_____________________________________________________________________________________________________
                    do d=1,dim
                        h = -(p(i)/rho(i)**2 + p(j)/rho(j)**2)*dwdx(d,k)
                        dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
                        dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
                        !_____________________________________________________________________________________________________

                        if(isNaN(dwdx(d,k))) then
                            write(*,'(A)') 'SPH pressure:dwdx is NaN!'
                            pause
                        endif
                        if(isNaN(p(i)).or.isNaN(p(j))) then
                            write(*,'(A)') 'SPH pressure:p is NaN!'
                            write(*,'(A,I5)') 'itype(i):  ',itype(i)
                            write(*,'(A,I5)') 'itype(j):  ',itype(j)
                            pause
                        endif
                        if(isNaN(rho(i))) then
                            write(*,'(A)') 'SPH pressure:rho is NaN!'
                            pause
                        endif
                        if(isNaN(rhoij)) then
                            write(*,'(A)') 'SPH pressure:rhoij is NaN!'
                            pause
                        endif
                        if(isNaN(h)) then
                            write(*,'(A)') 'SPH pressure:h is NaN!'
                            pause
                        endif
                    enddo

                endif
            endif
        endif
124    enddo
    


    end