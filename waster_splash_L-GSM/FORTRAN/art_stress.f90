      subroutine art_stress(ntotal,hsml,mass,x,vx,niac,rho,c, pair_i,pair_j,itype,w,dwdx,stress,p,dvxdt)
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
      real(8) hsml(maxn), mass(maxn), x(dim,maxn),vx(dim,maxn),rho(maxn), c(maxn), w(max_interaction),dwdx(dim,max_interaction), dvxdt(dim,maxn), stress(maxn, dim, dim),p(maxn),totalstress(maxn, dim, dim)
      real(8) i,j,k,d,theta(ntotal),d1,d2,eps,h,ratio(ntotal),n,f,rr,dd
      real(8) dx, dvx(dim),rRx(ntotal),rRy(ntotal),rostressx(ntotal),rostressy(ntotal),R(maxn, dim, dim),itype(maxn)

      !parameters needed to be assigned

      eps=0.5
      
      n=2.55
      dd=0.04

       do i=1,ntotal
        do d=1,dim
          dvxdt(d,i) = 0.e0
        enddo
      enddo

      do i=1,ntotal
        do d1=1,dim
          do d2=1,dim
          if (d1.ne.d2)then
          totalstress(i,d1,d2)=stress(i,d1,d2)
          else
          totalstress(i,d1,d2)=stress(i,d1,d2)-p(i)
          endif
          end do
        end do
    end do





     do i=1,ntotal
      if(totalstress(i,1,1).eq.totalstress(i,2,2))then
      theta(i)= 0.
      else
      theta(i)= 0.5*atan(2*totalstress(i,1,2)/(totalstress(i,1,1)-totalstress(i,2,2)))
      endif

      !rostress-xx
      rostressx(i)=cos(theta(i))*cos(theta(i))*totalstress(i,1,1)+2*cos(theta(i))*sin(theta(i))*totalstress(i,1,2)+sin(theta(i))*sin(theta(i))*totalstress(i,2,2)
      rostressy(i)=sin(theta(i))*sin(theta(i))*totalstress(i,1,1)-2*cos(theta(i))*sin(theta(i))*totalstress(i,1,2)+cos(theta(i))*cos(theta(i))*totalstress(i,2,2)

      if(rostressx(i).gt.0)then
      rRx(i)=-eps*rostressx(i)/rho(i)**2
      else
      rRx(i)=0.
      end if

      if(rostressy(i).gt.0)then
      rRy(i)=-eps*rostressy(i)/rho(i)**2
      else
      rRy(i)=0.
      end if


     R(i,1,1)=rRx(i)*cos(theta(i))*cos(theta(i))+rRy(i)*sin(theta(i))*sin(theta(i))
     R(i,2,2)=rRx(i)*sin(theta(i))*sin(theta(i))+rRy(i)*cos(theta(i))*cos(theta(i))
     R(i,1,2)=(rRx(i)-rRy(i))*sin(theta(i))*cos(theta(i))
     R(i,2,1)=(rRx(i)-rRy(i))*sin(theta(i))*cos(theta(i))

     enddo


      do k=1,niac
        i = pair_i(k)
        j = pair_j(k) 
        h=0.

        rr=sqrt((x(1,i)-x(1,j))**2+(x(2,i)-x(2,j))**2)

        if(rr.lt.dd)then

        ratio(i)=15./(7.*pi*hsml(i)**2)*0.26157407
        f=(w(k)/ratio(i))**n
!         f=11.

         do d=1,dim                
            if (d.eq.1) then
!     x-coordinate of acceleration
              h = (R(i,1,1)+ R(j,1,1))*f*dwdx(1,k)
                    
                    if(isNaN(h)) then
                          write(*,'(A)') 'viscous force:h is NaN!!'
                     pause
                     endif

                    if (dim.ge.2) then
                           h = h + (R(i,1,2) +R(j,1,2))*f*dwdx(2,k)
!                            
                        if(isNaN(h)) then
                           write(*,'(A)') 'viscous force:h is NaN!!'
                           pause
                        endif
!                
                    endif
            endif
           
           
           if (d.eq.2) then
!     y-coordinate of acceleration
! 
               h = (R(i,2,1)+R(j,2,1))*f*dwdx(1,k)+(R(i,2,2)+R(j,2,2))*f*dwdx(2,k)
!
                  if(isNaN(h)) then
                       write(*,'(A)') 'viscous force:h is NaN!!'
                       pause
                  endif
                   

            endif

            dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
            dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
!            if(itype(i).ge.100)then
!             dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
!             else
!             dvxdt(d,i) =0.
!             endif
!
!             if(itype(j).ge.100)then
!             dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
!             else
!             dvxdt(d,j) =0.
!             endif


             enddo

             endif
             enddo

    end