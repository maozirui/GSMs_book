      subroutine viscous_forces(ntotal, itype, niac, pair_i,pair_j, dwdx, mass, rho, vx, stress, dvxdt) 
!     Calculate SPH sum for viscous force (eta Tab),b/rho
!     stress : viscous shear stress                  [in]
!     dvxdt  : acceleration dv/dt                   [out]
      use config_parameter
      implicit none
!
      real(8) ntotal, itype(maxn),niac,pair_i(max_interaction),pair_j(max_interaction)
      real(8) mass(maxn), rho(maxn),vx(dim,maxn), dwdx(dim, max_interaction),stress(maxn, dim, dim), dvxdt(dim, maxn)       
     real(8) i, j, k, d
	  real(8) h, he, rhoij
!     Initialization 
      do i = 1, ntotal
         do d = 1, dim
            dvxdt(d, i) = 0.e0
	     enddo
	  enddo 
!
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)  
        he = 0.e0
        ! Define the artificial stress
        if (itype(i)>0.and.itype(j)<0) then
            if (itype(j).eq.(-2.0))    then
               stress(j,1:dim,1:dim)=stress(i,1:dim,1:dim)
            elseif (itype(j).eq.(-1.0)) then
               stress(j,1,1)=stress(i,1,1)
               stress(j,2,2)=stress(i,2,2)
               stress(j,1,2)=-stress(i,1,2)
               stress(j,2,1)=-stress(i,2,1)
            endif
        endif

!----     For SPH algorithm 1
        if(pa_sph.eq.1) then  
          rhoij = 1.e0/(rho(i)*rho(j))
          do d=1,dim    
             if (d.eq.1) then
!     x-coordinate of acceleration
                h = (stress(i,1,1) + stress(j,1,1))*dwdx(1,k)
                if (dim.ge.2) then
                   h = h + (stress(i,1,2)+ stress(j,1,2))*dwdx(2,k)
!                   if (dim.eq.3) then
!                      h = h + (stress(i,1,3) + stress(j,1,3))*dwdx(3,k)
!                   endif
                endif            
             elseif (d.eq.2) then
!     y-coordinate of acceleration
               h =    (stress(i,1,2) + stress(j,1,2))*dwdx(1,k)+ (stress(i,2,2) + stress(j,2,2))*dwdx(2,k)
!               if (dim.eq.3) then
!                 h = h + (stress(i,2,3) + stress(j,2,3))*dwdx(3,k)
!               endif    
!             elseif (d.eq.3) then
!!     z-coordinate of acceleration
!               h = (stress(i,1,3) + stress(j,1,3))*dwdx(1,k) + (stress(i,2,3) + stress(j,2,3))*dwdx(2,k)+ (stress(i,3,3) + stress(j,3,3))*dwdx(3,k)
!             endif    
!		          
             h = h*rhoij
             if(isNaN(h)) then
                write(*,'(A)') 'viscous force:h is NaN!!'
                pause
             endif
             dvxdt(d,i) = dvxdt(d,i) + mass(j)*h !dvxdt is beyond limit of number
             if(isNaN(dvxdt(d,i))) then
                write(*,'(A)') 'viscous force:dvxdt(d,i) is NaN!!'
                pause
             endif
             dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
             endif
          enddo     
!----     For SPH algorithm 2
        else if (pa_sph.eq.2) then 
          do d=1,dim                
                if (d.eq.1) then
!     x-coordinate of acceleration
                     h = (stress(i,1,1)/rho(i)**2+ stress(j,1,1)/rho(j)**2)*dwdx(1,k)
                    if(isNaN(h)) then
                          write(*,'(A)') 'viscous force:h is NaN!!'
                     pause
                     endif

                    if (dim.ge.2) then
                           h = h + (stress(i,1,2)/rho(i)**2 + stress(j,1,2)/rho(j)**2)*dwdx(2,k)
!                            h = (stress(i,1,2)/rho(i)**2 + stress(j,1,2)/rho(j)**2)*dwdx(2,k)
                    if(isNaN(h)) then
                           write(*,'(A)') 'viscous force:h is NaN!!'
                    pause
                    endif
!                 if (dim.eq.3) then
!                   h = h + (stress(i,1,3)/rho(i)**2 + stress(j,1,3)/rho(j)**2)*dwdx(3,k)
!                   if(isNaN(h)) then
!                   write(*,'(A)') 'viscous force:h is NaN!!'
!                   pause
!                 endif
!                 endif
                    endif
                   endif
           
           
              if (d.eq.2) then
!     y-coordinate of acceleration
   
               h = (stress(i,2,1)/rho(i)**2 +  stress(j,2,1)/rho(j)**2)*dwdx(1,k) +(stress(i,2,2)/rho(i)**2 + stress(j,2,2)/rho(j)**2)*dwdx(2,k)

                  if(isNaN(h)) then
                       write(*,'(A)') 'viscous force:h is NaN!!'
                       pause
                  endif
!               if (dim.eq.3) then
!                 h = h + (stress(i,2,3)/rho(i)**2+  stress(j,2,3)/rho(j)**2)*dwdx(3,k)
!                 if(isNaN(h)) then
!                  write(*,'(A)') 'viscous force:h is NaN!!'
!                  pause
!               endif   
!               endif           
!             elseif (d.eq.3) then
!!     z-coordinate of acceleration
!               h =     (stress(i,1,3)/rho(i)**2 + stress(j,1,3)/rho(j)**2)*dwdx(1,k)+ (stress(i,2,3)/rho(i)**2 +stress(j,2,3)/rho(j)**2)*dwdx(2,k)+ (stress(i,3,3)/rho(i)**2 +stress(j,3,3)/rho(j)**2)*dwdx(3,k) 
!               if(isNaN(h)) then
!                  write(*,'(A)') 'viscous force:h is NaN!!'
!                  pause
!               endif           
!             endif                       

                endif
!             endif
                  dvxdt(d,i) = dvxdt(d,i) + mass(j)*h
                  dvxdt(d,j) = dvxdt(d,j) - mass(i)*h
          enddo
!
!        endif
        endif
!
      enddo
!
	end