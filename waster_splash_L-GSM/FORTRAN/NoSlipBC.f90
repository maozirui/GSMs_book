  subroutine no_slipbc(ntotal,nvirt,x,hsml,vx,p,stress,itype,rho)
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
      real(8) ntotal, itype(maxn), niac,pair_i(max_interaction), pair_j(max_interaction),nvirt,stress(ntotal+nvirt, dim, dim),p(ntotal+nvirt)
      real(8) mass(maxn), x(dim,maxn), hsml(maxn),dvxdt(dim,maxn),vx(dim,maxn),indvxdt(dim,maxn),rho(maxn)
      real(8) i, j, k, d, kn, Pr, E, Fn, dn, mu, Ve, delta1, rr1,delta2,Xe,deltax,deltay
      real(8) dx(dim), rr(ntotal), f, rr0, dd, n1, n2,p1,p2,index,min,beta
!      real(8) rij,near_r(ntotal+nvirt),near_number(ntotal+nvirt)
      integer(4) d1,d2
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!虚粒子速度为反
         beta=1.5
!        

         
         do i=ntotal+1, ntotal+nvirt
              index=0.
               rr0=2.0*hsml(i)
               min=rr0

            do j=1, ntotal
               rr(j)=sqrt((x(1,i)-x(1,j))*(x(1,i)-x(1,j))+(x(2,i)-x(2,j))*(x(2,i)-x(2,j)))

               if(rr(j).le.rr0)then
                if(rr(j).le.min)then
                min=rr(j)
                index=j
                endif
               endif

              enddo
                 

                 
                if (index.ne.0.) then
                
!                if(i.ge.(ntotal+609)) then    !底面边界
                 p(i)=p(index)
   
             do d1=1,dim
                do d2=1,dim

                  stress(i,d1,d2)=stress(index,d1,d2)


                  enddo
                enddo

                  endif
         
     enddo

      
       
     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!速度为零 
      !do i=ntotal+1, ntotal+nvirt
      !    index=0.
      !    rr0=1.1*hsml(i)
      !    min=rr0
      !
      !    do j=1, ntotal
      !        rr(j)=sqrt((x(1,i)-x(1,j))*(x(1,i)-x(1,j))+(x(2,i)-x(2,j))*(x(2,i)-x(2,j)))
      !
      !        if(rr(j).le.rr0)then
      !            if(rr(j).le.min)then
      !                min=rr(j)
      !                index=j
      !            endif
      !        endif
      !
      !    enddo
      !
      !
      !
      !    if (index.ne.0.) then
      !
      !        p(i)=p(index)
      !
      !        if(i.le.(ntotal+609))then
      !            vx(1,index)=0
      !            vx(2,index)=0
      !        endif
      !
      !        do d1=1,dim
      !            do d2=1,dim
      !
      !                stress(i,d1,d2)=stress(index,d1,d2)
      !
      !            enddo
      !        enddo
      !
      !    endif
      !
      !enddo

end           
