      subroutine art_heat(ntotal,hsml,mass,x,vx,niac,rho,u,c,pair_i,pair_j,w,dwdx,dedt)

!c----------------------------------------------------------------------
!c     Subroutine to calculate the artificial heat(Fulk, 1994, p, a-17) 
!
!c     ntotal : Number of particles                                  [in]
!c     hsml   : Smoothing Length                                     [in]
!c     mass   : Particle masses                                      [in]
!c     x      : Coordinates of all particles                         [in]
!c     vx     : Velocities of all particles                          [in]
!c     rho    : Density                                              [in]
!c     u      : specific internal energy                             [in]
!c     c      : Sound veolcity                                       [in]
!c     niac   : Number of interaction pairs                          [in]
!c     pair_i : List of first partner of interaction pair            [in]
!c     pair_j : List of second partner of interaction pair           [in]
!c     w      : Kernel for all interaction pairs                     [in]
!c     dwdx   : Derivative of kernel with respect to x, y and z      [in]
!c     dedt   : produced artificial heat, adding to energy Eq.      [out]
      use config_parameter
      implicit none

   
      real(8) ntotal,niac,pair_i(max_interaction),pair_j(max_interaction)
      real(8) hsml(maxn), mass(maxn), x(dim,maxn),vx(dim,maxn),rho(maxn), u(maxn), c(maxn),w(max_interaction),dwdx(dim,max_interaction), dedt(maxn)
      real(8) i,j,k,d
      real(8) dx,dvx(dim),vr,rr,h,mc,mrho,mhsml,vcc(maxn), hvcc, mui,muj,muij,rdwdx,g1,g2,energy(maxn),lmass(maxn)
      parameter( g1 = 0.1)
!     Shear viscosity
      parameter( g2 = 0.  )
                     
!---  Parameter for the artificial heat conduction:
     
     mc=600
      do i=1,ntotal
        vcc(i) = 0.e0         
        dedt(i) = 0.e0
      enddo
    
      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        do d=1,dim
          dvx(d) = vx(d,j) - vx(d,i) 
        enddo        
        hvcc = dvx(1)*dwdx(1,k)
        do d=2,dim
          hvcc = hvcc + dvx(d)*dwdx(d,k)
        enddo    
         vcc(i) = vcc(i) + mass(j)*hvcc/rho(i)
         vcc(j) = vcc(j) + mass(i)*hvcc/rho(j)
      enddo  
   

      do k=1,niac
        i = pair_i(k)
        j = pair_j(k)
        mhsml= (hsml(i)+hsml(j))/2.
        mrho = 0.5e0*(rho(i) + rho(j))                          
        rr = 0.e0
        rdwdx = 0.e0
        do d=1,dim
          dx = x(d,i) -  x(d,j)
          rr = rr + dx*dx
          rdwdx  = rdwdx + dx*dwdx(d,k)            
        enddo             
        mui=g1*hsml(i)*mc*abs(vcc(i)) + g2*hsml(i)**2*(vcc(i))**2
        muj=g1*hsml(j)*mc*abs(vcc(j)) + g2*hsml(j)**2*(vcc(j))**2
        muij= 0.5*(mui+muj)     
        h = muij/(mrho*(rr+0.01*mhsml**2))*rdwdx

        dedt(i) = dedt(i) + mass(j)*h*(u(i)-u(j))
        dedt(j) = dedt(j) + mass(i)*h*(u(j)-u(i))
      enddo
  
      do i=1,ntotal
        dedt(i) = 2.0e0*dedt(i)          
      enddo

      end