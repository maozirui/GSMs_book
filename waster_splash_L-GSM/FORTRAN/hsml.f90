      subroutine h_upgrade(dt,ntotal,countiac, mass, vx, rho, niac, pair_i, pair_j, dwdx, hsml)
!-----------------------------------------------------------------------
!     Subroutine to evolve smoothing length

!     dt     : time step                                            [in]
!     ntotal : Number of particles                                  [in]
!     mass   : Particle masses                                      [in]
!     vx     : Velocities of all particles                          [in]
!     rho    : Density                                              [in]
!     niac   : Number of interaction pairs                          [in]
!     pair_i : List of first partner of interaction pair            [in]
!     pair_j : List of second partner of interaction pair           [in]
!     dwdx   : Derivative of kernel with respect to x, y and z      [in]
!     hsml   : Smoothing Length                                 [in/out]
      use config_parameter
      implicit none
!
      real(8) ntotal, niac, pair_i(max_interaction),pair_j(max_interaction)
      real(8) mass(maxn), vx(dim, maxn), rho(maxn),dwdx(dim, max_interaction), hsml(maxn)
      real(8) i,j,k,d

      real(8) dt, fac, dvx(dim), hvcc, vcc(maxn), dhsml(maxn)
      real(8) countiac(maxn)
!
      if (sle.eq.0 ) then
!---  Keep smoothing length unchanged.
      return
!
      else if (sle.eq.2) then
!---  dh/dt = (-1/dim)*(h/rho)*(drho/dt). (Benz 1989)
        do i=1,ntotal
          vcc(i) = 0.e0
        enddo
!
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
          vcc(i) = vcc(i) + mass(j)*hvcc/rho(j)
          vcc(j) = vcc(j) + mass(i)*hvcc/rho(i)
        enddo
!
        do i = 1, ntotal
          dhsml(i) = (hsml(i)/dim)*vcc(i)
          hsml(i) = hsml(i) + dt*dhsml(i)
          if (hsml(i).le.0) hsml(i) = hsml(i) - dt*dhsml(i)

        enddo
!
      else if(sle.eq.1) then
!---  h = fac * (m/rho)^(1/dim). (Monaghan 1985)
        fac = 1.0
        do i = 1, ntotal
          hsml(i) = fac * (mass(i)/rho(i))**(1./dim)

        enddo
        

!
      endif
!
      end