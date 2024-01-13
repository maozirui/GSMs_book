      subroutine viscosity(ntotal, itype, rho, eta)
!----------------------------------------------------------------------
!   Subroutine to define the FLUID particle viscosity
 
!     ntotal  : Number of particles                                 [in]
!     itype    : Type of particle                                   [in]
!     x       : Coordinates of all particles                        [in]
!     rho     : Density                                             [in]
!     eta     : Dynamic viscosity                                  [out]
!     Note: it is likely,
!     kinetic viscosity = dynamic viscosity / density
      use config_parameter
      implicit none
!      
      real(8) ntotal, itype(maxn)
      real(8) rho(maxn), eta(maxn)
	real(8) i
!
      do i=1,ntotal
        if (abs(itype(i)).eq.1) then
          eta(i)=0.
        else if (abs(itype(i)).eq.2) then
          eta(i)=1.0e-3 !For water in general situation
	  else if (abs(itype(i)).eq.100) then
          eta(i)=0.
          else if (abs(itype(i)).eq.200) then
          eta(i)=0.5

        else if (abs(itype(i)).eq.7) then
          eta(i)=0. !For explosive gas
        endif
      enddo  
! 
      end