      subroutine kernel(r,dx,hsml,w,dwdx)
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing kernel wij and its 
!   derivatives dwdxij.
!     if skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
!            = 2, Gauss kernel   (Gingold and Monaghan 1981) 
!            = 3, Quintic kernel (Morris 1997)

!     r    : Distance between particles i and j                     [in]
!     dx   : x-, y- and z-distance between i and j                  [in]  
!     hsml : Smoothing length                                       [in]
!     w    : Kernel for all interaction pairs                      [out]
!     dwdx : Derivative of kernel with respect to x, y and z       [out]
      use config_parameter
      implicit none
!      
      real(8) r, dx(dim), hsml, w, dwdx(dim)
      
      real(8) i, j, d
      real(8) q, dw, factor
      real(8) a, k
!
      q = r/hsml 
      w = 0.e0
      do d=1,dim         
        dwdx(d) = 0.e0
      enddo
!----
      if (skf.eq.1) then  !Cubic spline kernel
!      
        if (dim.eq.1) then
          factor = 1.e0/hsml
        elseif (dim.eq.2) then
          factor = 15.e0/(7.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 3.e0/(2.e0*pi*hsml*hsml*hsml)
        else
         print *,' >>> Error <<< : Wrong dimension: Dim =',dim
         pause
        endif     
                                              
        if (q.ge.0.and.q.le.1.e0) then          
          w = factor * (2./3. - q*q + q**3 / 2.)
          do d = 1, dim
            dwdx(d) = factor * (-2.+3./2.*q)/hsml**2 * dx(d)       
          enddo   
        else if (q.gt.1.e0.and.q.le.2) then          
             w = factor * 1.e0/6.e0 * (2.-q)**3 
             do d = 1, dim
            dwdx(d) =-factor * 1.e0/6.e0 * 3.*(2.-q)**2/hsml * (dx(d)/r)
             enddo              
	  else
	    w=0.
          do d= 1, dim
            dwdx(d) = 0.
          enddo             
        endif     
!---                                    
      else if (skf.eq.2) then  !Gauss kernel
!      
        factor = 1.e0 / (hsml**dim * pi**(dim/2.))      
	if(q.ge.0.and.q.le.3) then
	  w = factor * exp(-q*q)
          do d = 1, dim
            dwdx(d) = w * ( -2.* dx(d)/hsml/hsml)
          enddo 
	else
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo 	   
	endif	       
!----	
      else if (skf.eq.3) then  !Quintic kernel
!      
        if (dim.eq.1) then
          factor = 1.e0 / (120.e0*hsml)
        elseif (dim.eq.2) then
          factor = 7.e0 / (478.e0*pi*hsml*hsml)
        elseif (dim.eq.3) then
          factor = 1.e0 / (120.e0*pi*hsml*hsml*hsml)
        else
         write(*,'(A,2x,I2)') ' >>> Error <<< : Wrong dimension: Dim =', dim
         pause
        endif              
	if(q.ge.0.and.q.le.1) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 + 15*(1-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * ( -5*(3-q)**4 + 30*(2-q)**4 - 75*(1-q)**4 ) / hsml * (dx(d)/r) 
            if(r.le.1e-5) then
               write(*,'(A)') 'kernel: two particles are too close!!'
            endif        
         enddo 
	else if(q.gt.1.and.q.le.2) then
          w = factor * ( (3-q)**5 - 6*(2-q)**5 )
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4 + 30*(2-q)**4) / hsml * (dx(d)/r)
            if(r.le.1e-5) then
               write(*,'(A)') 'kernel: two particles are too close!!'
            endif  
          enddo 
        else if(q.gt.2.and.q.le.3) then
          w = factor * (3-q)**5 
          do d= 1, dim
            dwdx(d) = factor * (-5*(3-q)**4) / hsml * (dx(d)/r)
            if(r.le.1e-5) then
               write(*,'(A)') 'kernel: two particles are too close!!'
            endif  
          enddo 
        else   
	  w = 0.
          do d = 1, dim
            dwdx(d) = 0.
          enddo  
        endif
!---
      else if (skf.eq.4) then  !Yang kernel
        a = -10.0/((k+31)*pi*hsml*hsml)
        k = -1.
        if(q<1) then
         w = a*(k*q**3 - 3.0*(k+1)*q*q + 3.0*(k+3)*q - k - 7)
         a = a*(3.0*k*q*q - 6.0*(k+1)*q + 3.0*(k+3))/hsml/r
         dwdx = a*dx
        else
         w = a*(q-2.0)**3
         a = a*3.0*(q-2.0)**2/hsml/r
         dwdx = a*dx
        end if
                           
!                
      endif
      end