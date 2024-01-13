!==========================================================================
!                  ºËÌÝ¶ÈÐÞÕý(KGC)_correct
!--------------------------------------------------------------------------
      subroutine kernelgrad_correct(ntotal,niac,x,hsml,pair_i,pair_j,itype,rho,mass,dwdx)
     
!     ntotal    : Number of particles                               [in]
!     niac      : Number of interaction pairs                       [in]
!     x         : Coordinates of all particles                      [in]
!     hsml      : Smoothing Length                                  [in]
!     pair_i : List of first partner of interaction pair            [in]
!     pair_j : List of second partner of interaction pair           [in]
!     mass      : Particle masses                                   [in]
!     rho       : Density                                           [in]
!     x       : Coordinates of all particles                        [in]
!     itype     : type of particles                                 [in]
!     dwdx      : Derivative of kernel with respect to x, y and z  [out]
      use config_parameter
      implicit none

      real(8) i, j, k, d, ntotal,niac,maxp
      real(8) itype(maxn)
      real(8) rr2, hsml(maxn),mhsml, vol_i,x(dim,maxn),vol_j,detM,rho(maxn),mass(maxn)
      real(8) ndwdx(dim),drxx(dim), dwdx_n(dim),one_over_detM 
      real(8) ww 
      real(8) aM_a11(maxn),aM_a12(maxn),aM_a13(maxn)
      real(8) aM_a21(maxn),aM_a22(maxn),aM_a23(maxn)
      real(8) aM_a31(maxn),aM_a32(maxn),aM_a33(maxn)
      real(8) aL_a11(maxn),aL_a12(maxn),aL_a13(maxn)
      real(8) aL_a21(maxn),aL_a22(maxn),aL_a23(maxn)
      real(8) aL_a31(maxn),aL_a32(maxn),aL_a33(maxn)
      real(8) dwdxi(dim,max_interaction),dwdxj(dim,max_interaction)
      real(8)  dwdx(dim,max_interaction)
      real(8) pair_i(max_interaction),pair_j(max_interaction)

       !dwdx(1:dim,1:max_interaction) = 0.0
      do i=1,ntotal
!       Kernel ReNormalisation ARRAYS --
      aM_a11(i)=0.    !Mi(1,1)
      aM_a12(i)=0.    !Mi(1,2)
      aM_a13(i)=0.    !Mi(1,3)
      aM_a21(i)=0.    !Mi(2,1)
      aM_a22(i)=0.    !Mi(2,2)
      aM_a23(i)=0.    !Mi(2,3)
      aM_a31(i)=0.    !Mi(3,1)
      aM_a32(i)=0.    !Mi(3,2)
      aM_a33(i)=0.    !Mi(3,3)
      aL_a11(i)=0.    !Li(1,1)
      aL_a12(i)=0.    !Li(1,2)
      aL_a13(i)=0.    !Li(1,3)
      aL_a21(i)=0.    !Li(2,1)
      aL_a22(i)=0.    !Li(2,2)
      aL_a23(i)=0.    !Li(2,3)
      aL_a31(i)=0.    !Li(3,1)
      aL_a32(i)=0.    !Li(3,2)
      aL_a33(i)=0.    !Li(3,3)
      end do

      do k=1,niac
      i = pair_i(k)
      j = pair_j(k)

      drxx(1) = x(1,i) - x(1,j)
      drxx(2) = x(2,i) - x(2,j)
!      drxx(3) = x(3,i) - x(3,j)
      rr2 =sqrt(drxx(1)*drxx(1) + drxx(2)*drxx(2))
      mhsml = (hsml(i)+hsml(j))/2
      call kernel(rr2, drxx, mhsml,ww ,ndwdx)

      vol_j = mass(j)/rho(j) 


      do d=1,dim				                         
      dwdx_n(d) = vol_j*ndwdx(d)
      end do

      aM_a11(i) = aM_a11(i) - dwdx_n(1)*drxx(1)
      aM_a12(i) = aM_a12(i) - dwdx_n(1)*drxx(2)
!      aM_a13(i) = aM_a13(i) - dwdx_n(1)*drxx(3)
      aM_a21(i) = aM_a21(i) - dwdx_n(2)*drxx(1)
      aM_a22(i) = aM_a22(i) - dwdx_n(2)*drxx(2)
!      aM_a23(i) = aM_a23(i) - dwdx_n(2)*drxx(3)
!      aM_a31(i) = aM_a31(i) - dwdx_n(3)*drxx(1)
!      aM_a32(i) = aM_a32(i) - dwdx_n(3)*drxx(2)
!      aM_a33(i) = aM_a33(i) - dwdx_n(3)*drxx(3)

      vol_i = mass(i)/rho(i) 

      do d=1,dim				                         
      dwdx_n(d) = vol_i*ndwdx(d)
      end do

      aM_a11(j) = aM_a11(j) - dwdx_n(1)*drxx(1)
      aM_a12(j) = aM_a12(j) - dwdx_n(1)*drxx(2)
!      aM_a13(j) = aM_a13(j) - dwdx_n(1)*drxx(3)
      aM_a21(j) = aM_a21(j) - dwdx_n(2)*drxx(1)
      aM_a22(j) = aM_a22(j) - dwdx_n(2)*drxx(2)
!      aM_a23(j) = aM_a23(j) - dwdx_n(2)*drxx(3)
!      aM_a31(j) = aM_a31(j) - dwdx_n(3)*drxx(1)
!      aM_a32(j) = aM_a32(j) - dwdx_n(3)*drxx(2)
!      aM_a33(j) = aM_a33(j) - dwdx_n(3)*drxx(3)

      end do

      do i=1,ntotal
      if(i.ge.1.and.i.le.ntotal)then
!     M should be symmetric
	  aM_a12(i)=0.5*(aM_a12(i)+aM_a21(i))
      aM_a13(i)=0.5*(aM_a13(i)+aM_a31(i))
      aM_a23(i)=0.5*(aM_a23(i)+aM_a32(i))

	  aM_a21(i)=aM_a12(i)
      aM_a31(i)=aM_a13(i)
      aM_a32(i)=aM_a23(i)
!     Determinant of matrix M det(M)
      detM=aM_a11(i)*aM_a22(i)*aM_a33(i)+aM_a12(i)*aM_a23(i)*aM_a31(i)+aM_a21(i)*aM_a32(i)*aM_a13(i)-aM_a13(i)*aM_a22(i)*aM_a31(i)-aM_a12(i)*aM_a21(i)*aM_a33(i)-aM_a23(i)*aM_a32(i)*aM_a11(i)

      if(abs(detM).gt.0.01.and.abs(aM_a11(i)).gt.0.25.and.abs(aM_a22(i)).gt.0.25.and.abs(aM_a33(i)).gt.0.25) then
	
	  one_over_detM = 1.0/detM
!     Get the inversion of aM
      aL_a11(i)=(aM_a22(i)*aM_a33(i)-aM_a23(i)*aM_a32(i))*one_over_detM
      aL_a22(i)=(aM_a11(i)*aM_a33(i)-aM_a13(i)*aM_a31(i))*one_over_detM
      aL_a33(i)=(aM_a11(i)*aM_a22(i)-aM_a12(i)*aM_a21(i))*one_over_detM
      aL_a12(i)=-(aM_a21(i)*aM_a33(i)-aM_a23(i)*aM_a31(i))*one_over_detM
      aL_a13(i)=(aM_a21(i)*aM_a32(i)-aM_a22(i)*aM_a31(i))*one_over_detM
      aL_a21(i)=-(aM_a12(i)*aM_a33(i)-aM_a13(i)*aM_a32(i))*one_over_detM
      aL_a23(i)=-(aM_a11(i)*aM_a32(i)-aM_a12(i)*aM_a31(i))*one_over_detM
      aL_a31(i)=(aM_a12(i)*aM_a23(i)-aM_a13(i)*aM_a22(i))*one_over_detM
      aL_a32(i)=-(aM_a11(i)*aM_a23(i)-aM_a13(i)*aM_a21(i))*one_over_detM

      else
!     No correction 
      aL_a11(i) = 1.
      aL_a12(i) = 0.
      aL_a13(i) = 0.
      aL_a21(i) = 0.
	  aL_a22(i) = 1.
      aL_a23(i) = 0.
      aL_a31(i) = 0.
	  aL_a32(i) = 0.
      aL_a33(i) = 1.
      endif	  

      else
!	 No correction 
      aL_a11(i) = 1.
      aL_a12(i) = 0.
      aL_a13(i) = 0.
      aL_a21(i) = 0.
	  aL_a22(i) = 1.
      aL_a23(i) = 0.
      aL_a31(i) = 0.
	  aL_a32(i) = 0.
      aL_a33(i) = 1.
      endif
      enddo

      do k=1,niac
      i = pair_i(k)
      j = pair_j(k)

      if(i.le.ntotal.and.j.le.ntotal)then
!     dwdx(1,k)=aL_a11(i)*dwdx(1,k)+aL_a12(i)*dwdx(2,k)
!     dwdx(1,k)=aL_a11(i)*dwdx(1,k)+aL_a12(i)*dwdx(2,k)

         dwdxi(1,k)=aL_a11(i)*dwdx(1,k)+aL_a12(i)*dwdx(2,k)
         dwdxi(2,k)=aL_a21(i)*dwdx(1,k)+aL_a22(i)*dwdx(2,k)
!         dwdxi(3,k)=aL_a31(i)*dwdx(1,k)+aL_a32(i)*dwdx(2,k)+aL_a33(i)*dwdx(3,k)

	  dwdxj(1,k)=aL_a11(j)*dwdx(1,k)+aL_a12(j)*dwdx(2,k)
         dwdxj(2,k)=aL_a21(j)*dwdx(1,k)+aL_a22(j)*dwdx(2,k)
!         dwdxj(3,k)=aL_a31(j)*dwdx(1,k)+aL_a32(j)*dwdx(2,k)+aL_a33(j)*dwdx(3,k)

	  dwdx(1,k)=(dwdxi(1,k)+dwdxj(1,k))/2
	  dwdx(2,k)=(dwdxi(2,k)+dwdxj(2,k))/2
!         dwdx(3,k)=(dwdxi(3,k)+dwdxj(3,k))/2
      end if

      end do
      end subroutine

