      subroutine link_list(itimestep,ntotal,hsml,x,niac,pair_i,pair_j,w,dwdx,countiac,hsmlvector,itype)
     
!----------------------------------------------------------------------
!   Subroutine to calculate the smoothing funciton for each particle and
!   the interaction parameters used by the SPH algorithm. Interaction 
!   pairs are determined by using a sorting grid linked list  

!     itimestep : Current time step                                 [in]
!     ntotal    : Number of particles                               [in]
!     hsml      : Smoothing Length, same for all particles          [in]
!     x         : Coordinates of all particles                      [in]
!     niac      : Number of interaction pairs                      [out]
!     pair_i    : List of first partner of interaction pair        [out]
!     pair_j    : List of second partner of interaction pair       [out]
!     w         : Kernel for all interaction pairs                 [out]
!     dwdx      : Derivative of kernel with respect to x, y and z  [out]
!     countiac  : Number of neighboring particles                  [out]
      use config_parameter
      implicit none

!     Parameter used for sorting grid cells in the link list algorithm
!     maxngx  : Maximum number of sorting grid cells in x-direction
!     maxngy  : Maximum number of sorting grid cells in y-direction
!     maxngz  : Maximum number of sorting grid cells in z-direction
!     Determining maximum number of sorting grid cells:
!     (For an homogeneous particle distribution:)
!     1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
!     2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
!     3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)
      integer maxngx,maxngy,maxngz
      parameter ( maxngx  = 300 , maxngy  = 300   , maxngz  = 1 )      
      real(8) x_t(dim)
      integer i, j, sumiac, maxiac,miniac, maxp,minp
      integer grid(maxngx,maxngy),xgcell(dim,maxn),gcell(dim),xcell,ycell,zcell,celldata(maxn),minxcell(dim),maxxcell(dim),dnxgcell(dim),dpxgcell(dim),ngridx(dim),ghsmlx(dim),n_particle(maxngx,maxngy)
      real(8) hsml2,dr,r,dx(dim),mingridx(dim),maxgridx(dim),tdwdx(dim), dgeomx(dim),x_max(dim),x_min(dim)
      real(8) mhsml,hsmlvector(maxn),x_length,itype(maxn)
      real(8)  ntotal,niac,pair_i(max_interaction),pair_j(max_interaction), countiac(maxn)
      real(8) hsml(maxn), x(dim,maxn),scale_k, w(max_interaction),dwdx(dim,max_interaction),d,noiac
      integer(4) itimestep
     


      if (skf.eq.1) then 
          scale_k = 2 
      else if (skf.eq.2) then 
          scale_k = 3 
      else if (skf.eq.3) then 
          scale_k = 3 
      else if (skf.eq.4) then 
          scale_k = 2     
      endif 
      x_length=x_maxgeom-x_mingeom
!------------------------------------------------------------------------
! calculate the domain size 
          x_max(1)=x_maxgeom
          x_min(1)=x_mingeom
             if (dim.ge.2) then
                x_max(2) = y_maxgeom
                x_min(2) = y_mingeom
                   if (dim.eq.3) then
                       x_max(dim) = z_maxgeom
                       x_min(dim) = z_mingeom
                   endif
            endif
!-----------------------------------------------------------------------------      
      do i=1,ntotal
          countiac(i) = 0
      enddo

!     Initialize grid:  
      
      call init_grid(ntotal,hsml*scale_k,grid,ngridx,ghsmlx,x_max,&
     x_min,maxgridx,mingridx,dgeomx,n_particle)
!    Position particles on grid and create linked list:
      do i=1,ntotal
         if (itype(i).ne.3) then 
            x_t(1)=x(1,i)
            x_t(2)=x(2,i)
!           x_t(3)=x(3,i)
         !   call grid_geom(i,x_t,ngridx,maxgridx,mingridx,dgeomx,gcell)
            do d=1,dim
               xgcell(d,i) = gcell(d)
            enddo
            celldata(i) = grid(gcell(1),gcell(2))
            grid(gcell(1),gcell(2)) = i
          endif 
      enddo
!     Determine interaction parameters:

      niac = 0
      do i=1,ntotal
        if (itype(i).ne.3) then 
!     Determine range of grid to go through:
        do d=1,dim
          minxcell(d) = 1
          maxxcell(d) = 1
        enddo
        do d=1,dim
          dnxgcell(d) = xgcell(d,i) - ghsmlx(d)
          dpxgcell(d) = xgcell(d,i) + ghsmlx(d)
          !minxcell(d) = max(dnxgcell(d),1)
          !maxxcell(d) = min(dpxgcell(d),ngridx(d))
        enddo
        minxcell(1) = dnxgcell(1)
        minxcell(2) = max(dnxgcell(2),1)
        maxxcell(1)= dpxgcell(1)
        maxxcell(2) = min(dpxgcell(2),ngridx(2))
!     Search grid:
!        do zcell=minxcell(3),maxxcell(3)
          do ycell=minxcell(2),maxxcell(2)
           do xcell=minxcell(1),maxxcell(1)
               
            if (xcell.le.0) then
               j = grid(xcell+ngridx(1),ycell) 
            elseif (xcell.gt.ngridx(1)) then
               j = grid(xcell-ngridx(1),ycell) 
            else    
               j = grid(xcell,ycell)
            endif   

1           if (j.ne.0) then
             if (j.gt.i) then
                 
                  dx(1) = x(1,i) - x(1,j)
                  if (dx(1).ge.x_length/2) then
                       dx(1)=x(1,i) - x(1,j) - x_length
                  elseif (dx(1).le.-x_length/2) then
                       dx(1) = x(1,i) - x(1,j) + x_length
                  endif

                  dr    = dx(1)*dx(1)
                  do d=2,dim
                      dx(d) = x(d,i) - x(d,j)
                      dr    = dr + dx(d)*dx(d)
                  enddo
                  mhsml = (hsmlvector(i)+hsmlvector(j))/2.
                 if (sqrt(dr).lt.scale_k*mhsml) then
                  if (niac.lt.max_interaction) then
!     Neighboring pair list, and totalinteraction number and
!     the interaction number for each particle 
                    niac = niac + 1
                    pair_i(niac) = i
                    pair_j(niac) = j
                    r = sqrt(dr)
                    countiac(i) = countiac(i) + 1
                    countiac(j) = countiac(j) + 1
!--- Kernel and derivations of kernel
                    call kernel(r,dx,mhsml,w(niac),tdwdx)
	                do d = 1, dim
	                    dwdx(d,niac)=tdwdx(d)
                    enddo                  
                  else
                    print *,  ' >>> Error <<< : too many interactions'
                    stop
                  endif
                endif
                endif
                j = celldata(j)
                goto 1
              endif
            enddo
          enddo
          
         endif
        enddo
!      enddo


!     Statistics for the interaction
      sumiac = 0
      maxiac = 0
      miniac = 1000
      noiac  = 0
      do i=1,ntotal
        sumiac = sumiac + countiac(i)
        if (countiac(i).gt.maxiac) then
	  maxiac = countiac(i)
	  maxp = i
	endif
	if (countiac(i).lt.miniac) then 
	  miniac = countiac(i)
          minp = i
	endif
        if (countiac(i).eq.0)      noiac  = noiac + 1
      enddo
 
      if (mod(itimestep,print_step).eq.0) then
        if (int_stat) then
 !         print *,' >> Statistics: interactions per particle:'
 !         print *,'**** Particle:',maxp, ' maximal interactions:',maxiac
 !         print *,'**** Particle:',minp, ' minimal interactions:',miniac
          print *,'**** Average :',real(sumiac)/real(ntotal)
          print *,'**** Total pairs : ',niac
!          print *,'**** Particles with no interactions:',noiac
        endif     
      endif
!      write(*,*) '**** Information for particle ****',dgeomx(1)/ngridx(1)/0.04
!      write(*,*) '**** Information for particle ****',dgeomx(2)/ngridx(2)/0.04
!      write(*,*) '**** Information for particle ****',dgeomx(3)/ngridx(3)/0.04
      end
