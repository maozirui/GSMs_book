      subroutine init_grid(ntotal,hsml,grid,ngridx,ghsmlx,x_max,x_min,maxgridx,mingridx,dgeomx,n_particle)

!c----------------------------------------------------------------------      
!c   Subroutine to established a pair linked list by sorting grid cell.
!c   It is suitable for a homogeneous particle distribution with the 
!c   same smoothing length in an instant. A fixed number of particles
!c   lie in each cell. 
!
!c     ntotal   : Number of particles                                [in]
!c     hsml     : Smoothing Length                                   [in]
!c     grid     : array of grid cells                               [out]
!c     ngridx   : Number of sorting grid cells in x, y, z-direction [out]
!c     ghsmlx   : Smoothing length measured in cells of the grid    [out]
!c     maxgridx : Maximum x-, y- and z-coordinate of grid range     [out]
!c     mingridx : Minimum x-, y- and z-coordinate of grid range     [out]
!c     dgeomx   : x-, y- and z-expansion of grid range              [out]


!c     Parameter used for sorting grid cells in the link list algorithm
!c     maxngx  : Maximum number of sorting grid cells in x-direction
!c     maxngy  : Maximum number of sorting grid cells in y-direction
!c     maxngz  : Maximum number of sorting grid cells in z-direction
!c     Determining maximum number of sorting grid cells:
!c     (For an homogeneous particle distribution:)
!c     1-dim. problem: maxngx = maxn ,  maxngy = maxngz = 1
!c     2-dim. problem: maxngx = maxngy ~ sqrt(maxn) ,  maxngz = 1
!c     3-dim. problem: maxngx = maxngy = maxngz ~ maxn^(1/3)

      use config_parameter
      implicit none
      integer maxngx,maxngy,maxngz
      parameter ( maxngx  = 300   ,maxngy  = 300  ,maxngz  = 1   )
      integer grid(maxngx,maxngy), ngridx(dim),ghsmlx(dim),n_particle(maxngx,maxngy)
      real(8) maxgridx(dim), mingridx(dim), dgeomx(dim)
      integer i, j, k, maxng(dim), ngrid(dim),mp,sort_free
      real(8) nppg,x1,dx_int,x_max(dim),x_min(dim)
      real(8) ntotal,hsml,d

!     Averaged number of particles per grid cell

!      parameter( nppg = 2.01e0 )
!      
!       dx_int = 0.02/1000*1.2

!     Initialize parameters: Maximum number of grid cells
      maxng(1) = maxngx
      if (dim.ge.2) then
        maxng(2) = maxngy
        if (dim.eq.3) then
           maxng(dim) = maxngz
        endif
      endif
      
      do d=1,dim
          ngrid(d) = 1
      enddo
      
!     Range of sorting grid

      maxgridx(1) = x_max(1)
      mingridx(1) = x_min(1)
      if (dim.ge.2) then
        maxgridx(2) = x_max(2)
        mingridx(2) = x_min(2)
        if (dim.eq.3) then
          maxgridx(dim) = x_max(dim)
          mingridx(dim) = x_min(dim)
        endif
      endif
      
      do d=1,dim
         dgeomx(d) = maxgridx(d) - mingridx(d)
      enddo
      
!     Number of grid cells in x-, y- and z-direction:
      if (dim.eq.1) then
        ngridx(1) = min(int(ntotal/nppg) + 1,maxng(1))
      else if (dim.eq.2) then
        ngridx(1) = min(int(dgeomx(1)/hsml),maxng(1))
        ngridx(2) = min(int(dgeomx(2)/hsml),maxng(2))
!        ngridx(3)=1
      endif
      
      !do d=1,dim
      !   dgeomx(d) = ngridx(d)*hsml
      !enddo
      
      maxgridx(1) = mingridx(1)+dgeomx(1)
      if (dim.ge.2) then
        maxgridx(2) = mingridx(2)+dgeomx(2)
        if (dim.eq.3) then
          maxgridx(dim) = mingridx(dim)+dgeomx(dim)
        endif
      endif
      
      
!    Smoothing Length measured in grid cells:

      do d=1,dim
         ghsmlx(d) = int(real(ngridx(d))*hsml/dgeomx(d)) + 1
      enddo

      do d=1,dim
        ngrid(d) = ngridx(d)
      enddo

!     Initialize grid
          
      do i=1,ngrid(1)
        do j=1,ngrid(2)
!          do k=1,ngrid(3)
              grid(i,j) = 0
              n_particle(i,j)=0
!          enddo
        enddo
      enddo
     
!      write(*,*) '**** Information for particle ****',dgeomx(1)/ngridx(1)/dx_int
      end
