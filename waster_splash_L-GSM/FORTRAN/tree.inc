c--------------------------------------------------------------------------------

c     Parameters for the TREE algorithm (S. Pfalzner, see tree.f) and tree search:
c     npartm  : Maximum number of particles
c     ncelm   : Maximum number of cells
c     ntwigm  : Maximum number of twigs
c     nsnapm  :
c     nlevmax : Maximum number of levels
c     plottree: Print the tree structure (output: fort.15)
      integer npartm,ncelm,ntwigm,nsnapm,nlevmax
      logical plottree
        
      parameter ( npartm  = maxn     )
      parameter ( ncelm   = 2*npartm )
      parameter ( ntwigm  = 2*npartm )
      parameter ( nsnapm  = 100      )
      parameter ( nlevmax = 50       )
      parameter ( plottree = .false. )