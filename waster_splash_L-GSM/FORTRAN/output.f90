    subroutine output(x, vx, mass, rho, p, c,itype, hsml, ntotal,itimestep, dt)
    !----------------------------------------------------------------------
    !     Subroutine for saving particle information to external disk file

    !     x-- coordinates of particles                                  [in]
    !     vx-- veloc          ities of particles                        [in]
    !     mass-- mass of particles                                      [in]
    !     rho-- dnesities of particles                                  [in]
    !     p-- pressure  of particles                                    [in]
    !     u-- internal energy of particles                              [in]
    !     c-- sound velocity of particles                               [in]
    !     itype-- types of particles                                    [in]
    !     hsml-- smoothing lengths of particles                         [in]
    !     ntotal-- total particle number                                [in]
    use config_parameter
    implicit none
    !
    real(8) itype(maxn), ntotal,temptype
    real(8) x(dim, maxn), vx(dim, maxn), mass(maxn),rho(maxn),p(maxn), u(maxn), c(maxn), hsml(maxn),dt,f_R(maxn),ahdudt(maxn),stress(maxn, dim, dim),accumulated_strain(maxn)
    real(8) i, d, explosive, shell, obj, gg(maxn),HT,SR,V0,nt
    character (40) :: xname, sname, oname, dname
    integer(4) itimestep



    !      npart1 = 1892
    !      npart2 = 2783
    !      npart3 = 1000
    !      npart4 = 1892
    !      npart5 = 2783
    !      npart6 = 1000

    !      output in paraview data format
    !      if(itimestep<10) then
    !         write(xname, '(A, I1, A)') 'data1/t=00000', itimestep,'step_xv.vtk'
    !      else if(itimestep>9.and.itimestep<100) then
    !         write(xname, '(A, I2, A)') 'data1/t=0000', itimestep,'step_xv.vtk'
    !      else if(itimestep>99.and.itimestep<1000) then
    !         write(xname, '(A, I3, A)') 'data1/t=000', itimestep,'step_xv.vtk'
    !      else if(itimestep>999.and.itimestep<10000) then
    !         write(xname, '(A, I4, A)') 'data1/t=00', itimestep,'step_xv.vtk'
    !      else if(itimestep>9999.and.itimestep<100000) then
    !         write(xname, '(A, I5, A)') 'data1/t=0', itimestep,'step_xv.vtk'
    !      endif

    !      open (1, file = xname)

    !      write (1, '(A/A/A//A)') '# vtk DataFile Version 2.0', 'step_xv', 'ASCII','DATASET UNSTRUCTURED_GRID'

    !      if (ntotal.ne.0) then
    !       write (1,'(A,2X,I6,2X,A)') 'POINTS',ntotal, 'float'
    !        do i=1,ntotal
    !            write(1,'(F16.8,2X,F16.8,2X,F16.8)') (x(d, i), d=1,dim)
    !        enddo

    !        write (1,'(/A,2X,I6,2X,I8)') 'CELLS',ntotal,2*ntotal
    !        do i=1,ntotal
    !            write(1,'(I4,2X,I6)')1, i
    !        enddo

    !        write (1,'(/A,2X,I6)') 'CELL_TYPES',ntotal
    !        do i=1,ntotal
    !            write(1,'(I2)') 2
    !        enddo

    !        write (1,'(/A,2X,I6/,A,2X,A,2X,A)') 'POINT_DATA',ntotal,'VECTORS','velocity','float'
    !        do i=1,ntotal
    !           write(1,'(F16.8,2X,F16.8,2X,F16.8)') (vx(d, i), d = 1, dim)
    !        enddo

    !      write (1,'(/A,2X,A,2X,A,2X,I2/,A,2X,A)') 'SCALARS','rho','float',1, 'LOOKUP_TABLE','default'
    !        do i=1,ntotal
    !            write(1,'(F16.8)') rho(i)
    !        enddo

    !      write (1,'(/A,2X,A,2X,A,2X,I2/,A,2X,A)') 'SCALARS','mass','float',1,'LOOKUP_TABLE','default'
    !        do i=1,ntotal
    !            write(1,'(F16.8)') mass(i)
    !        enddo

    !        write (1,'(/A,2X,A,2X,A,2X,I2/,A,2X,A)') 'SCALARS','p','float',1,'LOOKUP_TABLE','default'
    !        do i=1,ntotal
    !            write(1,'(E16.8)') p(i)
    !        enddo

    !      write (1,'(/A,2X,A,2X,A,2X,I2/,A,2X,A)') 'SCALARS','itype','float',1, 'LOOKUP_TABLE','default'
    !        do i=1,ntotal
    !            write(1,'(I3)') itype(i)
    !        enddo

    !      write (1,'(/A,2X,A,2X,A,2X,I2/,A,2X,A)') 'SCALARS','hsml','float',1,'LOOKUP_TABLE','default'
    !        do i=1,ntotal
    !            write(1,'(F16.8)') hsml(i)
    !        enddo
    !      end if

    !------------------------------------------------------------------------
    !     output result in tecplot format
    !if(itimestep>0.and.itimestep<10) then
    !    write(xname, '(A, I1, A)') 'data1/t=000000', itimestep,'step_xv.dat'
    !else if(itimestep>9.and.itimestep<100) then
    !    write(xname, '(A, I2, A)') 'data1/t=00000', itimestep,'step_xv.dat'
    !else if(itimestep>99.and.itimestep<1000) then
    !    write(xname, '(A, I3, A)') 'data1/t=0000', itimestep,'step_xv.dat'
    !else if(itimestep>999.and.itimestep<10000) then
    !    write(xname, '(A, I4, A)') 'data1/t=000', itimestep,'step_xv.dat'
    !else if(itimestep>9999.and.itimestep<100000) then
    !    write(xname, '(A, I5, A)') 'data1/t=00', itimestep,'step_xv.dat'
    !else if(itimestep>99999.and.itimestep<1000000) then
    !    write(xname, '(A, I6, A)') 'data1/t=0', itimestep,'step_xv.dat'
    !else if(itimestep>999999.and.itimestep<10000000) then
    !    write(xname, '(A, I7, A)') 'data1/t=', itimestep,'step_xv.dat'
    !endif
    
    write(xname, '(A, F10.7, A)') 'data1/t=', dt,'sec_xv.dat'
    
    open (1, file = xname)

    write (1, '(A, I7, A)') ' title="n=', itimestep, ' s"'

    ! for 3D output
    !      write(1,*) 'variables = x, y, z, vx, vy, vz, rho, mass, itype, p, u, hsml, f_R'
    ! for 2D output
    write(1,*) 'variables = x[m], y[m], vx[m/s], vy[m/x], rho, p,itype'
    ! for 1D output
    !      write(1,*) 'variables = x, vx, rho, mass, itype, p, u, hsml, f_R'

    !      do i=1,ntotal
    !         if(x(1,i).ge.0.01) then
    !            explosive=explosive+1
    !         elseif (itype(i).eq.102) then
    !            shell=shell+1
    !         elseif (itype(i).eq.104) then
    !            obj=obj+1
    !         endif
    !      enddo
    nt=ntotal
    do i=1,ntotal
        if (itype(i)==50) nt=nt-1
    enddo
    
    if (ntotal.ne.0) then
        write (1,'(/A,(F9.0),A10)') 'ZONE i =',nt, 'f = point'
        do i=1,ntotal
            !           if(itype(i)==7) then
            !              temptype=7
            !           else if(itype(i)==102) then
            !              temptype=
            
            !           else
            !              temptype=104
            !           end if

            if (itype(i).ne.50) write(1,1001) (x(d, i), d=1,dim), (vx(d, i), d = 1, dim), rho(i), p(i),itype(i)
            !            write(2,1002) mass(i), rho(i), p(i), u(i)
            !            write(3,1003) itype(i), hsml(i)
        enddo
    end if
        write(xname, '(A, F10.7, A)') 'data2/t=', dt,'sec_xv.dat'
    
    open (1, file = xname)

    write (1, '(A, I7, A)') ' title="n=', itimestep, ' s"'

    write(1,*) 'variables = x[m], y[m], vx[m/s], vy[m/x], rho, p,itype'
      d=0
    do i=1,ntotal
        if (mass(i).eq.200) d=d+1
    enddo
    
    if (ntotal.ne.0) then
        write (1,'(/A,(F9.0),A10)') 'ZONE i =',d, 'f = point'
        do i=1,ntotal
            if (mass(i).eq.200) write(1,1001) (x(d, i), d=1,dim), (vx(d, i), d = 1, dim), rho(i), p(i),itype(i)
        enddo
    end if
    
    !----------------------------------------------------------------------------------
    !        write (1,'(/A,I6,A10)') 'ZONE i =',shell, 'f = point'
    !        do i=1,ntotal
    !            if(itype(i).eq.102) then
    !            write(1,1001) (x(d, i), d=1,dim), (vx(d, i), d = 1, dim),p(i)
    !            write(2,1002) mass(i), rho(i), p(i), u(i)
    !            write(3,1003) itype(i), hsml(i)
    !            endif
    !        enddo

    !        write (1,'(/A,I6,A10)') 'ZONE i =',obj, 'f = point'
    !        do i=1,ntotal
    !            if(itype(i).eq.104) then
    !            write(1,1001) (x(d, i), d=1,dim), (vx(d, i), d = 1, dim), p(i)
    !            write(2,1002) mass(i), rho(i), p(i), u(i)
    !            write(3,1003) itype(i), hsml(i)
    !            endif
    !        enddo
    !      end if


    !1001  format(8(e22.10,2x),I5,4(2x,e22.10))
1001 format(7(e22.10,2x))
    !1001  format(13(e22.10,2x))
    !1002  format(10(2x, e14.8))
    !1003  format(10(2x, e14.8))

    close(1)

    ! output the SR and HT in the selected step
    HT=0
    SR=0
    V0=1000
    do i=1,ntotal
        if (itype(i)>50.and.x(2,i)>HT) then
            HT=x(2,i)
        endif

        if (itype(i)>50.and.x(1,i)>SR) then
            SR=x(1,i)
        endif
        
        if (flowtype==3) then
            if (itype(i)>50.and.abs(x(1,i)-12.75)<abs(x(1,V0)-12.75)) V0=i
        endif 
    enddo
    V0=vx(1,V0)
    
    open(unit=10, file = "0HS.txt")
    if (flowtype<3) then
        write(10,1001) dt, HT, SR
    elseif (flowtype==3) then
        write(10,1001) dt, HT, V0
    endif
1002 format(1x,3(2x,e15.8))

    !      open (1, file='data1/shelldev.txt',status='old',position='append')
    !      write(1,'(e22.10,2x,e22.10)') dt*itimestep, maxvy
    !      close(1)

    end
