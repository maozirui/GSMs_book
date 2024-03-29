    subroutine triangulation_order3(x, itimestep,node_num, itype, hsml, countiac, numb )

    !*****************************************************************************80
    !
    !! MAIN is the main program for TABLE_DELAUNAY.
    !
    !  Discussion:
    !
    !    TABLE_DELAUNAY computes the Delaunay triangulation of a TABLE dataset.
    !
    !    The dataset is simply a set of points in the plane.
    !
    !    Thus, given a set of points V1, V2, ..., VN, we apply a standard
    !    Delaunay triangulation.  The Delaunay triangulation is an organization
    !    of the data into triples, forming a triangulation of the data, with
    !    the property that the circumcircle of each triangle never contains
    !    another data point.
    !
    !  Usage:
    !
    !    table_delaunay prefix
    !
    !    where:
    !
    !    'prefix' is the common prefix for the node and triangle files:
    !
    !    * prefix_nodes.txt,     the node coordinates (input).
    !    * prefix_elements.txt,  the nodes that make up each triangle (output).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 October 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Usage:
    !
    !    table_delaunay input_filename
    !
    use config_parameter
    implicit none
    character (40) :: xname
    character ( len = 80 ) :: file_name = 'triangulation_order3_plot.eps'
    integer ( kind = 4 ) arg_num
    integer ( kind = 4 ) iarg,i
    integer ( kind = 4 ) iargc
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) node_dim
    character ( len = 255 ) node_filename
    integer ( kind = 4 ) node_num
    real ( kind = 8 ), allocatable, dimension ( :, : ) :: node_xy
    character ( len = 255 ) prefix
    character ( len = 255 ) element_filename
    integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_neighbor
    integer ( kind = 4 ), allocatable, dimension ( :, : ) :: triangle_node
    integer ( kind = 4 )  nabes_node(node_num,15)
    integer ( kind = 4 )  nabes_num(node_num)
    integer ( kind = 4 ) triangle_num
    integer ( kind = 4 ) nabes_max
    integer ( kind = 4 ) triangle_order
    integer ( kind = 4 ) :: element_show = 1
    integer ( kind = 4 ) :: node_show = 1
    real(8) x(dim,maxn),itype(maxn),itypee(node_num),hsml
    integer(4) countiac(maxn),numb(maxn,15),itimestep


    node_dim=dim
    allocate ( node_xy(1:node_dim,1:node_num) )

    node_xy(1:node_dim,1:node_num)=x(1:node_dim,1:node_num)

    triangle_order = 3

    allocate ( triangle_node(triangle_order,3*node_num) )
    allocate ( triangle_neighbor(triangle_order,3*node_num) )
    itypee(1:node_num)=itype(1:node_num)

    call dtris2 ( node_num, node_xy,itypee,hsml, triangle_num, triangle_node, &
        triangle_neighbor )

    !if (plot.and.mod(itimestep,save_step).eq.0) call triangulation_order3_plot ( file_name, node_num, node_xy, triangle_num, &
    !    triangle_node, node_show, element_show )

    nabes_max=15*node_num


    call triangulation_order3_neighbor_nodes ( node_num, triangle_num, &
        nabes_max, triangle_node, nabes_num, nabes_node )

    countiac(1:maxn)=0
    countiac(1:node_num)=nabes_num(1:node_num)
    numb(1:maxn,1:15)=0
    numb(1:node_num,1:15)=nabes_node(1:node_num,1:15)
    !
    !  Free memory.
    !
    deallocate ( node_xy )
    deallocate ( triangle_node )
    deallocate ( triangle_neighbor )

    return
    end



    subroutine triangulation_order3_neighbor_nodes ( node_num, element_num, &
        nabes_max, element_node, nabes_num, nabes_node )

    !*****************************************************************************80
    !
    !! TRIANGULATION_ORDER3_NEIGHBOR_NODES determines triangulation neighbor nodes.
    !
    !  Example:
    !
    !    On input, the triangle data structure is:
    !
    !    Triangle  Nodes
    !    --------  ----------
    !     1        3,   4,   1
    !     2        3,   1,   2
    !     3        3,   2,   6
    !     4        2,   1,   5
    !     5        6,   2,   5
    !
    !  On output, the auxilliary neighbor arrays are:
    !
    !    Node  Num  First
    !    ----  ---  -----
    !     1     4     1
    !     2     4     5
    !     3     4     9
    !     4     2    13
    !     5     3    15
    !     6     3    18
    !
    !  and the neighbor array is:
    !
    !    Position  Node
    !    --------  ----
    !
    !     1        2
    !     2        3
    !     3        4
    !     4        5
    !    -----------
    !     5        1
    !     6        3
    !     7        5
    !     8        6
    !    -----------
    !     9        1
    !    10        2
    !    11        4
    !    12        6
    !    -----------
    !    13        1
    !    14        3
    !    -----------
    !    15        1
    !    16        2
    !    17        6
    !    -----------
    !    18        2
    !    19        3
    !    20        5
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 July 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
    !
    !    Input, integer ( kind = 4 ) NABES_MAX, the maximum dimension of NABES.
    !
    !    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), the nodes
    !    that make up each triangle.
    !
    !    Output, integer ( kind = 4 ) NABES_FIRST(NODE_NUM), the index in NABES
    !    of the first neighbor in the list for each node.
    !
    !    Output, integer ( kind = 4 ) NABES_NUM(NODE_NUM), the number of neighbors
    !    of each node.
    !
    !    Output, integer ( kind = 4 ) NABES_DIM, the dimension of NABES.
    !
    !    Output, integer ( kind = 4 ) NABES(NABES_DIM), a list of the neighbors
    !    of all the nodes.  Neighbors of node 1 are listed first, and so on.
    !
    implicit none

    integer ( kind = 4 ) nabes_max
    integer ( kind = 4 ) node_num
    integer ( kind = 4 ) element_num
    integer ( kind = 4 ), parameter :: element_order = 3

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i_current
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) nabe
    integer ( kind = 4 ) nabes(nabes_max)
    integer ( kind = 4 ) nabes1(nabes_max)
    integer ( kind = 4 ) nabes_dim
    integer ( kind = 4 ) nabes_first(node_num)
    integer ( kind = 4 ) nabes_num(node_num)
    integer ( kind = 4 ) tri
    integer ( kind = 4 ) element_node(element_order,element_num)
    integer ( kind = 4 ) unique_num
    integer ( kind = 4 ) nabes_node(node_num,15)

    !
    !  Step 1.  From the triangle list (I,J,K)
    !  construct the neighbor relations: (I,J), (J,K), (K,I), (J,I), (K,J), (I,K).
    !
    nabes_dim = 0
    do tri = 1, element_num
        i = element_node(1,tri)
        j = element_node(2,tri)
        k = element_node(3,tri)
        nabes1(nabes_dim+1:nabes_dim+6) = (/ i, i, j, j, k, k /)
        nabes(nabes_dim+1:nabes_dim+6) = (/ j, k, i, k, i, j /)
        nabes_dim = nabes_dim + 6
    end do
    !
    !  Step 2. Dictionary sort the neighbor relations.
    !
    call i4vec2_sort_a ( nabes_dim, nabes1, nabes )
    !
    !  Step 3. Remove duplicate entries.
    !
    call i4vec2_sorted_unique ( nabes_dim, nabes1, nabes, unique_num )

    nabes_dim = unique_num
    !
    !  Step 4. Construct the NABES_NUM and NABES_FIRST data.
    !
    nabes_num(1:node_num) = 0
    nabes_first(1:node_num) = 0
    nabes_node(1:node_num,1:15) = 0

    i_current = 0
    do nabe = 1, nabes_dim
        i = nabes1(nabe)
        if ( i == i_current ) then
            nabes_num(i) = nabes_num(i) + 1
        else
            i_current = i
            nabes_first(i) = nabe
            nabes_num(i) = 1
        end if
    end do

    do i=1,node_num
        do j=1,nabes_num(i)
            nabes_node(i,j)=nabes(nabes_first(i)+j-1)
        end do
    end do


    return
    end


    subroutine triangulation_order3_plot ( file_name, node_num, node_xy, &
        element_num, element_node, node_show, element_show )

    !*****************************************************************************80
    !
    !! TRIANGULATION_ORDER3_PLOT plots a 3-node triangulation of a set of nodes.
    !
    !  Discussion:
    !
    !    The triangulation is most usually a Delaunay triangulation,
    !    but this is not necessary.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 June 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) FILE_NAME, the name of the output file.
    !
    !    Input, integer ( kind = 4 ) NODE_NUM, the number of nodes.
    !
    !    Input, real ( kind = 8 ) NODE_XY(2,NODE_NUM), the coordinates of the nodes.
    !
    !    Input, integer ( kind = 4 ) ELEMENT_NUM, the number of triangles.
    !
    !    Input, integer ( kind = 4 ) ELEMENT_NODE(3,ELEMENT_NUM), lists, for
    !    each triangle, the indices of the nodes that form the vertices of the
    !    triangle.
    !
    !    Input, integer ( kind = 4 ) NODE_SHOW,
    !    0, do not show nodes;
    !    1, show nodes;
    !    2, show nodes and label them.
    !
    !    Input, integer ( kind = 4 ) ELEMENT_SHOW,
    !    0, do not show triangles;
    !    1, show triangles;
    !    2, show triangles and label them.
    !
    implicit none

    integer ( kind = 4 ) node_num
    integer ( kind = 4 ) element_num
    integer ( kind = 4 ), parameter :: element_order = 3

    real ( kind = 8 ) ave_x
    real ( kind = 8 ) ave_y
    real(8):: circle_size
    real(8) delta
    real(8) e
    character ( len = * )  file_name
    integer ( kind = 4 ) file_unit
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) node
    integer ( kind = 4 ) node_show
    real ( kind = 8 ) node_xy(2,node_num)
    character ( len = 40 ) string
    integer ( kind = 4 ) triangle
    integer ( kind = 4 ) element_node(element_order,element_num)
    integer ( kind = 4 ) element_show
    real ( kind = 8 ) x_max
    real ( kind = 8 ) x_min
    real(8) x_ps
    integer ( kind = 4 ) :: x_ps_max = 576
    integer ( kind = 4 ) :: x_ps_max_clip = 594
    integer ( kind = 4 ) :: x_ps_min = 36
    integer ( kind = 4 ) :: x_ps_min_clip = 18
    real ( kind = 8 ) x_scale
    real ( kind = 8 ) y_max
    real ( kind = 8 ) y_min
    real(8) y_ps
    integer ( kind = 4 ) :: y_ps_max = 666
    integer ( kind = 4 ) :: y_ps_max_clip = 684
    integer ( kind = 4 ) :: y_ps_min = 126
    integer ( kind = 4 ) :: y_ps_min_clip = 108
    real ( kind = 8 ) y_scale
    !
    !  We need to do some figuring here, so that we can determine
    !  the range of the data, and hence the height and width
    !  of the piece of paper.
    !

    x_max = maxval ( node_xy(1,1:node_num) )
    x_min = minval ( node_xy(1,1:node_num) )
    x_scale = x_max - x_min

    x_max = x_max + 0.05D+00 * x_scale
    x_min = x_min - 0.05D+00 * x_scale
    x_scale = x_max - x_min

    y_max = maxval ( node_xy(2,1:node_num) )
    y_min = minval ( node_xy(2,1:node_num) )
    y_scale = y_max - y_min

    y_max = y_max + 0.05D+00 * y_scale
    y_min = y_min - 0.05D+00 * y_scale
    y_scale = y_max - y_min

    if ( x_scale < y_scale ) then

        delta = nint ( real ( x_ps_max - x_ps_min, kind = 8 ) &
            * ( y_scale - x_scale ) / ( 2.0D+00 * y_scale ) )

        x_ps_max = x_ps_max - delta
        x_ps_min = x_ps_min + delta

        x_ps_max_clip = x_ps_max_clip - delta
        x_ps_min_clip = x_ps_min_clip + delta

        x_scale = y_scale

    else if ( y_scale < x_scale ) then

        delta = nint ( real ( y_ps_max - y_ps_min, kind = 8 ) &
            * ( x_scale - y_scale ) / ( 2.0D+00 * x_scale ) )

        y_ps_max      = y_ps_max - delta
        y_ps_min      = y_ps_min + delta

        y_ps_max_clip = y_ps_max_clip - delta
        y_ps_min_clip = y_ps_min_clip + delta

        y_scale = x_scale

    end if

    call get_unit ( file_unit )

    open ( unit = file_unit, file = file_name, status = 'replace', &
        iostat = ios )

    if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TRIANGULATION_ORDER3_PLOT - Fatal error!'
        write ( *, '(a)' ) '  Can not open output file.'
        return
    end if

    write ( file_unit, '(a)' ) '%!PS-Adobe-3.0 EPSF-3.0'
    write ( file_unit, '(a)' ) '%%Creator: triangulation_order3_plot.f90'
    write ( file_unit, '(a)' ) '%%Title: ' // trim ( file_name )
    write ( file_unit, '(a)' ) '%%Pages: 1'
    write ( file_unit, '(a,i3,2x,i3,2x,i3,2x,i3)' ) '%%BoundingBox: ', &
        x_ps_min, y_ps_min, x_ps_max, y_ps_max
    write ( file_unit, '(a)' ) '%%Document-Fonts: Times-Roman'
    write ( file_unit, '(a)' ) '%%LanguageLevel: 1'
    write ( file_unit, '(a)' ) '%%EndComments'
    write ( file_unit, '(a)' ) '%%BeginProlog'
    write ( file_unit, '(a)' ) '/inch {72 mul} def'
    write ( file_unit, '(a)' ) '%%EndProlog'
    write ( file_unit, '(a)' ) '%%Page: 1 1'
    write ( file_unit, '(a)' ) 'save'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Increase line width from default 0.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '1 setlinewidth'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB line color to very light gray.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.900  0.900  0.900 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Draw a gray border around the page.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) 'newpath'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' moveto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_min, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_max, y_ps_max, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_max, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', x_ps_min, y_ps_min, ' lineto'
    write ( file_unit, '(a)' ) 'stroke'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the RGB color to black.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '0.000  0.000  0.000 setrgbcolor'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Set the font and its size.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '/Times-Roman findfont'
    write ( file_unit, '(a)' ) '0.50 inch scalefont'
    write ( file_unit, '(a)' ) 'setfont'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Print a title.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  210  702  moveto'
    write ( file_unit, '(a)' ) '%  (Triangulation)  show'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  Define a clipping polygon.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) 'newpath'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
        x_ps_min_clip, y_ps_min_clip, ' moveto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
        x_ps_max_clip, y_ps_min_clip, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
        x_ps_max_clip, y_ps_max_clip, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
        x_ps_min_clip, y_ps_max_clip, ' lineto'
    write ( file_unit, '(a,i3,2x,i3,2x,a)' ) '  ', &
        x_ps_min_clip, y_ps_min_clip, ' lineto'
    write ( file_unit, '(a)' ) 'clip newpath'
    !
    !  Draw the nodes.
    !
    if ( node_num <= 200 ) then
        circle_size = 0.2
    else if ( node_num <= 500 ) then
        circle_size = 0.2
    else if ( node_num <= 1000 ) then
        circle_size = 0.2
    else if ( node_num <= 5000 ) then
        circle_size = 0.2
    else
        circle_size = 0.2
    end if

    if ( 1 <= node_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw filled dots at the nodes.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.150  0.750 setrgbcolor'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

            x_ps = ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
                + (         node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
                / ( x_max                   - x_min )

            y_ps = ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
                + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
                / ( y_max                   - y_min )

            write ( file_unit, '(a,e11.4,2x,e11.4,2x,E11.4,2x,a)' ) 'newpath ', x_ps, y_ps, &
                circle_size, '0 360 arc closepath fill'

        end do

    end if
    !
    !  Label the nodes.
    !
    if ( 2 <= node_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the nodes:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to darker blue.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.000  0.250  0.850 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.010 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do node = 1, node_num

            x_ps = ( ( x_max - node_xy(1,node)         ) * real ( x_ps_min, kind = 8 )   &
                + (       + node_xy(1,node) - x_min ) * real ( x_ps_max, kind = 8 ) ) &
                / ( x_max                   - x_min )

            y_ps =  ( ( y_max - node_xy(2,node)         ) * real ( y_ps_min, kind = 8 )   &
                + (         node_xy(2,node) - y_min ) * real ( y_ps_max, kind = 8 ) ) &
                / ( y_max                   - y_min )

            write ( string, '(i4)' ) node
            string = adjustl ( string )

            write ( file_unit, '(e11.4,2x,e11.4,a)' ) x_ps+0.01, y_ps+0.01, &
                ' moveto (' // trim ( string ) // ') show'

        end do

    end if
    !
    !  Draw the triangles.
    !
    if ( 1 <= element_show ) then
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.900  0.200  0.100 setrgbcolor'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Draw the triangles.'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, element_num

            write ( file_unit, '(a)' ) 'newpath'

            do i = 1, 4

                e = i4_wrap ( i, 1, 3 )

                node = element_node(e,triangle)

                x_ps = ( ( x_max - node_xy(1,node)         ) &
                    * real ( x_ps_min, kind = 8 )   &
                    + (         node_xy(1,node) - x_min ) &
                    * real ( x_ps_max, kind = 8 ) ) &
                    / ( x_max                   - x_min )

                y_ps = ( ( y_max - node_xy(2,node)         ) &
                    * real ( y_ps_min, kind = 8 )   &
                    + (         node_xy(2,node) - y_min ) &
                    * real ( y_ps_max, kind = 8 ) ) &
                    / ( y_max                   - y_min )

                if ( i == 1 ) then
                    write ( file_unit, '(e11.4,2x,e11.4,2x,a)' ) x_ps, y_ps, ' moveto'
                else
                    write ( file_unit, '(e11.4,2x,e11.4,2x,a)' ) x_ps, y_ps, ' lineto'
                end if

            end do

            write ( file_unit, '(a)' ) 'stroke'

        end do

    end if
    !
    !  Label the triangles.
    !
    if ( 2 <= element_show ) then

        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Label the triangles:'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '%  Set the RGB color to darker red.'
        write ( file_unit, '(a)' ) '%'
        write ( file_unit, '(a)' ) '0.950  0.250  0.150 setrgbcolor'
        write ( file_unit, '(a)' ) '/Times-Roman findfont'
        write ( file_unit, '(a)' ) '0.010 inch scalefont'
        write ( file_unit, '(a)' ) 'setfont'
        write ( file_unit, '(a)' ) '%'

        do triangle = 1, element_num

            ave_x = 0.0D+00
            ave_y = 0.0D+00

            do i = 1, 3

                node = element_node(i,triangle)

                ave_x = ave_x + node_xy(1,node)
                ave_y = ave_y + node_xy(2,node)

            end do

            ave_x = ave_x / 3.0D+00
            ave_y = ave_y / 3.0D+00

            x_ps =  ( ( x_max - ave_x         ) * real ( x_ps_min, kind = 8 )   &
                + (       + ave_x - x_min ) * real ( x_ps_max, kind = 8 ) ) &
                / ( x_max         - x_min )

            y_ps =  ( ( y_max - ave_y         ) * real ( y_ps_min, kind = 8 )   &
                + (         ave_y - y_min ) * real ( y_ps_max, kind = 8 ) ) &
                / ( y_max         - y_min )

            write ( string, '(i4)' ) triangle
            string = adjustl ( string )

            write ( file_unit, '(e11.4,2x,e11.4,a)' ) x_ps, y_ps, ' moveto (' &
                // trim ( string ) // ') show'

        end do

    end if

    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) 'restore  showpage'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%  End of page.'
    write ( file_unit, '(a)' ) '%'
    write ( file_unit, '(a)' ) '%%Trailer'
    write ( file_unit, '(a)' ) '%%EOF'
    close ( unit = file_unit )

    return
    end


    subroutine ch_cap ( c )

    !*****************************************************************************80
    !
    !! CH_CAP capitalizes a single character.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character C, the character to capitalize.
    !
    implicit none

    character c
    integer ( kind = 4 ) itemp

    itemp = ichar ( c )

    if ( 97 <= itemp .and. itemp <= 122 ) then
        c = char ( itemp - 32 )
    end if

    return
    end
    function ch_eqi ( c1, c2 )

    !*****************************************************************************80
    !
    !! CH_EQI is a case insensitive comparison of two characters for equality.
    !
    !  Example:
    !
    !    CH_EQI ( 'A', 'a' ) is .TRUE.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    28 July 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character C1, C2, the characters to compare.
    !
    !    Output, logical CH_EQI, the result of the comparison.
    !
    implicit none

    logical ch_eqi
    character c1
    character c1_cap
    character c2
    character c2_cap

    c1_cap = c1
    c2_cap = c2

    call ch_cap ( c1_cap )
    call ch_cap ( c2_cap )

    if ( c1_cap == c2_cap ) then
        ch_eqi = .true.
    else
        ch_eqi = .false.
    end if

    return
    end
    subroutine ch_to_digit ( c, digit )

    !*****************************************************************************80
    !
    !! CH_TO_DIGIT returns the value of a base 10 digit.
    !
    !  Example:
    !
    !     C   DIGIT
    !    ---  -----
    !    '0'    0
    !    '1'    1
    !    ...  ...
    !    '9'    9
    !    ' '    0
    !    'X'   -1
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    04 August 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character C, the decimal digit, '0' through '9' or blank
    !    are legal.
    !
    !    Output, integer ( kind = 4 ) DIGIT, the corresponding value.
    !    If C was 'illegal', then DIGIT is -1.
    !
    implicit none

    character c
    integer ( kind = 4 ) digit

    if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

    else if ( c == ' ' ) then

        digit = 0

    else

        digit = -1

    end if

    return
    end
    function diaedg ( x0, y0, x1, y1, x2, y2, x3, y3 )

    !*****************************************************************************80
    !
    !! DIAEDG chooses a diagonal edge.
    !
    !  Discussion:
    !
    !    The routine determines whether 0--2 or 1--3 is the diagonal edge
    !    that should be chosen, based on the circumcircle criterion, where
    !    (X0,Y0), (X1,Y1), (X2,Y2), (X3,Y3) are the vertices of a simple
    !    quadrilateral in counterclockwise order.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 February 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X0, Y0, X1, Y1, X2, Y2, X3, Y3, the
    !    coordinates of the vertices of a quadrilateral, given in
    !    counter clockwise order.
    !
    !    Output, integer ( kind = 4 ) DIAEDG, chooses a diagonal:
    !    +1, if diagonal edge 02 is chosen;
    !    -1, if diagonal edge 13 is chosen;
    !     0, if the four vertices are cocircular.
    !
    implicit none

    real ( kind = 8 ) ca
    real ( kind = 8 ) cb
    integer ( kind = 4 ) diaedg
    real ( kind = 8 ) dx10
    real ( kind = 8 ) dx12
    real ( kind = 8 ) dx30
    real ( kind = 8 ) dx32
    real ( kind = 8 ) dy10
    real ( kind = 8 ) dy12
    real ( kind = 8 ) dy30
    real ( kind = 8 ) dy32
    real ( kind = 8 ) s
    real ( kind = 8 ) tol
    real ( kind = 8 ) tola
    real ( kind = 8 ) tolb
    real ( kind = 8 ) x0
    real ( kind = 8 ) x1
    real ( kind = 8 ) x2
    real ( kind = 8 ) x3
    real ( kind = 8 ) y0
    real ( kind = 8 ) y1
    real ( kind = 8 ) y2
    real ( kind = 8 ) y3

    tol = 1.0D-100 * epsilon ( tol )

    dx10 = x1 - x0
    dy10 = y1 - y0
    dx12 = x1 - x2
    dy12 = y1 - y2
    dx30 = x3 - x0
    dy30 = y3 - y0
    dx32 = x3 - x2
    dy32 = y3 - y2

    tola = tol * max ( abs ( dx10 ), abs ( dy10 ), abs ( dx30 ), abs ( dy30 ) )
    tolb = tol * max ( abs ( dx12 ), abs ( dy12 ), abs ( dx32 ), abs ( dy32 ) )

    ca = dx10 * dx30 + dy10 * dy30
    cb = dx12 * dx32 + dy12 * dy32

    if ( tola < ca .and. tolb < cb ) then

        diaedg = -1

    else if ( ca < -tola .and. cb < -tolb ) then

        diaedg = 1

    else

        tola = max ( tola, tolb )
        s = ( dx10 * dy30 - dx30 * dy10 ) * cb + ( dx32 * dy12 - dx12 * dy32 ) * ca

        if ( tola < s ) then
            diaedg = -1
        else if ( s < -tola ) then
            diaedg = 1
        else
            diaedg = 0
        end if

    end if

    return
    end
    subroutine dtris2 ( point_num, point_xy, itype,hsml, tri_num, tri_vert, tri_nabe )

    !*****************************************************************************80
    !
    !! DTRIS2 constructs a Delaunay triangulation of 2D vertices.
    !
    !  Discussion:
    !
    !    The routine constructs the Delaunay triangulation of a set of 2D vertices
    !    using an incremental approach and diagonal edge swaps.  Vertices are
    !    first sorted in lexicographically increasing (X,Y) order, and
    !    then are inserted one at a time from outside the convex hull.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 August 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of vertices.
    !
    !    Input/output, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
    !    of the vertices.  On output, the vertices have been sorted into
    !    dictionary order.
    !
    !    Output, integer ( kind = 4 ) TRI_NUM, the number of triangles in the
    !    triangulation; TRI_NUM is equal to 2*POINT_NUM - NB - 2, where NB is the
    !    number of boundary vertices.
    !
    !    Output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the nodes that make up
    !    each triangle.  The elements are indices of POINT_XY.  The vertices of the
    !    triangles are in counter clockwise order.
    !
    !    Output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor
    !    list.  Positive elements are indices of TIL; negative elements are used
    !    for links of a counter clockwise linked list of boundary edges;
    !    LINK = -(3*I + J-1) where I, J = triangle, edge index; TRI_NABE(J,I) refers
    !    to the neighbor along edge from vertex J to J+1 (mod 3).
    !
    use config_parameter
    implicit none

    integer ( kind = 4 ) point_num

    real ( kind = 8 ) cmax
    integer ( kind = 4 ) e
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) indx(point_num)
    integer ( kind = 4 ) j
    integer ( kind = 4 ) k
    integer ( kind = 4 ) l
    integer ( kind = 4 ) ledg
    integer ( kind = 4 ) lr
    integer ( kind = 4 ) lrline
    integer ( kind = 4 ) ltri
    integer ( kind = 4 ) m
    integer ( kind = 4 ) m1
    integer ( kind = 4 ) m2
    integer ( kind = 4 ) n
    real ( kind = 8 ) point_xy(2,point_num)
    integer ( kind = 4 ) redg
    integer ( kind = 4 ) rtri
    integer ( kind = 4 ) stack(point_num)
    integer ( kind = 4 ) t
    real ( kind = 8 ) tol
    integer ( kind = 4 ) top
    integer ( kind = 4 ) tri_nabe(3,point_num*2)
    integer ( kind = 4 ) tri_num
    integer ( kind = 4 ) tri_vert(3,point_num*2)
    real(8) mean_len, a, b, c, angle_1, angle_2, angle_3, min_angle,max_len,itype(point_num),hsml,min_len,sign12,sign23,sign31

    tol = 1.0D-100 * epsilon ( tol )

    ierr = 0
    !
    !  Sort the vertices by increasing (x,y).
    !
    call r82vec_sort_heap_index_a ( point_num, point_xy, indx )

    call r82vec_permute ( point_num, indx, point_xy )
    !
    !  Make sure that the data points are "reasonably" distinct.
    !
    m1 = 1

    do i = 2, point_num

        m = m1
        m1 = i

        k = 0

        do j = 1, 2

            cmax = max ( abs ( point_xy(j,m) ), abs ( point_xy(j,m1) ) )

            if ( tol * ( cmax + 1.0D+00 ) &
                < abs ( point_xy(j,m) - point_xy(j,m1) ) ) then
            k = j
            exit
            end if

        end do

        if ( k == 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
            write ( *, '(a,i8)' ) '  Fails for point number I = ', i
            write ( *, '(a,i8)' ) '  M = ', m
            write ( *, '(a,i8)' ) '  M1 = ', m1
            write ( *, '(a,2g14.6)' ) '  X,Y(M)  = ', point_xy(1,m), point_xy(2,m)
            write ( *, '(a,2g14.6)' ) '  X,Y(M1) = ', point_xy(1,m1), point_xy(2,m1)
            ierr = 224
            return
        end if

    end do
    !
    !  Starting from points M1 and M2, search for a third point M that
    !  makes a "healthy" triangle (M1,M2,M)
    !
    m1 = 1
    m2 = 2
    j = 3

    do

        if ( point_num < j ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
            ierr = 225
            return
        end if

        m = j

        lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
            point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

        if ( lr /= 0 ) then
            exit
        end if

        j = j + 1

    end do
    !
    !  Set up the triangle information for (M1,M2,M), and for any other
    !  triangles you created because points were collinear with M1, M2.
    !
    tri_num = j - 2

    if ( lr == -1 ) then

        tri_vert(1,1) = m1
        tri_vert(2,1) = m2
        tri_vert(3,1) = m
        tri_nabe(3,1) = -3

        do i = 2, tri_num

            m1 = m2
            m2 = i+1
            tri_vert(1,i) = m1
            tri_vert(2,i) = m2
            tri_vert(3,i) = m
            tri_nabe(1,i-1) = -3 * i
            tri_nabe(2,i-1) = i
            tri_nabe(3,i) = i - 1

        end do

        tri_nabe(1,tri_num) = -3 * tri_num - 1
        tri_nabe(2,tri_num) = -5
        ledg = 2
        ltri = tri_num

    else

        tri_vert(1,1) = m2
        tri_vert(2,1) = m1
        tri_vert(3,1) = m
        tri_nabe(1,1) = -4

        do i = 2, tri_num
            m1 = m2
            m2 = i+1
            tri_vert(1,i) = m2
            tri_vert(2,i) = m1
            tri_vert(3,i) = m
            tri_nabe(3,i-1) = i
            tri_nabe(1,i) = -3 * i - 3
            tri_nabe(2,i) = i - 1
        end do

        tri_nabe(3,tri_num) = -3 * tri_num
        tri_nabe(2,1) = -3 * tri_num - 2
        ledg = 2
        ltri = 1

    end if
    !
    !  Insert the vertices one at a time from outside the convex hull,
    !  determine visible boundary edges, and apply diagonal edge swaps until
    !  Delaunay triangulation of vertices (so far) is obtained.
    !
    top = 0

    do i = j+1, point_num

        m = i
        m1 = tri_vert(ledg,ltri)

        if ( ledg <= 2 ) then
            m2 = tri_vert(ledg+1,ltri)
        else
            m2 = tri_vert(1,ltri)
        end if

        lr = lrline ( point_xy(1,m), point_xy(2,m), point_xy(1,m1), &
            point_xy(2,m1), point_xy(1,m2), point_xy(2,m2), 0.0D+00 )

        if ( 0 < lr ) then
            rtri = ltri
            redg = ledg
            ltri = 0
        else
            l = -tri_nabe(ledg,ltri)
            rtri = l / 3
            redg = mod(l,3) + 1
        end if

        call vbedg ( point_xy(1,m), point_xy(2,m), point_num, point_xy, tri_num, &
            tri_vert, tri_nabe, ltri, ledg, rtri, redg )

        n = tri_num + 1
        l = -tri_nabe(ledg,ltri)

        do

            t = l / 3
            e = mod ( l, 3 ) + 1
            l = -tri_nabe(e,t)
            m2 = tri_vert(e,t)

            if ( e <= 2 ) then
                m1 = tri_vert(e+1,t)
            else
                m1 = tri_vert(1,t)
            end if

            tri_num = tri_num + 1
            tri_nabe(e,t) = tri_num
            tri_vert(1,tri_num) = m1
            tri_vert(2,tri_num) = m2
            tri_vert(3,tri_num) = m
            tri_nabe(1,tri_num) = t
            tri_nabe(2,tri_num) = tri_num - 1
            tri_nabe(3,tri_num) = tri_num + 1
            top = top + 1

            if ( point_num < top ) then
                ierr = 8
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
                write ( *, '(a)' ) '  Stack overflow.'
                return
            end if

            stack(top) = tri_num

            if ( t == rtri .and. e == redg ) then
                exit
            end if

        end do

        tri_nabe(ledg,ltri) = -3 * n - 1
        tri_nabe(2,n) = -3 * tri_num - 2
        tri_nabe(3,tri_num) = -l
        ltri = n
        ledg = 2

        call swapec ( m, top, ltri, ledg, point_num, point_xy, tri_num, &
            tri_vert, tri_nabe, stack, ierr )

        if ( ierr /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'DTRIS2 - Fatal error!'
            write ( *, '(a)' ) '  Error return from SWAPEC.'
            return
        end if

    end do
    !
    !  Now account for the sorting that we did.
    !
    do i = 1, 3
        do j = 1, tri_num
            tri_vert(i,j) = indx ( tri_vert(i,j) )
        end do
    end do

    call perm_inverse ( point_num, indx )

    call r82vec_permute ( point_num, indx, point_xy )


    mean_len=0.0;

    do i=1,tri_num
101     a=sqrt((point_xy(1,tri_vert(1,i))-point_xy(1,tri_vert(2,i)))**2+&
            (point_xy(2,tri_vert(1,i))-point_xy(2,tri_vert(2,i)))**2)
        b=sqrt((point_xy(1,tri_vert(2,i))-point_xy(1,tri_vert(3,i)))**2+&
            (point_xy(2,tri_vert(2,i))-point_xy(2,tri_vert(3,i)))**2)
        c=sqrt((point_xy(1,tri_vert(3,i))-point_xy(1,tri_vert(1,i)))**2+&
            (point_xy(2,tri_vert(3,i))-point_xy(2,tri_vert(1,i)))**2)

        angle_1=acos((a**2+b**2-c**2)/2.0/a/b)
        angle_2=acos((a**2+c**2-b**2)/2.0/a/c)
        angle_3=acos((c**2+b**2-a**2)/2.0/c/b)
        min_angle=min(angle_1,angle_2,angle_3)

        if (min_angle.gt.12.0*3.1415926/180.0) then
            mean_len=mean_len+a+b+c
        else
            if (i.lt.tri_num) then
                do j=i,tri_num-1
                    tri_vert(1,j)=tri_vert(1,j+1)
                    tri_vert(2,j)=tri_vert(2,j+1)
                    tri_vert(3,j)=tri_vert(3,j+1)
                end do

                tri_vert(1,tri_num)=0
                tri_vert(2,tri_num)=0
                tri_vert(3,tri_num)=0
                tri_num=tri_num-1
                goto 101
            elseif (i.eq.tri_num) then
                tri_vert(1,tri_num)=0
                tri_vert(2,tri_num)=0
                tri_vert(3,tri_num)=0
                tri_num=tri_num-1
                goto 300
            endif
        endif
    end do
300 mean_len=mean_len/(3.0*tri_num)

    do i=1,tri_num
601     sign12=itype(tri_vert(1,i))*itype(tri_vert(2,i))
        sign31=itype(tri_vert(1,i))*itype(tri_vert(3,i))
        sign23=itype(tri_vert(2,i))*itype(tri_vert(3,i))
        if ((sign12.lt.0.0).or.(sign23.lt.0.0).or.(sign31.lt.0.0)) then
            a=1000.0
            b=1000.0
            c=1000.0

            if (sign12.lt.0.0) a=sqrt((point_xy(1,tri_vert(1,i))-point_xy(1,tri_vert(2,i)))**2+(point_xy(2,tri_vert(1,i))-point_xy(2,tri_vert(2,i)))**2)
            if (sign31.lt.0.0) c=sqrt((point_xy(1,tri_vert(1,i))-point_xy(1,tri_vert(3,i)))**2+(point_xy(2,tri_vert(1,i))-point_xy(2,tri_vert(3,i)))**2)
            if (sign23.lt.0.0) b=sqrt((point_xy(1,tri_vert(3,i))-point_xy(1,tri_vert(2,i)))**2+(point_xy(2,tri_vert(3,i))-point_xy(2,tri_vert(2,i)))**2)

            min_len=min(a,b,c)

            if (min_len.gt.hsml) goto 600

        endif

        goto 700

600     if (i.lt.tri_num) then
            do j=i,tri_num-1
                tri_vert(1,j)=tri_vert(1,j+1)
                tri_vert(2,j)=tri_vert(2,j+1)
                tri_vert(3,j)=tri_vert(3,j+1)
            end do

            tri_vert(1,tri_num)=0
            tri_vert(2,tri_num)=0
            tri_vert(3,tri_num)=0
            tri_num=tri_num-1
            goto 601
        elseif (i.eq.tri_num) then
            tri_vert(1,tri_num)=0
            tri_vert(2,tri_num)=0
            tri_vert(3,tri_num)=0
            tri_num=tri_num-1
            exit
        endif


700 enddo


    do i=1,tri_num
200     a=sqrt((point_xy(1,tri_vert(1,i))-point_xy(1,tri_vert(2,i)))**2+&
            (point_xy(2,tri_vert(1,i))-point_xy(2,tri_vert(2,i)))**2)
        b=sqrt((point_xy(1,tri_vert(2,i))-point_xy(1,tri_vert(3,i)))**2+&
            (point_xy(2,tri_vert(2,i))-point_xy(2,tri_vert(3,i)))**2)
        c=sqrt((point_xy(1,tri_vert(3,i))-point_xy(1,tri_vert(1,i)))**2+&
            (point_xy(2,tri_vert(3,i))-point_xy(2,tri_vert(1,i)))**2)
        max_len=max (a,b,c)
        if (max_len>=2.0*hsml) then
            if (i.lt.tri_num) then
                do j=i,tri_num-1
                    tri_vert(1,j)=tri_vert(1,j+1)
                    tri_vert(2,j)=tri_vert(2,j+1)
                    tri_vert(3,j)=tri_vert(3,j+1)
                end do
                tri_vert(1,tri_num)=0
                tri_vert(2,tri_num)=0
                tri_vert(3,tri_num)=0
                tri_num=tri_num-1
                goto 200
            elseif (i.eq.tri_num) then
                tri_vert(1,tri_num)=0
                tri_vert(2,tri_num)=0
                tri_vert(3,tri_num)=0
                tri_num=tri_num-1
                exit
            endif
        endif
    end do

    return
    end
    subroutine file_column_count ( input_file_name, column_num )

    !*****************************************************************************80
    !
    !! FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
    !
    !  Discussion:
    !
    !    The file is assumed to be a simple text file.
    !
    !    Most lines of the file is presumed to consist of COLUMN_NUM words,
    !    separated by spaces.  There may also be some blank lines, and some
    !    comment lines,
    !    which have a "#" in column 1.
    !
    !    The routine tries to find the first non-comment non-blank line and
    !    counts the number of words in that line.
    !
    !    If all lines are blanks or comments, it goes back and tries to analyze
    !    a comment line.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    21 June 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) INPUT_FILE_NAME, the name of the file.
    !
    !    Output, integer ( kind = 4 ) COLUMN_NUM, the number of columns in
    !    the file.
    !
    implicit none

    integer ( kind = 4 ) column_num
    logical got_one
    character ( len = * ) input_file_name
    integer ( kind = 4 ) input_status
    integer ( kind = 4 ) input_unit
    character ( len = 255 ) line
    !
    !  Open the file.
    !
    call get_unit ( input_unit )

    open ( unit = input_unit, file = input_file_name, status = 'old', &
        form = 'formatted', access = 'sequential', iostat = input_status )

    if ( input_status /= 0 ) then
        column_num = -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not open the input file "' &
            // trim ( input_file_name ) // '" on unit ', input_unit
        return
    end if
    !
    !  Read one line, but skip blank lines and comment lines.
    !
    got_one = .false.

    do

        read ( input_unit, '(a)', iostat = input_status ) line

        if ( input_status /= 0 ) then
            exit
        end if

        if ( len_trim ( line ) == 0 ) then
            cycle
        end if

        if ( line(1:1) == '#' ) then
            cycle
        end if

        got_one = .true.
        exit

    end do

    if ( .not. got_one ) then

        rewind ( input_unit )

        do

            read ( input_unit, '(a)', iostat = input_status ) line

            if ( input_status /= 0 ) then
                exit
            end if

            if ( len_trim ( line ) == 0 ) then
                cycle
            end if

            got_one = .true.
            exit

        end do

    end if

    close ( unit = input_unit )

    if ( .not. got_one ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning!'
        write ( *, '(a)' ) '  The file does not seem to contain any data.'
        column_num = -1
        return
    end if

    call s_word_count ( line, column_num )

    return
    end
    subroutine file_row_count ( input_file_name, row_num )

    !*****************************************************************************80
    !
    !! FILE_ROW_COUNT counts the number of row records in a file.
    !
    !  Discussion:
    !
    !    It does not count lines that are blank, or that begin with a
    !    comment symbol '#'.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    06 March 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) INPUT_FILE_NAME, the name of the input file.
    !
    !    Output, integer ( kind = 4 ) ROW_NUM, the number of rows found.
    !
    implicit none

    integer ( kind = 4 ) bad_num
    integer ( kind = 4 ) comment_num
    integer ( kind = 4 ) ierror
    character ( len = * ) input_file_name
    integer ( kind = 4 ) input_status
    integer ( kind = 4 ) input_unit
    character ( len = 255 ) line
    integer ( kind = 4 ) record_num
    integer ( kind = 4 ) row_num

    call get_unit ( input_unit )

    open ( unit = input_unit, file = input_file_name, status = 'old', &
        iostat = input_status )

    if ( input_status /= 0 ) then
        row_num = -1;
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_ROW_COUNT - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not open the input file "' // &
            trim ( input_file_name ) // '" on unit ', input_unit
        stop
    end if

    comment_num = 0
    row_num = 0
    record_num = 0
    bad_num = 0

    do

        read ( input_unit, '(a)', iostat = input_status ) line

        if ( input_status /= 0 ) then
            ierror = record_num
            exit
        end if

        record_num = record_num + 1

        if ( line(1:1) == '#' ) then
            comment_num = comment_num + 1
            cycle
        end if

        if ( len_trim ( line ) == 0 ) then
            comment_num = comment_num + 1
            cycle
        end if

        row_num = row_num + 1

    end do

    close ( unit = input_unit )

    return
    end
    subroutine get_unit ( iunit )

    !*****************************************************************************80
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5, 6 and 9, which
    !    are commonly reserved for console I/O).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer ( kind = 4 ) IUNIT, the free unit number.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) ios
    integer ( kind = 4 ) iunit
    logical lopen

    iunit = 0

    do i = 1, 99

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

            inquire ( unit = i, opened = lopen, iostat = ios )

            if ( ios == 0 ) then
                if ( .not. lopen ) then
                    iunit = i
                    return
                end if
            end if

        end if

    end do

    return
    end
    function i4_modp ( i, j )

    !*****************************************************************************80
    !
    !! I4_MODP returns the nonnegative remainder of I4 division.
    !
    !  Discussion:
    !
    !    If
    !      NREM = I4_MODP ( I, J )
    !      NMULT = ( I - NREM ) / J
    !    then
    !      I = J * NMULT + NREM
    !    where NREM is always nonnegative.
    !
    !    The MOD function computes a result with the same sign as the
    !    quantity being divided.  Thus, suppose you had an angle A,
    !    and you wanted to ensure that it was between 0 and 360.
    !    Then mod(A,360) would do, if A was positive, but if A
    !    was negative, your result would be between -360 and 0.
    !
    !    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !  Example:
    !
    !        I     J     MOD I4_MODP    Factorization
    !
    !      107    50       7       7    107 =  2 *  50 + 7
    !      107   -50       7       7    107 = -2 * -50 + 7
    !     -107    50      -7      43   -107 = -3 *  50 + 43
    !     -107   -50      -7      43   -107 =  3 * -50 + 43
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the number to be divided.
    !
    !    Input, integer ( kind = 4 ) J, the number that divides I.
    !
    !    Output, integer ( kind = 4 ) I4_MODP, the nonnegative remainder when I is
    !    divided by J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ) i4_modp
    integer ( kind = 4 ) j
    integer ( kind = 4 ) value

    if ( j == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
    end if

    value = mod ( i, j )

    if ( value < 0 ) then
        value = value + abs ( j )
    end if

    i4_modp = value

    return
    end
    function i4_sign ( x )

    !*****************************************************************************80
    !
    !! I4_SIGN evaluates the sign of an I4.
    !
    !  Discussion:
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    27 March 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) X, the number whose sign is desired.
    !
    !    Output, integer ( kind = 4 ) I4_SIGN, the sign of X:
    !
    implicit none

    integer ( kind = 4 ) i4_sign
    integer ( kind = 4 ) x

    if ( x < 0 ) then
        i4_sign = -1
    else
        i4_sign = +1
    end if

    return
    end
    function i4_wrap ( ival, ilo, ihi )

    !*****************************************************************************80
    !
    !! I4_WRAP forces an I4 to lie between given limits by wrapping.
    !
    !  Discussion:
    !
    !    An I4 is an integer ( kind = 4 ) value.
    !
    !  Example:
    !
    !    ILO = 4, IHI = 8
    !
    !    I  Value
    !
    !    -2     8
    !    -1     4
    !     0     5
    !     1     6
    !     2     7
    !     3     8
    !     4     4
    !     5     5
    !     6     6
    !     7     7
    !     8     8
    !     9     4
    !    10     5
    !    11     6
    !    12     7
    !    13     8
    !    14     4
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    19 August 2003
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) IVAL, a value.
    !
    !    Input, integer ( kind = 4 ) ILO, IHI, the desired bounds.
    !
    !    Output, integer ( kind = 4 ) I4_WRAP, a "wrapped" version of the value.
    !
    implicit none

    integer ( kind = 4 ) i4_modp
    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) ival
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    integer ( kind = 4 ) value
    integer ( kind = 4 ) wide

    jlo = min ( ilo, ihi )
    jhi = max ( ilo, ihi )

    wide = jhi - jlo + 1

    if ( wide == 1 ) then
        value = jlo
    else
        value = jlo + i4_modp ( ival - jlo, wide )
    end if

    i4_wrap = value

    return
    end
    subroutine i4mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! I4MAT_TRANSPOSE_PRINT_SOME prints some of the transpose of an I4MAT.
    !
    !  Discussion:
    !
    !    An I4MAT is a rectangular array of I4 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 February 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    Input, integer ( kind = 4 ) A(M,N), an M by N matrix to be printed.
    !
    !    Input, integer ILO, JLO, the first row and column to print.
    !
    !    Input, integer IHI, JHI, the last row and column to print.
    !
    !    Input, character ( len = * ) TITLE, an optional title.
    !
    implicit none

    integer ( kind = 4 ), parameter :: incx = 10
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) a(m,n)
    character ( len = 8 ) ctemp(incx)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i2hi
    integer ( kind = 4 ) i2lo
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) inc
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j2hi
    integer ( kind = 4 ) j2lo
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    character ( len = * ) title

    if ( 0 < len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
    end if

    do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
            i2 = i + 1 - i2lo
            write ( ctemp(i2), '(i8)' ) i
        end do

        write ( *, '(''  Row '',10a8)' ) ctemp(1:inc)
        write ( *, '(a)' ) '  Col'
        write ( *, '(a)' ) ' '

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

            do i2 = 1, inc

                i = i2lo - 1 + i2

                write ( ctemp(i2), '(i8)' ) a(i,j)

            end do

            write ( *, '(i5,1x,10a8)' ) j, ( ctemp(i), i = 1, inc )

        end do

    end do

    return
    end
    subroutine i4mat_write ( output_filename, m, n, table )

    !*****************************************************************************80
    !
    !! I4MAT_WRITE writes an I4MAT file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 August 2009
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) OUTPUT_FILENAME, the output file name.
    !
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    Input, integer ( kind = 4 ) N, the number of points.
    !
    !    Input, integer ( kind = 4 ) TABLE(M,N), the table data.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) j
    character ( len = * ) output_filename
    integer ( kind = 4 ) output_status
    integer ( kind = 4 ) output_unit
    character ( len = 30 ) string
    integer ( kind = 4 ) table(m,n)
    !
    !  Open the file.
    !
    call get_unit ( output_unit )

    open ( unit = output_unit, file = output_filename, &
        status = 'replace', iostat = output_status )

    if ( output_status /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4MAT_WRITE - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not open the output file "' // &
            trim ( output_filename ) // '" on unit ', output_unit
        output_unit = -1
        stop
    end if
    !
    !  Create a format string.
    !
    if ( 0 < m .and. 0 < n ) then

        write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
        !
        !  Write the data.
        !
        do j = 1, n
            write ( output_unit, string ) table(1:m,j)
        end do

    end if
    !
    !  Close the file.
    !
    close ( unit = output_unit )

    return
    end
    subroutine i4vec_indicator ( n, a )

    !*****************************************************************************80
    !
    !! I4VEC_INDICATOR sets an I4VEC to the indicator vector.
    !
    !  Discussion:
    !
    !    An I4VEC is a vector of I4's.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    01 May 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of elements of A.
    !
    !    Output, integer ( kind = 4 ) A(N), the array to be initialized.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) a(n)
    integer ( kind = 4 ) i

    do i = 1, n
        a(i) = i
    end do

    return
    end
    function lrline ( xu, yu, xv1, yv1, xv2, yv2, dv )

    !*****************************************************************************80
    !
    !! LRLINE determines if a point is left of, right or, or on a directed line.
    !
    !  Discussion:
    !
    !    The directed line is parallel to, and at a signed distance DV from
    !    a directed base line from (XV1,YV1) to (XV2,YV2).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 July 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) XU, YU, the coordinates of the point whose
    !    position relative to the directed line is to be determined.
    !
    !    Input, real ( kind = 8 ) XV1, YV1, XV2, YV2, the coordinates of two points
    !    that determine the directed base line.
    !
    !    Input, real ( kind = 8 ) DV, the signed distance of the directed line
    !    from the directed base line through the points (XV1,YV1) and (XV2,YV2).
    !    DV is positive for a line to the left of the base line.
    !
    !    Output, integer ( kind = 4 ) LRLINE, the result:
    !    +1, the point is to the right of the directed line;
    !     0, the point is on the directed line;
    !    -1, the point is to the left of the directed line.
    !
    implicit none

    real ( kind = 8 ) dv
    real ( kind = 8 ) dx
    real ( kind = 8 ) dxu
    real ( kind = 8 ) dy
    real ( kind = 8 ) dyu
    integer ( kind = 4 ) lrline
    real ( kind = 8 ) t
    real ( kind = 8 ) tol
    real ( kind = 8 ) tolabs
    real ( kind = 8 ) xu
    real ( kind = 8 ) xv1
    real ( kind = 8 ) xv2
    real ( kind = 8 ) yu
    real ( kind = 8 ) yv1
    real ( kind = 8 ) yv2

    tol = 1.0D-100 * epsilon ( tol )

    dx = xv2 - xv1
    dy = yv2 - yv1
    dxu = xu - xv1
    dyu = yu - yv1

    tolabs = tol * max ( abs ( dx ), abs ( dy ), abs ( dxu ), &
        abs ( dyu ), abs ( dv ) )

    t = dy * dxu - dx * dyu + dv * sqrt ( dx * dx + dy * dy )

    if ( tolabs < t ) then
        lrline = 1
    else if ( -tolabs <= t ) then
        lrline = 0
    else
        lrline = -1
    end if

    return
    end
    subroutine perm_check ( n, p, base, ierror )

    !*****************************************************************************80
    !
    !! PERM_CHECK checks that a vector represents a permutation.
    !
    !  Discussion:
    !
    !    The routine verifies that each of the integers from BASE to
    !    to BASE+N-1 occurs among the N entries of the permutation.
    !
    !    Set the input quantity BASE to 0, if P is a 0-based permutation,
    !    or to 1 if P is a 1-based permutation.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    31 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries.
    !
    !    Input, integer ( kind = 4 ) P(N), the array to check.
    !
    !    Input, integer ( kind = 4 ) BASE, the index base.
    !
    !    Output, integer ( kind = 4 ) IERROR, error flag.
    !    0, the array represents a permutation.
    !    nonzero, the array does not represent a permutation.  The smallest
    !    missing value is equal to IERROR.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) base
    integer ( kind = 4 ) find
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) p(n)
    integer ( kind = 4 ) seek

    ierror = 0

    do seek = base, base + n - 1

        ierror = 1

        do find = 1, n
            if ( p(find) == seek ) then
                ierror = 0
                exit
            end if
        end do

        if ( ierror /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PERM_CHECK - Fatal error!'
            write ( *, '(a)' ) '  The input array does not represent'
            write ( *, '(a)' ) '  a proper permutation.'
            stop
        end if

    end do

    return
    end
    subroutine perm_inverse ( n, p )

    !*****************************************************************************80
    !
    !! PERM_INVERSE inverts a permutation "in place".
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    02 January 2006
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of objects being permuted.
    !
    !    Input/output, integer ( kind = 4 ) P(N), the permutation, in standard
    !    index form.  On output, P describes the inverse permutation
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ), parameter :: base = 1
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i0
    integer ( kind = 4 ) i1
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i4_sign
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) is
    integer ( kind = 4 ) p(n)

    if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of N = ', n
        stop
    end if

    call perm_check ( n, p, base, ierror )

    if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PERM_INVERSE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
    end if


    is = 1

    do i = 1, n

        i1 = p(i)

        do while ( i < i1 )
            i2 = p(i1)
            p(i1) = -i2
            i1 = i2
        end do

        is = - i4_sign ( p(i) )
        p(i) = is * abs ( p(i) )

    end do

    do i = 1, n

        i1 = - p(i)

        if ( 0 <= i1 ) then

            i0 = i

            do

                i2 = p(i1)
                p(i1) = i0

                if ( i2 < 0 ) then
                    exit
                end if

                i0 = i1
                i1 = i2

            end do

        end if

    end do

    return
    end
    subroutine r82vec_permute ( n, p, a )

    !*****************************************************************************80
    !
    !! R82VEC_PERMUTE permutes an R82VEC in place.
    !
    !  Discussion:
    !
    !    An R82VEC is an array of pairs of R8 values.
    !
    !    The same logic can be used to permute an array of objects of any
    !    arithmetic type, or an array of objects of any complexity.  The only
    !    temporary storage required is enough to store a single object.  The number
    !    of data movements made is N + the number of cycles of order 2 or more,
    !    which is never more than N + N/2.
    !
    !  Example:
    !
    !    Input:
    !
    !      N = 5
    !      P = (   2,    4,    5,    1,    3 )
    !      A = ( 1.0,  2.0,  3.0,  4.0,  5.0 )
    !          (11.0, 22.0, 33.0, 44.0, 55.0 )
    !
    !    Output:
    !
    !      A    = (  2.0,  4.0,  5.0,  1.0,  3.0 )
    !             ( 22.0, 44.0, 55.0, 11.0, 33.0 ).
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    13 March 2005
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of objects.
    !
    !    Input, integer ( kind = 4 ) P(N), the permutation.  P(I) = J means
    !    that the I-th element of the output array should be the J-th
    !    element of the input array.
    !
    !    Input/output, real ( kind = 8 ) A(2,N), the array to be permuted.
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ), parameter :: dim_num = 2

    real ( kind = 8 ) a(dim_num,n)
    real ( kind = 8 ) a_temp(dim_num)
    integer ( kind = 4 ), parameter :: base = 1
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) iget
    integer ( kind = 4 ) iput
    integer ( kind = 4 ) istart
    integer ( kind = 4 ) p(n)

    call perm_check ( n, p, base, ierror )

    if ( ierror /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
        write ( *, '(a)' ) '  PERM_CHECK rejects this permutation.'
        stop
    end if
    !
    !  Search for the next element of the permutation that has not been used.
    !
    do istart = 1, n

        if ( p(istart) < 0 ) then

            cycle

        else if ( p(istart) == istart ) then

            p(istart) = - p(istart)
            cycle

        else

            a_temp(1:dim_num) = a(1:dim_num,istart)
            iget = istart
            !
            !  Copy the new value into the vacated entry.
            !
            do

                iput = iget
                iget = p(iget)

                p(iput) = - p(iput)

                if ( iget < 1 .or. n < iget ) then
                    write ( *, '(a)' ) ' '
                    write ( *, '(a)' ) 'R82VEC_PERMUTE - Fatal error!'
                    write ( *, '(a)' ) '  A permutation index is out of range.'
                    write ( *, '(a,i8,a,i8)' ) '  P(', iput, ') = ', iget
                    stop
                end if

                if ( iget == istart ) then
                    a(1:dim_num,iput) = a_temp(1:dim_num)
                    exit
                end if

                a(1:dim_num,iput) = a(1:dim_num,iget)

            end do

        end if

    end do
    !
    !  Restore the signs of the entries.
    !
    p(1:n) = - p(1:n)

    return
    end
    subroutine r82vec_sort_heap_index_a ( n, a, indx )

    !*****************************************************************************80
    !
    !! R82VEC_SORT_HEAP_INDEX_A ascending index heaps an R82VEC.
    !
    !  Discussion:
    !
    !    An R82VEC is an array of R82's.
    !
    !    The sorting is not actually carried out.  Rather an index array is
    !    created which defines the sorting.  This array may be used to sort
    !    or index the array, or to sort or index related arrays keyed on the
    !    original array.
    !
    !    Once the index array is computed, the sorting can be carried out
    !    "implicitly:
    !
    !      A(1:2,INDX(1:N)) is sorted,
    !
    !    or explicitly, by the call
    !
    !      call r82vec_permute ( n, indx, a )
    !
    !    after which A(1:2,I), I = 1 to N is sorted.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    08 December 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of entries in the array.
    !
    !    Input, real ( kind = 8 ) A(2,N), an array to be index-sorted.
    !
    !    Output, integer ( kind = 4 ) INDX(N), the sort index.  The
    !    I-th element of the sorted array is A(1:2,INDX(I)).
    !
    implicit none

    integer ( kind = 4 ) n
    integer ( kind = 4 ), parameter :: dim_num = 2

    real ( kind = 8 ) a(dim_num,n)
    real ( kind = 8 ) aval(dim_num)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) indx(n)
    integer ( kind = 4 ) indxt
    integer ( kind = 4 ) ir
    integer ( kind = 4 ) j
    integer ( kind = 4 ) l

    if ( n < 1 ) then
        return
    end if

    do i = 1, n
        indx(i) = i
    end do

    if ( n == 1 ) then
        return
    end if

    l = n / 2 + 1
    ir = n

    do

        if ( 1 < l ) then

            l = l - 1
            indxt = indx(l)
            aval(1:dim_num) = a(1:dim_num,indxt)

        else

            indxt = indx(ir)
            aval(1:dim_num) = a(1:dim_num,indxt)
            indx(ir) = indx(1)
            ir = ir - 1

            if ( ir == 1 ) then
                indx(1) = indxt
                exit
            end if

        end if

        i = l
        j = l + l

        do while ( j <= ir )

            if ( j < ir ) then
                if (   a(1,indx(j)) <  a(1,indx(j+1)) .or. &
                    ( a(1,indx(j)) == a(1,indx(j+1)) .and. &
                    a(2,indx(j)) <  a(2,indx(j+1)) ) ) then
                j = j + 1
                end if
            end if

            if (   aval(1) <  a(1,indx(j)) .or. &
                ( aval(1) == a(1,indx(j)) .and. &
                aval(2) <  a(2,indx(j)) ) ) then
            indx(i) = indx(j)
            i = j
            j = j + j
            else
                j = ir + 1
            end if

        end do

        indx(i) = indxt

    end do

    return
    end
    subroutine r8mat_data_read ( input_filename, m, n, table )

    !*****************************************************************************80
    !
    !! R8MAT_DATA_READ reads data from an R8MAT file.
    !
    !  Discussion:
    !
    !    The file may contain more than N points, but this routine will
    !    return after reading N of them.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 October 2008
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
    !
    !    Input, integer ( kind = 4 ) M, the spatial dimension.
    !
    !    Input, integer ( kind = 4 ) N, the number of points.
    !
    !    Output, real ( kind = 8 ) TABLE(M,N), the table data.
    !
    implicit none

    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    integer ( kind = 4 ) ierror
    character ( len = * ) input_filename
    integer ( kind = 4 ) input_status
    integer ( kind = 4 ) input_unit
    integer ( kind = 4 ) j
    character ( len = 255 ) line
    real ( kind = 8 ) table(m,n)
    real ( kind = 8 ) x(m)

    ierror = 0

    call get_unit ( input_unit )

    open ( unit = input_unit, file = input_filename, status = 'old', &
        iostat = input_status )

    if ( input_status /= 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
        write ( *, '(a,i8)' ) '  Could not open the input file "' // &
            trim ( input_filename ) // '" on unit ', input_unit
        stop
    end if

    j = 0

    do while ( j < n )

        read ( input_unit, '(a)', iostat = input_status ) line

        if ( input_status /= 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8MAT_DATA_READ - Fatal error!'
            write ( *, '(a)' ) '  Error while reading lines of data.'
            write ( *, '(a,i8)' ) '  Number of values expected per line M = ', m
            write ( *, '(a,i8)' ) '  Number of data lines read, J =         ', j
            write ( *, '(a,i8)' ) '  Number of data lines needed, N =       ', n
            stop
        end if

        if ( line(1:1) == '#' .or. len_trim ( line ) == 0 ) then
            cycle
        end if

        call s_to_r8vec ( line, m, x, ierror )

        if ( ierror /= 0 ) then
            cycle
        end if

        j = j + 1

        table(1:m,j) = x(1:m)

    end do

    close ( unit = input_unit )

    return
    end
    subroutine r8mat_header_read ( input_filename, m, n )

    !*****************************************************************************80
    !
    !! R8MAT_HEADER_READ reads the header from an R8MAT file.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
    !
    !    Output, integer ( kind = 4 ) M, spatial dimension.
    !
    !    Output, integer ( kind = 4 ) N, the number of points.
    !
    implicit none

    character ( len = * ) input_filename
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    call file_column_count ( input_filename, m )

    if ( m <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
        write ( *, '(a)' ) '  to count the number of data columns in'
        write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
        stop
    end if

    call file_row_count ( input_filename, n )

    if ( n <= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  There was some kind of I/O problem while trying'
        write ( *, '(a)' ) '  to count the number of data rows in'
        write ( *, '(a)' ) '  the file "' // trim ( input_filename ) // '".'
        stop
    end if

    return
    end
    subroutine r8mat_transpose_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

    !*****************************************************************************80
    !
    !! R8MAT_TRANSPOSE_PRINT_SOME prints some of an R8MAT, transposed.
    !
    !  Discussion:
    !
    !    An R8MAT is an array of R8 values.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 June 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) M, N, the number of rows and columns.
    !
    !    Input, real ( kind = 8 ) A(M,N), an M by N matrix to be printed.
    !
    !    Input, integer ( kind = 4 ) ILO, JLO, the first row and column to print.
    !
    !    Input, integer ( kind = 4 ) IHI, JHI, the last row and column to print.
    !
    !    Input, character ( len = * ) TITLE, an optional title.
    !
    implicit none

    integer ( kind = 4 ), parameter :: incx = 5
    integer ( kind = 4 ) m
    integer ( kind = 4 ) n

    real ( kind = 8 ) a(m,n)
    character ( len = 14 ) ctemp(incx)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) i2
    integer ( kind = 4 ) i2hi
    integer ( kind = 4 ) i2lo
    integer ( kind = 4 ) ihi
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) inc
    integer ( kind = 4 ) j
    integer ( kind = 4 ) j2hi
    integer ( kind = 4 ) j2lo
    integer ( kind = 4 ) jhi
    integer ( kind = 4 ) jlo
    character ( len = * ) title

    if ( 0 < len_trim ( title ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) trim ( title )
    end if

    do i2lo = max ( ilo, 1 ), min ( ihi, m ), incx

        i2hi = i2lo + incx - 1
        i2hi = min ( i2hi, m )
        i2hi = min ( i2hi, ihi )

        inc = i2hi + 1 - i2lo

        write ( *, '(a)' ) ' '

        do i = i2lo, i2hi
            i2 = i + 1 - i2lo
            write ( ctemp(i2), '(i8,6x)' ) i
        end do

        write ( *, '(''  Row   '',5a14)' ) ctemp(1:inc)
        write ( *, '(a)' ) '  Col'
        write ( *, '(a)' ) ' '

        j2lo = max ( jlo, 1 )
        j2hi = min ( jhi, n )

        do j = j2lo, j2hi

            do i2 = 1, inc
                i = i2lo - 1 + i2
                write ( ctemp(i2), '(g14.6)' ) a(i,j)
            end do

            write ( *, '(i5,1x,5a14)' ) j, ( ctemp(i), i = 1, inc )

        end do

    end do

    return
    end
    subroutine s_blank_delete ( s )

    !*****************************************************************************80
    !
    !! S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
    !
    !  Discussion:
    !
    !    All TAB characters are also removed.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    26 July 1998
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !
    implicit none

    character c
    integer ( kind = 4 ) get
    integer ( kind = 4 ) put
    integer ( kind = 4 ) nchar
    character ( len = * ) s
    character, parameter :: TAB = char ( 9 )

    put = 0
    nchar = len_trim ( s )

    do get = 1, nchar

        c = s(get:get)

        if ( c /= ' ' .and. c /= TAB ) then
            put = put + 1
            s(put:put) = c
        end if

    end do

    s(put+1:nchar) = ' '

    return
    end
    subroutine s_to_r8 ( s, dval, ierror, length )

    !*****************************************************************************80
    !
    !! S_TO_R8 reads an R8 from a string.
    !
    !  Discussion:
    !
    !    The routine will read as many characters as possible until it reaches
    !    the end of the string, or encounters a character which cannot be
    !    part of the number.
    !
    !    Legal input is:
    !
    !       1 blanks,
    !       2 '+' or '-' sign,
    !       2.5 blanks
    !       3 integer part,
    !       4 decimal point,
    !       5 fraction part,
    !       6 'E' or 'e' or 'D' or 'd', exponent marker,
    !       7 exponent sign,
    !       8 exponent integer part,
    !       9 exponent decimal point,
    !      10 exponent fraction part,
    !      11 blanks,
    !      12 final comma or semicolon,
    !
    !    with most quantities optional.
    !
    !  Example:
    !
    !    S                 DVAL
    !
    !    '1'               1.0
    !    '     1   '       1.0
    !    '1A'              1.0
    !    '12,34,56'        12.0
    !    '  34 7'          34.0
    !    '-1E2ABCD'        -100.0
    !    '-1X2ABCD'        -1.0
    !    ' 2E-1'           0.2
    !    '23.45'           23.45
    !    '-4.2E+2'         -420.0
    !    '17d2'            1700.0
    !    '-14e-2'         -0.14
    !    'e2'              100.0
    !    '-12.73e-9.23'   -12.73 * 10.0^(-9.23)
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, the string containing the
    !    data to be read.  Reading will begin at position 1 and
    !    terminate at the end of the string, or when no more
    !    characters can be read to form a legal real.  Blanks,
    !    commas, or other nonnumeric data will, in particular,
    !    cause the conversion to halt.
    !
    !    Output, real ( kind = 8 ) DVAL, the value read from the string.
    !
    !    Output, integer ( kind = 4 ) IERROR, error flag.
    !    0, no errors occurred.
    !    1, 2, 6 or 7, the input number was garbled.  The
    !    value of IERROR is the last type of input successfully
    !    read.  For instance, 1 means initial blanks, 2 means
    !    a plus or minus sign, and so on.
    !
    !    Output, integer ( kind = 4 ) LENGTH, the number of characters read
    !    to form the number, including any terminating
    !    characters such as a trailing comma or blanks.
    !
    implicit none

    character c
    logical ch_eqi
    real ( kind = 8 ) dval
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) ihave
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) iterm
    integer ( kind = 4 ) jbot
    integer ( kind = 4 ) jsgn
    integer ( kind = 4 ) jtop
    integer ( kind = 4 ) length
    integer ( kind = 4 ) nchar
    integer ( kind = 4 ) ndig
    real ( kind = 8 ) rbot
    real ( kind = 8 ) rexp
    real ( kind = 8 ) rtop
    character ( len = * ) s

    nchar = len_trim ( s )

    ierror = 0
    dval = 0.0D+00
    length = -1
    isgn = 1
    rtop = 0
    rbot = 1
    jsgn = 1
    jtop = 0
    jbot = 1
    ihave = 1
    iterm = 0

    do

        length = length + 1

        if ( nchar < length+1 ) then
            exit
        end if

        c = s(length+1:length+1)
        !
        !  Blank character.
        !
        if ( c == ' ' ) then

            if ( ihave == 2 ) then

            else if ( ihave == 6 .or. ihave == 7 ) then
                iterm = 1
            else if ( 1 < ihave ) then
                ihave = 11
            end if
            !
            !  Comma.
            !
        else if ( c == ',' .or. c == ';' ) then

            if ( ihave /= 1 ) then
                iterm = 1
                ihave = 12
                length = length + 1
            end if
            !
            !  Minus sign.
            !
        else if ( c == '-' ) then

            if ( ihave == 1 ) then
                ihave = 2
                isgn = -1
            else if ( ihave == 6 ) then
                ihave = 7
                jsgn = -1
            else
                iterm = 1
            end if
            !
            !  Plus sign.
            !
        else if ( c == '+' ) then

            if ( ihave == 1 ) then
                ihave = 2
            else if ( ihave == 6 ) then
                ihave = 7
            else
                iterm = 1
            end if
            !
            !  Decimal point.
            !
        else if ( c == '.' ) then

            if ( ihave < 4 ) then
                ihave = 4
            else if ( 6 <= ihave .and. ihave <= 8 ) then
                ihave = 9
            else
                iterm = 1
            end if
            !
            !  Scientific notation exponent marker.
            !
        else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

            if ( ihave < 6 ) then
                ihave = 6
            else
                iterm = 1
            end if
            !
            !  Digit.
            !
        else if (  ihave < 11 .and. lle ( '0', c ) .and. lle ( c, '9' ) ) then

            if ( ihave <= 2 ) then
                ihave = 3
            else if ( ihave == 4 ) then
                ihave = 5
            else if ( ihave == 6 .or. ihave == 7 ) then
                ihave = 8
            else if ( ihave == 9 ) then
                ihave = 10
            end if

            call ch_to_digit ( c, ndig )

            if ( ihave == 3 ) then
                rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
            else if ( ihave == 5 ) then
                rtop = 10.0D+00 * rtop + real ( ndig, kind = 8 )
                rbot = 10.0D+00 * rbot
            else if ( ihave == 8 ) then
                jtop = 10 * jtop + ndig
            else if ( ihave == 10 ) then
                jtop = 10 * jtop + ndig
                jbot = 10 * jbot
            end if
            !
            !  Anything else is regarded as a terminator.
            !
        else
            iterm = 1
        end if
        !
        !  If we haven't seen a terminator, and we haven't examined the
        !  entire string, go get the next character.
        !
        if ( iterm == 1 ) then
            exit
        end if

    end do
    !
    !  If we haven't seen a terminator, and we have examined the
    !  entire string, then we're done, and LENGTH is equal to NCHAR.
    !
    if ( iterm /= 1 .and. length+1 == nchar ) then
        length = nchar
    end if
    !
    !  Number seems to have terminated.  Have we got a legal number?
    !  Not if we terminated in states 1, 2, 6 or 7!
    !
    if ( ihave == 1 .or. ihave == 2 .or. ihave == 6 .or. ihave == 7 ) then
        ierror = ihave
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Illegal or nonnumeric input:'
        write ( *, '(a)' ) '    ' // trim ( s )
        return
    end if
    !
    !  Number seems OK.  Form it.
    !
    if ( jtop == 0 ) then
        rexp = 1.0D+00
    else
        if ( jbot == 1 ) then
            rexp = 10.0D+00 ** ( jsgn * jtop )
        else
            rexp = 10.0D+00 ** ( real ( jsgn * jtop, kind = 8 ) &
                / real ( jbot, kind = 8 ) )
        end if
    end if

    dval = real ( isgn, kind = 8 ) * rexp * rtop / rbot

    return
    end
    subroutine s_to_r8vec ( s, n, rvec, ierror )

    !*****************************************************************************80
    !
    !! S_TO_R8VEC reads an R8VEC from a string.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    07 September 2004
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, the string to be read.
    !
    !    Input, integer ( kind = 4 ) N, the number of values expected.
    !
    !    Output, real ( kind = 8 ) RVEC(N), the values read from the string.
    !
    !    Output, integer ( kind = 4 ) IERROR, error flag.
    !    0, no errors occurred.
    !    -K, could not read data for entries -K through N.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierror
    integer ( kind = 4 ) ilo
    integer ( kind = 4 ) lchar
    real ( kind = 8 ) rvec(n)
    character ( len = * ) s

    i = 0
    ierror = 0
    ilo = 1

    do while ( i < n )

        i = i + 1

        call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

        if ( ierror /= 0 ) then
            ierror = -i
            exit
        end if

        ilo = ilo + lchar

    end do

    return
    end
    subroutine s_word_count ( s, nword )

    !*****************************************************************************80
    !
    !! S_WORD_COUNT counts the number of "words" in a string.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 April 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, character ( len = * ) S, the string to be examined.
    !
    !    Output, integer ( kind = 4 ) NWORD, the number of "words" in the string.
    !    Words are presumed to be separated by one or more blanks.
    !
    implicit none

    logical blank
    integer ( kind = 4 ) i
    integer ( kind = 4 ) lens
    integer ( kind = 4 ) nword
    character ( len = * ) s

    nword = 0
    lens = len ( s )

    if ( lens <= 0 ) then
        return
    end if

    blank = .true.

    do i = 1, lens

        if ( s(i:i) == ' ' ) then
            blank = .true.
        else if ( blank ) then
            nword = nword + 1
            blank = .false.
        end if

    end do

    return
    end
    subroutine swapec ( i, top, btri, bedg, point_num, point_xy, tri_num, &
        tri_vert, tri_nabe, stack, ierr )

    !*****************************************************************************80
    !
    !! SWAPEC swaps diagonal edges until all triangles are Delaunay.
    !
    !  Discussion:
    !
    !    The routine swaps diagonal edges in a 2D triangulation, based on
    !    the empty circumcircle criterion, until all triangles are Delaunay,
    !    given that I is the index of the new vertex added to the triangulation.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    14 July 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) I, the index of the new vertex.
    !
    !    Input/output, integer ( kind = 4 ) TOP, the index of the top of the stack.
    !    On output, TOP is zero.
    !
    !    Input/output, integer ( kind = 4 ) BTRI, BEDG; on input, if positive, are
    !    the triangle and edge indices of a boundary edge whose updated indices
    !    must be recorded.  On output, these may be updated because of swaps.
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
    !
    !    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
    !    of the points.
    !
    !    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
    !
    !    Input/output, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle
    !    incidence list.  May be updated on output because of swaps.
    !
    !    Input/output, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle
    !    neighbor list; negative values are used for links of the counter-clockwise
    !    linked list of boundary edges;  May be updated on output because of swaps.
    !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    !
    !    Workspace, integer ( kind = 4 ) STACK(MAXST); on input, entries 1 through
    !    TOP contain the indices of initial triangles (involving vertex I)
    !    put in stack; the edges opposite I should be in interior;  entries
    !    TOP+1 through MAXST are used as a stack.
    !
    !    Output, integer ( kind = 4 ) IERR is set to 8 for abnormal return.
    !
    implicit none

    integer ( kind = 4 ) point_num
    integer ( kind = 4 ) tri_num

    integer ( kind = 4 ) a
    integer ( kind = 4 ) b
    integer ( kind = 4 ) bedg
    integer ( kind = 4 ) btri
    integer ( kind = 4 ) c
    integer ( kind = 4 ) diaedg
    integer ( kind = 4 ) e
    integer ( kind = 4 ) ee
    integer ( kind = 4 ) em1
    integer ( kind = 4 ) ep1
    integer ( kind = 4 ) f
    integer ( kind = 4 ) fm1
    integer ( kind = 4 ) fp1
    integer ( kind = 4 ) i
    integer ( kind = 4 ) ierr
    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) l
    integer ( kind = 4 ) r
    integer ( kind = 4 ) s
    integer ( kind = 4 ) stack(point_num)
    integer ( kind = 4 ) swap
    integer ( kind = 4 ) t
    integer ( kind = 4 ) top
    integer ( kind = 4 ) tri_nabe(3,tri_num)
    integer ( kind = 4 ) tri_vert(3,tri_num)
    integer ( kind = 4 ) tt
    integer ( kind = 4 ) u
    real ( kind = 8 ) point_xy(2,point_num)
    real ( kind = 8 ) x
    real ( kind = 8 ) y
    !
    !  Determine whether triangles in stack are Delaunay, and swap
    !  diagonal edge of convex quadrilateral if not.
    !
    x = point_xy(1,i)
    y = point_xy(2,i)

    do

        if ( top <= 0 ) then
            exit
        end if

        t = stack(top)
        top = top - 1

        if ( tri_vert(1,t) == i ) then
            e = 2
            b = tri_vert(3,t)
        else if ( tri_vert(2,t) == i ) then
            e = 3
            b = tri_vert(1,t)
        else
            e = 1
            b = tri_vert(2,t)
        end if

        a = tri_vert(e,t)
        u = tri_nabe(e,t)

        if ( tri_nabe(1,u) == t ) then
            f = 1
            c = tri_vert(3,u)
        else if ( tri_nabe(2,u) == t ) then
            f = 2
            c = tri_vert(1,u)
        else
            f = 3
            c = tri_vert(2,u)
        end if

        swap = diaedg ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,c), &
            point_xy(2,c), point_xy(1,b), point_xy(2,b) )

        if ( swap == 1 ) then

            em1 = i4_wrap ( e - 1, 1, 3 )
            ep1 = i4_wrap ( e + 1, 1, 3 )
            fm1 = i4_wrap ( f - 1, 1, 3 )
            fp1 = i4_wrap ( f + 1, 1, 3 )

            tri_vert(ep1,t) = c
            tri_vert(fp1,u) = i
            r = tri_nabe(ep1,t)
            s = tri_nabe(fp1,u)
            tri_nabe(ep1,t) = u
            tri_nabe(fp1,u) = t
            tri_nabe(e,t) = s
            tri_nabe(f,u) = r

            if ( 0 < tri_nabe(fm1,u) ) then
                top = top + 1
                stack(top) = u
            end if

            if ( 0 < s ) then

                if ( tri_nabe(1,s) == u ) then
                    tri_nabe(1,s) = t
                else if ( tri_nabe(2,s) == u ) then
                    tri_nabe(2,s) = t
                else
                    tri_nabe(3,s) = t
                end if

                top = top + 1

                if ( point_num < top ) then
                    ierr = 8
                    return
                end if

                stack(top) = t

            else

                if ( u == btri .and. fp1 == bedg ) then
                    btri = t
                    bedg = e
                end if

                l = - ( 3 * t + e - 1 )
                tt = t
                ee = em1

                do while ( 0 < tri_nabe(ee,tt) )

                    tt = tri_nabe(ee,tt)

                    if ( tri_vert(1,tt) == a ) then
                        ee = 3
                    else if ( tri_vert(2,tt) == a ) then
                        ee = 1
                    else
                        ee = 2
                    end if

                end do

                tri_nabe(ee,tt) = l

            end if

            if ( 0 < r ) then

                if ( tri_nabe(1,r) == t ) then
                    tri_nabe(1,r) = u
                else if ( tri_nabe(2,r) == t ) then
                    tri_nabe(2,r) = u
                else
                    tri_nabe(3,r) = u
                end if

            else

                if ( t == btri .and. ep1 == bedg ) then
                    btri = u
                    bedg = f
                end if

                l = - ( 3 * u + f - 1 )
                tt = u
                ee = fm1

                do while ( 0 < tri_nabe(ee,tt) )

                    tt = tri_nabe(ee,tt)

                    if ( tri_vert(1,tt) == b ) then
                        ee = 3
                    else if ( tri_vert(2,tt) == b ) then
                        ee = 1
                    else
                        ee = 2
                    end if

                end do

                tri_nabe(ee,tt) = l

            end if

        end if

    end do

    return
    end
    subroutine timestamp ( )

    !*****************************************************************************80
    !
    !! TIMESTAMP prints the current YMDHMS date as a time stamp.
    !
    !  Example:
    !
    !    31 May 2001   9:45:54.872 AM
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    18 May 2013
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    None
    !
    implicit none

    character ( len = 8 ) ampm
    integer ( kind = 4 ) d
    integer ( kind = 4 ) h
    integer ( kind = 4 ) m
    integer ( kind = 4 ) mm
    character ( len = 9 ), parameter, dimension(12) :: month = (/ &
        'January  ', 'February ', 'March    ', 'April    ', &
        'May      ', 'June     ', 'July     ', 'August   ', &
        'September', 'October  ', 'November ', 'December ' /)
    integer ( kind = 4 ) n
    integer ( kind = 4 ) s
    integer ( kind = 4 ) values(8)
    integer ( kind = 4 ) y

    call date_and_time ( values = values )

    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    mm = values(8)

    if ( h < 12 ) then
        ampm = 'AM'
    else if ( h == 12 ) then
        if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
        else
            ampm = 'PM'
        end if
    else
        h = h - 12
        if ( h < 12 ) then
            ampm = 'PM'
        else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
                ampm = 'Midnight'
            else
                ampm = 'AM'
            end if
        end if
    end if

    write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
        d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

    return
    end
    subroutine vbedg ( x, y, point_num, point_xy, tri_num, tri_vert, tri_nabe, &
        ltri, ledg, rtri, redg )

    !*****************************************************************************80
    !
    !! VBEDG determines which boundary edges are visible to a point.
    !
    !  Discussion:
    !
    !    The point (X,Y) is assumed to be outside the convex hull of the
    !    region covered by the 2D triangulation.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 August 2001
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Barry Joe.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Barry Joe,
    !    GEOMPACK - a software package for the generation of meshes
    !    using geometric algorithms,
    !    Advances in Engineering Software,
    !    Volume 13, pages 325-331, 1991.
    !
    !  Parameters:
    !
    !    Input, real ( kind = 8 ) X, Y, the coordinates of a point outside
    !    the convex hull of the current triangulation.
    !
    !    Input, integer ( kind = 4 ) POINT_NUM, the number of points.
    !
    !    Input, real ( kind = 8 ) POINT_XY(2,POINT_NUM), the coordinates
    !    of the vertices.
    !
    !    Input, integer ( kind = 4 ) TRI_NUM, the number of triangles.
    !
    !    Input, integer ( kind = 4 ) TRI_VERT(3,TRI_NUM), the triangle incidence
    !    list.
    !
    !    Input, integer ( kind = 4 ) TRI_NABE(3,TRI_NUM), the triangle neighbor
    !    list; negative values are used for links of a counter clockwise linked
    !    list of boundary edges;
    !      LINK = -(3*I + J-1) where I, J = triangle, edge index.
    !
    !    Input/output, integer ( kind = 4 ) LTRI, LEDG.  If LTRI /= 0 then these
    !    values are assumed to be already computed and are not changed, else they
    !    are updated.  On output, LTRI is the index of boundary triangle to the
    !    left of the leftmost boundary triangle visible from (X,Y), and LEDG is
    !    the boundary edge of triangle LTRI to the left of the leftmost boundary
    !    edge visible from (X,Y).  1 <= LEDG <= 3.
    !
    !    Input/output, integer ( kind = 4 ) RTRI.  On input, the index of the
    !    boundary triangle to begin the search at.  On output, the index of the
    !    rightmost boundary triangle visible from (X,Y).
    !
    !    Input/output, integer ( kind = 4 ) REDG, the edge of triangle RTRI that
    !    is visible from (X,Y).  1 <= REDG <= 3.
    !
    implicit none

    integer ( kind = 4 ) point_num
    integer ( kind = 4 ) tri_num

    integer ( kind = 4 ) a
    integer ( kind = 4 ) b
    integer ( kind = 4 ) e
    integer ( kind = 4 ) i4_wrap
    integer ( kind = 4 ) l
    logical ldone
    integer ( kind = 4 ) ledg
    integer ( kind = 4 ) lr
    integer ( kind = 4 ) lrline
    integer ( kind = 4 ) ltri
    real ( kind = 8 ) point_xy(2,point_num)
    integer ( kind = 4 ) redg
    integer ( kind = 4 ) rtri
    integer ( kind = 4 ) t
    integer ( kind = 4 ) tri_nabe(3,tri_num)
    integer ( kind = 4 ) tri_vert(3,tri_num)
    real ( kind = 8 ) x
    real ( kind = 8 ) y
    !
    !  Find the rightmost visible boundary edge using links, then possibly
    !  leftmost visible boundary edge using triangle neighbor information.
    !
    if ( ltri == 0 ) then
        ldone = .false.
        ltri = rtri
        ledg = redg
    else
        ldone = .true.
    end if

    do

        l = -tri_nabe(redg,rtri)
        t = l / 3
        e = mod ( l, 3 ) + 1
        a = tri_vert(e,t)

        if ( e <= 2 ) then
            b = tri_vert(e+1,t)
        else
            b = tri_vert(1,t)
        end if

        lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
            point_xy(2,b), 0.0D+00 )

        if ( lr <= 0 ) then
            exit
        end if

        rtri = t
        redg = e

    end do

    if ( ldone ) then
        return
    end if

    t = ltri
    e = ledg

    do

        b = tri_vert(e,t)
        e = i4_wrap ( e-1, 1, 3 )

        do while ( 0 < tri_nabe(e,t) )

            t = tri_nabe(e,t)

            if ( tri_vert(1,t) == b ) then
                e = 3
            else if ( tri_vert(2,t) == b ) then
                e = 1
            else
                e = 2
            end if

        end do

        a = tri_vert(e,t)

        lr = lrline ( x, y, point_xy(1,a), point_xy(2,a), point_xy(1,b), &
            point_xy(2,b), 0.0D+00 )

        if ( lr <= 0 ) then
            exit
        end if

    end do

    ltri = t
    ledg = e

    return
    end

    subroutine i4vec2_sort_a ( n, a1, a2 )

    !*****************************************************************************80
    !
    !! I4VEC2_SORT_A ascending sorts a vector of pairs of integers.
    !
    !  Discussion:
    !
    !    Each item to be sorted is a pair of integers (I,J), with the I
    !    and J values stored in separate vectors A1 and A2.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    25 September 2001
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of items of data.
    !
    !    Input/output, integer ( kind = 4 ) A1(N), A2(N), the data to be sorted.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) a1(n)
    integer ( kind = 4 ) a2(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ) temp

    if ( n <= 1 ) then
        return
    end if
    !
    !  Initialize.
    !
    i = 0
    indx = 0
    isgn = 0
    j = 0
    !
    !  Call the external heap sorter.
    !
    do

        call sort_heap_external ( n, indx, i, j, isgn )
        !
        !  Interchange the I and J objects.
        !
        if ( 0 < indx ) then

            temp  = a1(i)
            a1(i) = a1(j)
            a1(j) = temp

            temp  = a2(i)
            a2(i) = a2(j)
            a2(j) = temp
            !
            !  Compare the I and J objects.
            !
        else if ( indx < 0 ) then

            call i4vec2_compare ( n, a1, a2, i, j, isgn )

        else if ( indx == 0 ) then

            exit

        end if

    end do

    return
    end
    subroutine i4vec2_sorted_unique ( n, a1, a2, unique_num )

    !*****************************************************************************80
    !
    !! I4VEC2_SORTED_UNIQUE gets the unique elements in a sorted I4VEC2.
    !
    !  Discussion:
    !
    !    Item I is stored as the pair A1(I), A2(I).
    !
    !    The items must have been sorted, or at least it must be the
    !    case that equal items are stored in adjacent vector locations.
    !
    !    If the items were not sorted, then this routine will only
    !    replace a string of equal values by a single representative.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    09 July 2000
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of items.
    !
    !    Input/output, integer ( kind = 4 ) A1(N), A2(N).
    !    On input, the array of N items.
    !    On output, an array of unique items.
    !
    !    Output, integer ( kind = 4 ) UNIQUE_NUM, the number of unique items.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) a1(n)
    integer ( kind = 4 ) a2(n)
    integer ( kind = 4 ) itest
    integer ( kind = 4 ) unique_num

    unique_num = 0

    if ( n <= 0 ) then
        return
    end if

    unique_num = 1

    do itest = 2, n

        if ( a1(itest) /= a1(unique_num) .or. a2(itest) /= a2(unique_num) ) then

            unique_num = unique_num + 1

            a1(unique_num) = a1(itest)
            a2(unique_num) = a2(itest)

        end if

    end do

    return
    end

    subroutine sort_heap_external ( n, indx, i, j, isgn )

    !*****************************************************************************80
    !
    !! SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
    !
    !  Discussion:
    !
    !    The actual list of data is not passed to the routine.  Hence this
    !    routine may be used to sort integers, real ( kind = 8 )s, numbers, names,
    !    dates, shoe sizes, and so on.  After each call, the routine asks
    !    the user to compare or interchange two items, until a special
    !    return value signals that the sorting is completed.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    05 February 2004
    !
    !  Author:
    !
    !    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
    !    FORTRAN90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Albert Nijenhuis, Herbert Wilf,
    !    Combinatorial Algorithms,
    !    Academic Press, 1978, second edition,
    !    ISBN 0-12-519260-6.
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of items to be sorted.
    !
    !    Input/output, integer ( kind = 4 ) INDX, the main communication signal.
    !
    !    The user must set INDX to 0 before the first call.
    !    Thereafter, the user should not change the value of INDX until
    !    the sorting is done.
    !
    !    On return, if INDX is
    !
    !      greater than 0,
    !      * interchange items I and J;
    !      * call again.
    !
    !      less than 0,
    !      * compare items I and J;
    !      * set ISGN = -1 if I < J, ISGN = +1 if J < I;
    !      * call again.
    !
    !      equal to 0, the sorting is done.
    !
    !    Output, integer ( kind = 4 ) I, J, the indices of two items.
    !    On return with INDX positive, elements I and J should be interchanged.
    !    On return with INDX negative, elements I and J should be compared, and
    !    the result reported in ISGN on the next call.
    !
    !    Input, integer ( kind = 4 ) ISGN, results of comparison of elements I
    !    and J.  (Used only when the previous call returned INDX less than 0).
    !    ISGN <= 0 means I is less than or equal to J;
    !    0 <= ISGN means I is greater than or equal to J.
    !
    implicit none

    integer ( kind = 4 ) i
    integer ( kind = 4 ), save :: i_save = 0
    integer ( kind = 4 ) indx
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j
    integer ( kind = 4 ), save :: j_save = 0
    integer ( kind = 4 ), save :: k = 0
    integer ( kind = 4 ), save :: k1 = 0
    integer ( kind = 4 ) n
    integer ( kind = 4 ), save :: n1 = 0
    !
    !  INDX = 0: This is the first call.
    !
    if ( indx == 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
        !
        !  INDX < 0: The user is returning the results of a comparison.
        !
    else if ( indx < 0 ) then

        if ( indx == -2 ) then

            if ( isgn < 0 ) then
                i_save = i_save + 1
            end if

            j_save = k1
            k1 = i_save
            indx = -1
            i = i_save
            j = j_save
            return

        end if

        if ( 0 < isgn ) then
            indx = 2
            i = i_save
            j = j_save
            return
        end if

        if ( k <= 1 ) then

            if ( n1 == 1 ) then
                i_save = 0
                j_save = 0
                indx = 0
            else
                i_save = n1
                n1 = n1 - 1
                j_save = 1
                indx = 1
            end if

            i = i_save
            j = j_save
            return

        end if

        k = k - 1
        k1 = k
        !
        !  0 < INDX, the user was asked to make an interchange.
        !
    else if ( indx == 1 ) then

        k1 = k

    end if

    do

        i_save = 2 * k1

        if ( i_save == n1 ) then
            j_save = k1
            k1 = i_save
            indx = -1
            i = i_save
            j = j_save
            return
        else if ( i_save <= n1 ) then
            j_save = i_save + 1
            indx = -2
            i = i_save
            j = j_save
            return
        end if

        if ( k <= 1 ) then
            exit
        end if

        k = k - 1
        k1 = k

    end do

    if ( n1 == 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
    else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
    end if

    return
    end


    subroutine i4vec2_compare ( n, a1, a2, i, j, isgn )

    !*****************************************************************************80
    !
    !! I4VEC2_COMPARE compares pairs of integers stored in two vectors.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    22 October 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Input, integer ( kind = 4 ) N, the number of data items.
    !
    !    Input, integer ( kind = 4 ) A1(N), A2(N), contain the two components
    !    of each item.
    !
    !    Input, integer ( kind = 4 ) I, J, the items to be compared.
    !
    !    Output, integer ( kind = 4 ) ISGN, the results of the comparison:
    !    -1, item I < item J,
    !     0, item I = item J,
    !    +1, item J < item I.
    !
    implicit none

    integer ( kind = 4 ) n

    integer ( kind = 4 ) a1(n)
    integer ( kind = 4 ) a2(n)
    integer ( kind = 4 ) i
    integer ( kind = 4 ) isgn
    integer ( kind = 4 ) j

    isgn = 0

    if ( a1(i) < a1(j) ) then

        isgn = -1

    else if ( a1(i) == a1(j) ) then

        if ( a2(i) < a2(j) ) then
            isgn = -1
        else if ( a2(i) < a2(j) ) then
            isgn = 0
        else if ( a2(j) < a2(i) ) then
            isgn = +1
        end if

    else if ( a1(j) < a1(i) ) then

        isgn = +1

    end if

    return
    end
