########################################
PolyClipper concepts
########################################

PolyClipper only has a few data types (``Vector2d``, ``Vector3d``, ``Plane2d``, ``Plane3d``, ``Vertex2d``, and ``Vertex3d``), all defined inside the C++ namespace ``PolyClipper``.  The main methods of PolyClipper are free functions for clipping and manipulating polytopes (i.e., polygons and polyhedra), which are represented as collections of ``Vertex2d`` (for polygons) and ``Vertex3d`` (for polyhedra).  Note that while PolyClipper provides native (``Vector2d``, ``Vector3d``) types, PolyClipper's functions are templated on the Vector type (actually a trait class describing the Vector operations), allowing the C++ user to use their own native geometric Vectors if so desired.  The built in Python module provided by PolyClipper is instantiated using PolyClipper's native Vector's, so if you want to provide Python versions with your own Vector types you must bind those yourself.  Further discussion and examples of this generalized C++ functionality can be found in the :ref:`Using PolyClipper with a user defined Vector type`.

Collections of vertices encapsulate all the necessary information to specify Polygons and Polyhedra.  This works because vertices contain their connectivity to neighbor vertices in the object.  For example, consider a hexagon as an example of a ``Polygon``:

.. image:: Vertex2d.*
   :width: 200
   :alt: You should see a labeled polygon here

The ``Vertex2d`` PolyClipper structure for vertex 2 includes the links to vertices (1, 3) as neighbors in counterclockwise order.

Similarly a cube can be represented as a collection of ``Vertex3d``:

.. image:: Cube.*
   :width: 200
   :alt: You should see a labeled cube here

In this case the neighbors of vertex 6 will be (5, 2, 7), listed counterclockwise as viewed from the exterior of the polyhedron.  Note while vertices of Polygons will always have 2 neighbors, vertices of Polyhedra can have 3 or more neighbors.  A pyramid is a simple example of this:

.. image:: Pyramid.*
   :width: 200
   :alt: You should see a labeld pyramid here

Here we can see that vertex 4 will have 4 neighbors: (0, 1, 2, 3).

Clipping operations
========================================

PolyClippers main functionality is the clipping of Polygons/Polyhedra by arbitray planes.  Planes are defined by a unit normal and distance from the plane to the origin, :math:`(\hat{n}, d)`.  Note given some point :math:`\vec{p_0}` in the plane, :math:`d = -\hat{n}\cdot\vec{p}_0`.

Points are defined as above a plane when they are on the side the plane normal :math:`\hat{n}` points toward, and below if on the opposite side.  This can be defined mathematically in terms of the signed distance of a point :math:`\vec{p}` from the plane:

.. math::
   d_s(\vec{p}) = (\vec{p} - \vec{p}_0) \cdot \hat{n} = d + \vec{p} \cdot \hat{n},

where again :math:`\vec{p}_0` is some point in the plane.  Now if :math:`d_s(\vec{p}) > 0` then :math:`\vec{p}` is above the plane, while :math:`d_s(\vec{p}) < 0` implies it is below.

Clipping a Polygon/Polyhedron in PolyClipper is defined such that the portion of the volume above the plane is retained, while that below is clipped and removed.  Consider for instance a non-convex polygon:

.. image:: notched_polygon.*
           :width: 400
           :alt: You should see a non-convex polygon

In this example the coordinates are shown in magenta on the interior of the polygon, and the vertex numbers in black on the exterior.  In PolyClippers Python interface this Polygon can be constructed using::

  notched_points = [Vector2d(*coords)
                    for coords in [(0,0), (4,0), (4,2), (3,2), (2,1), (1,2), (0,2)]]
  n = len(notched_points)
  notched_neighbors = [[(i - 1) % n, (i + 1) % n] for i in xrange(n)]
  poly = Polygon()
  initializePolygon(poly, notched_points, notched_neighbors)
  print "Starting poly: ", list(poly)

where we have used the fact the vertices are numbered sequentially in counter-clockwise order to generate the neighbor lists for each vertex.  Since Polygons are simply lists of ``Vertex2d``, the final print outputs something like::

  Starting poly:  [{pos=(0.000000 0.000000), neighbors=(6 1), ID=-1, clips=( )}, 
                   {pos=(4.000000 0.000000), neighbors=(0 2), ID=-1, clips=( )}, 
                   {pos=(4.000000 2.000000), neighbors=(1 3), ID=-1, clips=( )}, 
                   {pos=(3.000000 2.000000), neighbors=(2 4), ID=-1, clips=( )}, 
                   {pos=(2.000000 1.000000), neighbors=(3 5), ID=-1, clips=( )}, 
                   {pos=(1.000000 2.000000), neighbors=(4 6), ID=-1, clips=( )}, 
                   {pos=(0.000000 2.000000), neighbors=(5 0), ID=-1, clips=( )}]


Clipping by a single plane
----------------------------------------

We can clip this Polygon by a single plane defined by a {point, normal} of :math:`\{(3, 1), \widehat{{(-1, 0.5)}}\}` (where the wide-hat symbol implies constructing the unit vector) with::

  planes = [Plane2d(Vector2d(3, 1), Vector2d(-1, 0.5).unitVector(), 10)]
  clipPolygon(poly, planes)
  print "Single clip: ", list(poly)

resulting in

.. image:: notched_polygon_clip1.*
           :width: 400
           :alt: You should see a clipped polygon

and vertices now printed as::

  Single clip:  [{pos=(0.000000 0.000000), neighbors=(4 5), ID=0, clips=( )}, 
                 {pos=(3.000000 2.000000), neighbors=(6 2), ID=1, clips=( )}, 
                 {pos=(2.000000 1.000000), neighbors=(1 3), ID=2, clips=( )}, 
                 {pos=(1.000000 2.000000), neighbors=(2 4), ID=3, clips=( )}, 
                 {pos=(0.000000 2.000000), neighbors=(3 0), ID=4, clips=( )}, 
                 {pos=(2.500000 0.000000), neighbors=(0 6), ID=5, clips=( 10 )},
                 {pos=(3.500000 2.000000), neighbors=(5 1), ID=6, clips=( 10 )}]

Note the two new vertices (5 & 6) created by the clip plane have the ID of the plane that created them (10) listed in their ``clips`` parameter.  It is not required to construct Planes with unique ID's like this, in which case all Planes have the same ID and this ``clips`` parameter is less useful.

If we instead clip the original Polygon by the plane with the unit normal flipped in the opposite direction we get the other part of the Polygon::

  planes = [Plane2d(Vector2d(3, 1), Vector2d(1, -0.5).unitVector(), 20)]
  clipPolygon(poly, planes)
  print "Reverse clip: ", list(poly)

.. image:: notched_polygon_clip2.*
           :width: 400
           :alt: You should see a clipped polygon

and vertices::

  Reverse clip:  [{pos=(4.000000 0.000000), neighbors=(2 1), ID=0, clips=( )}, 
                  {pos=(4.000000 2.000000), neighbors=(0 3), ID=1, clips=( )}, 
                  {pos=(2.500000 0.000000), neighbors=(3 0), ID=2, clips=( 20 )}, 
                  {pos=(3.500000 2.000000), neighbors=(1 2), ID=3, clips=( 20 )}]

Clipping by multiple planes
----------------------------------------

Similarly we can clip by multiple planes simultaneously::

  planes = [Plane2d(Vector2d(3, 1), Vector2d(-1, 0.5).unitVector(), 10),
            Plane2d(Vector2d(2, 1.1), Vector2d(1, 5).unitVector(), 30)]
  clipPolygon(poly, planes)
  print "Double clip: ", list(poly)

.. image:: notched_polygon_clip3.*
           :width: 400
           :alt: You should see a clipped polygon

and the vertices are now::

  Double clip:  [{pos=(3.000000 2.000000), neighbors=(3 4), ID=0, clips=( )}, 
                 {pos=(1.000000 2.000000), neighbors=(5 2), ID=1, clips=( )}, 
                 {pos=(0.000000 2.000000), neighbors=(1 6), ID=2, clips=( )}, 
                 {pos=(3.500000 2.000000), neighbors=(7 0), ID=3, clips=( 10 )}, 
                 {pos=(2.083333 1.083333), neighbors=(0 7), ID=4, clips=( 30 )}, 
                 {pos=(1.875000 1.125000), neighbors=(6 1), ID=5, clips=( 30 )}, 
                 {pos=(0.000000 1.500000), neighbors=(2 5), ID=6, clips=( 30 )}, 
                 {pos=(2.954545 0.909091), neighbors=(4 3), ID=7, clips=( 10 30 )}]

Note in this case we have created two independent loop of vertices in our resulting Polygon.

3D Polyhedral clipping
-----------------------

The interface for clipping polyhedra in 3D is very similar to the 2D polygonal examples above.  For instance, if we extrude the non-convex polygonal example in in the :math:`z` direction for our initial polyhedron using the following Python code::

  notched_points = [Vector3d(*coords)
                    for coords in [(0,0,0), (4,0,0), (4,2,0), (3,2,0), (2,1,0), (1,2,0), (0,2,0),
                                   (0,0,1), (4,0,1), (4,2,1), (3,2,1), (2,1,1), (1,2,1), (0,2,1)]]
  notched_neighbors = [[7, 6, 1],   # 0
                       [0, 2, 8],   # 1
                       [1, 3, 9],   # 2
                       [4, 10, 2],  # 3
                       [5, 11, 3],  # 4
                       [6, 12, 4],  # 5
                       [13, 5, 0],  # 6
                       [8, 13, 0],  # 7
                       [1, 9, 7],   # 8
                       [2, 10, 8],  # 9
                       [9, 3, 11],  # 10
                       [10, 4, 12], # 11
                       [11, 5, 13], # 12
                       [7, 12, 6]]  # 13
  poly = Polyhedron()
  initializePolyhedron(poly, notched_points, notched_neighbors)

we get the following polyhedron:

.. image:: notched_polyhedron.*
           :width: 400
           :alt: You should see a clipped polyhedron

We can clip this example with a pair of planes::

  planes = [Plane3d(Vector3d(3, 1, 0), Vector3d(-1, 0.5, -1.5).unitVector(), 10),
            Plane3d(Vector3d(1, 1, 0), Vector3d(1, 0, -1).unitVector(), 30)]
  clipPolyhedron(poly, planes)
  print "Double clip: ", list(poly)

yielding a shape of:

.. image:: notched_polyhedron_clip3.*
           :width: 400
           :alt: You should see a clipped polyhedron

and vertex output::

  Double clip:  [{pos=(3.000000 2.000000 0.000000), neighbors=( 1 4 5 ), ID=0, clips=( )}, 
                 {pos=(2.000000 1.000000 0.000000), neighbors=( 2 6 0 ), ID=1, clips=( )}, 
                 {pos=(1.000000 2.000000 0.000000), neighbors=( 7 8 1 ), ID=2, clips=( )}, 
                 {pos=(2.500000 0.000000 0.000000), neighbors=( 5 9 10 ), ID=3, clips=( 10 )}, 
                 {pos=(3.000000 2.000000 0.333333), neighbors=( 6 5 0 ), ID=4, clips=( 10 )}, 
                 {pos=(3.500000 2.000000 0.000000), neighbors=( 4 3 0 ), ID=5, clips=( 10 )}, 
                 {pos=(2.000000 1.000000 0.666667), neighbors=( 11 4 1 ), ID=6, clips=( 10 )}, 
                 {pos=(1.000000 2.000000 0.000000), neighbors=( 10 8 2 ), ID=7, clips=( 30 )}, 
                 {pos=(1.000000 2.000000 0.000000), neighbors=( 7 11 2 ), ID=8, clips=( 30 )}, 
                 {pos=(1.600000 0.000000 0.600000), neighbors=( 11 10 3 ), ID=9, clips=( 10 30 )}, 
                 {pos=(1.000000 0.000000 0.000000), neighbors=( 9 7 3 ), ID=10, clips=( 30 )}, 
                 {pos=(1.833333 1.166667 0.833333), neighbors=( 8 9 6 ), ID=11, clips=( 10 30 )}]
