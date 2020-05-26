########################################
PolyClipper concepts
########################################

PolyClipper only has a few data types ``Vector2d``, ``Vector3d``, ``Plane2d``, ``Plane3d``, ``Vertex2d``, and ``Vertex3d``, all defined inside the C++ namespace ``PolyClipper``.  There are also ``Polygon`` and ``Polyhedron`` types defined, but these are simply aliases for ``std::vector<PolyClipper::Vertex2d>`` and ``std::vector<PolyClipper::Vertex3d>``.

Collections of vertices encapsulate all the necessary information to specify Polygons and Polyhedra.  This works because vertices contain their connectivity to neighbor vertices in the object.  For example, consider a hexagon as an example of a ``Polygon``:

.. image:: Vertex2d.*
   :width: 200
   :alt: You should see a labeled polygon here

The ``Vertex2d`` PolyClipper structure for vertex 2 includes the links to vertices (1, 3) as neighbors in counterclockwise order.

Similarly a cube can be represented as a collection of ``Vertex3d``:

.. image:: Cube.*
   :width: 200
   :alt: You should see a labeled cube here

In this case the neighbors of vertex 6 will be (5, 2, 7), listed counterclockwise as viewed from the exterior of the polyhedron.  Note vertices of Polyhedra can have 3 or more neighbors -- as a simple example consider a pyramid:

.. image:: Pyramid.*
   :width: 200
   :alt: You should see a labeld pyramid here

Here we can see that vertex 4 will have 4 neighbors: (0, 1, 2, 3).

Clipping operations
========================================

PolyClippers main functionality is the clipping of Polygons/Polyhedra by arbitray planes.  Planes are defined by a unit normal and distance from the plane to the origin, :math:`(\hat{n}, d)`.  Note given some point :math:`\vec{p_0}` in the plane, :math:`d = -\hat{n}\cdot\vec{p}_0`.

Points are defined as above a plane when they are on the side the plane normal :math:`\hat{n}` points toward, and below if on the opposite side.  This can be defined mathematically in terms of the signed distance of a point :math:`\vec{p}` from the plane:

.. math::
   d_s(\vec{p}) = (\vec{p} - \vec{p}_0) \cdot \hat{n} = d + \vec{p} \cdot \hat{n}.

Now if :math:`d_s(\vec{p}) > 0` then :math:`\vec{p}` is above the plane, while :math:`d_s(\vec{p}) < 0` implies it is below.

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

  Double clip:  [{pos=(3.000000 2.000000), neighbors=(3 4), ID=0, clips=( 20 )}, 
                 {pos=(1.000000 2.000000), neighbors=(5 2), ID=1, clips=( )}, 
                 {pos=(0.000000 2.000000), neighbors=(1 6), ID=2, clips=( )}, 
                 {pos=(3.500000 2.000000), neighbors=(7 0), ID=3, clips=( 10 20 )}, 
                 {pos=(2.083333 1.083333), neighbors=(0 7), ID=4, clips=( 30 )}, 
                 {pos=(1.875000 1.125000), neighbors=(6 1), ID=5, clips=( 30 )}, 
                 {pos=(0.000000 1.500000), neighbors=(2 5), ID=6, clips=( 30 )}, 
                 {pos=(2.954545 0.909091), neighbors=(4 3), ID=7, clips=( 10 30 )}]

Note in this case we have created two independent loop of vertices in our resulting Polygon.

