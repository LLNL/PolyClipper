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
--------------------

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

where we have used the fact the vertices are numbered sequentially in counter-clockwise order to generate the neighbor lists for each vertex.  We can clip this Polygon with a single plane in Python 
