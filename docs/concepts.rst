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

PolyClippers main functionality is the clipping of Polygons/Polyhedra by arbitray planes.  Planes are defined by a unit normal and distance from the plane to the origin, :math:`(\hat{n}, d)`.  Note given some point :math:`\vec{p}` in the plane, :math:`d = -\hat{n}\cdot\vec{p}`.
