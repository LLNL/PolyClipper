########################################
PolyClipper concepts and data types
########################################

PolyClipper only has a few data types ``Vector2d``, ``Vector3d``, ``Plane2d``, ``Plane3d``, ``Vertex2d``, and ``Vertex3d``, all defined inside the C++ namespace ``PolyClipper``.  There are also ``Polygon`` and ``Polyhedron`` types defined, but these are simply aliases for ``std::vector<PolyClipper::Vertex2d>`` and ``std::vector<PolyClipper::Vertex3d>``.

Collections of vertices encapsulate all the necessary information to specify Polygons and Polyhedra.  This works because vertices contain their connectivity to neighbor vertices in the object.  For example, consider a hexagon as an example of a ``Polygon``:

.. image:: Vertex2d.*
   :width: 200
   :alt: You should see a labeled polygon here

The ``Vertex2d`` PolyClipper structure for vertex 2 includes the links to vertices (1, 3) as neighbors in counterclockwise order.

----------
Vector2d
----------

Vector2d represents a 2D vector coordinate in space, :math:`(x,y)`.  It supports a small number of simple vector manipulation operations:

..
   .. module:: PolyClipper
   .. autoclass:: Vector2d
      :members:
