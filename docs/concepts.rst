########################################
PolyClipper concepts and data types
########################################

PolyClipper only has a few data types ``Vector2d``, ``Vector3d``, ``Plane2d``, ``Plane3d``, ``Vertex2d``, and ``Vertex3d``, all defined inside the C++ namespace ``PolyClipper``.  There are also ``Polygon`` and ``Polyhedron`` types defined, but these are simply aliases for ``std::vector<PolyClipper::Vertex2d>`` and ``std::vector<PolyClipper::Vertex3d>``.

Collections of vertices (``Vertex2d`` and ``Vertex3d``) encapsulate all the necessary information to specify Polygons and Polyhedra.  A Polygon

.. image:: Vertex2d.*
   :width: 200
   :alt: You should see a labelled polygon here

----------
Vector2d
----------

Vector2d represents a 2D vector coordinate in space, :math:`(x,y)`.  It supports a small number of simple vector manipulation operations:

..
   .. module:: PolyClipper
   .. autoclass:: Vector2d
      :members:
