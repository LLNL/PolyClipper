########################################
Polygon methods
########################################

These are the methods provided in the PolyClipper namespace (C++) and/or PolyClipper module (Python) for manipulating Polygons.  Recall that a Polygon is simply a ``std::vector<PolyClipper::Vertex2d>``, as described in :ref:`PolyClipper concepts`.  All methods here are included in the single header file ``polyclipper.hh``.

.. cpp:namespace:: PolyClipper

.. cpp:function:: void initializePolygon(Polygon& poly, const std::vector<PolyClipper::Vector2d>& positions, const std::vector<std::vector<int>>& neighbors)

   Initialize a PolyClipper::Polygon by constructing the necessary PolyClipper::Vertex2d objects described by the ``positions`` and ``neighbors``.  Note that the length of these arrays should be identical (i.e., the number of vertices in the resulting Polygon).  Each element of the ``neighbors`` array should be 2 elements long, listing the (clockwise, counterclockwise) neighbors for the vertex at the corresponding position in the ``positions`` array.
