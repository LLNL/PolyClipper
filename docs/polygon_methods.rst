########################################
Polygon methods
########################################

These are the methods provided in the PolyClipper namespace (C++) and/or PolyClipper module (Python) for manipulating Polygons.  Recall that a Polygon is simply a ``std::vector<PolyClipper::Vertex2d>``, as described in :ref:`PolyClipper concepts`.  All methods here are included in the single header file ``polyclipper.hh``.

.. cpp:namespace:: PolyClipper

.. cpp:function:: void initializePolygon(Polygon& poly, const std::vector<PolyClipper::Vector2d>& positions, const std::vector<std::vector<int>>& neighbors)

   Initialize a PolyClipper::Polygon by constructing the necessary PolyClipper::Vertex2d objects described by the ``positions`` and ``neighbors``.  Note that the length of these arrays should be identical (i.e., the number of vertices in the resulting Polygon).  Each element of the ``neighbors`` array should be 2 elements long, listing the (clockwise, counterclockwise) neighbors for the vertex at the corresponding index in the ``positions`` array.

   See the examples in :ref:`PolyClipper concepts`.

.. cpp:function:: std::string polygon2string(const Polygon& poly)

   Creates a human-readable string representation of the Polygon.

.. cpp:function:: void moments(double& zerothMoment, PolyClipper::Vector2d& firstMoment, const Polygon& polygon)

   Compute the zeroth (area) and first (centroid) moment of a Polygon.

   .. Note::

      In Python this method has a different signature, as the moments are returned as a Python tuple:

      .. py:function:: moments(polygon) -> (double, Vector2d)

   .. Note::

      While the area returned in this function is always correct, the centroid is only correct for convex polygons.  This is planned to be generalized to work for all polygons in a future release.

.. cpp:function:: void clipPolygon(Polygon& poly, const std::vector<Plane2d>& planes)

   Clip a Polygon by a set of planes in place.  Examples are shown in :ref:`Clipping operations`.  The region of the Polygon above the each plane (in the direction of the plane normal) is retained.

   After clipping, the ``Vertex2d::ID`` and ``Vertex2d::clips`` attribute of the vertices in the Polygon are modified, such that ID holds a unique identifier for each remaining vertex, and clips holds the ID's of any planes used to create the vertex.

.. cpp:function:: void collapseDegenerates(Polygon& poly, const double tol)

   Remove redundant vertices in the Polygon in place, such that any edge of the input Polygon with length less than ``tol`` is removed and their vertices combined.

.. cpp:function:: std::vector<std::vector<int>> extractFaces(const Polygon& poly)

   Return an array of the Vertex indices that represent the faces (or edges) of the Polygon.  The length of the returned array is the number of faces, and each element is of length 2 representing the face/edge.

.. cpp:function:: std::vector<std::set<int>> commonFaceClips(const Polygon& poly, const std::vector<std::vector<int>>& faces)

   Compute the unique Plane ID's responsible for creating each face in the Polygon.  Assumes the ``Vertex2d::clips`` attribute has been filled in by clipping the Polygon.

.. cpp:function:: std::vector<std::vector<int>> splitIntoTriangles(const Polygon& poly, const double tol = 0.0)

   Return a triangulation of the Polygon.  The result is an array of triples, with each triple the indices of the vertices making up each triangle.  The ``tol`` attribute is used to reject any triangles with areas less than ``tol``.

   .. Note::

      This method currently only works for convex Polygons, and raises an assertion if called with a non-convex Polygon.
