########################################
Functions for manipulating polygons
########################################

These are the methods provided in the ``PolyClipper`` namespace (C++) and/or ``PolyClipper`` module (Python) for manipulating polygons.  Recall that a polygon is simply a ``std::vector<PolyClipper::Vertex2d<>>``, as described in :ref:`PolyClipper concepts`.  All methods here are included in the single header file ``polyclipper2d.hh``.

.. cpp:namespace:: PolyClipper

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  void initializePolygon(std::vector<Vertex2d<VA>>& poly, const std::vector<typename VA::VECTOR>& positions, const std::vector<std::vector<int>>& neighbors)

   Initialize a polygon by constructing the necessary ``PolyClipper::Vertex2d`` objects described by the ``positions`` and ``neighbors``.  Note that the length of these arrays should be identical (i.e., the number of vertices in the resulting polygon).  Each element of the ``neighbors`` array should be 2 elements long, listing the (clockwise, counterclockwise) neighbors for the vertex at the corresponding index in the ``positions`` array.

   See the examples in :ref:`PolyClipper concepts`.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  std::string polygon2string(const std::vector<Vertex2d<VA>>& poly)

   Creates a human-readable string representation of the polygon.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  void moments(double& zerothMoment, typename VA::VECTOR& firstMoment, \
                               const std::vector<Vertex2d<VA>>& polygon)

   Compute the zeroth (area) and first (centroid) moment of a polygon.

   .. note::
      In Python this method has a different signature, as the moments are returned as a Python tuple:

      .. py:function:: moments(polygon) -> (double, Vector2d)

   .. warning::
      While the area returned in this function is always correct, the centroid is only correct for convex polygons.  This should be generalized to work for all polygons in a future release.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  void clipPolygon(std::vector<Vertex2d<VA>>& poly, \
                                   const std::vector<Plane<VA>>& planes)

   Clip a polygon by a set of planes in place.  Examples are shown in :ref:`Clipping operations`.  The region of the polygon above the each plane (in the direction of the plane normal) is retained.

   After clipping, the ``Vertex2d::ID`` and ``Vertex2d::clips`` attribute of the vertices in the polygon are modified, such that ID holds a unique identifier for each remaining vertex, and clips holds the ID's of any planes used to create the vertex.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  void collapseDegenerates(std::vector<Vertex2d<VA>>& poly, \
                                           const double tol)

   Remove redundant vertices in the polygon in place, such that any edge of the input polygon with length less than ``tol`` is removed and their vertices combined.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  std::vector<std::vector<int>> extractFaces(const std::vector<Vertex2d<VA>>& poly)

   Return an array of vertex indices that represent the faces (or edges) of the polygon.  The length of the returned array is the number of faces, and each element is of length 2 representing the face/edge.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  std::vector<std::set<int>> commonFaceClips(const std::vector<Vertex2d<VA>>& poly, \
                                                             const std::vector<std::vector<int>>& faces)

   Compute the unique plane ID's responsible for creating each face in the polygon.  Assumes the ``Vertex2d::clips`` attribute has been filled in by clipping the polygon.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  std::vector<std::vector<int>> splitIntoTriangles(const std::vector<Vertex2d<VA>>& poly, \
                                                                   const double tol = 0.0)

   Return a triangulation of the polygon.  The result is an array of triples, with each triple the indices of the vertices making up each triangle.  The ``tol`` attribute is used to reject any triangles with areas less than ``tol``.

   .. warning::
      This method currently only works for convex polygons, and raises an assertion if called with a non-convex polygon.
