########################################
Functions for manipulating polyhedra
########################################

These are the methods provided in the ``PolyClipper`` namespace (C++) and/or ``PolyClipper`` module (Python) for manipulating polyhedra.  Recall that a polyhedron is simply a ``std::vector<PolyClipper::Vertex3d>``, as described in :ref:`PolyClipper concepts`.  All methods here are included in the single header file ``polyclipper3d.hh``.

.. cpp:namespace:: PolyClipper

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector2d>> \
                  std::vector<std::vector<int>> splitIntoTriangles(const std::vector<Vertex2d<VA>>& poly, \
                                                                   const double tol = 0.0)

   Initialize a polyhedron by constructing the necessary ``PolyClipper::Vertex3d`` objects described by the ``positions`` and ``neighbors``.  Note that the length of these arrays should be identical (i.e., the number of vertices in the resulting polyhedron).  Each element of the ``neighbors`` array is an array of the neighbors for the vertex at the corresponding index in the ``positions`` array, where these neighbors are listed in counter-clockwise order around the vertex as viewed from the exterior of the polyhedron.
   
   See the examples in :ref:`PolyClipper concepts`.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  std::string polyhedron2string(const std::vector<Vertex3d<VA>>& poly)

   Creates a human-readable string representation of the polyhedron.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  void moments(double& zerothMoment, typename VA::VECTOR& firstMoment, \
                               const std::vector<Vertex3d<VA>>& polyhedron)

   Compute the zeroth (volume) and first (centroid) moment of a polyhedron.

   .. note::
      In Python this method has a different signature, as the moments are returned as a Python tuple:

      .. py:function:: moments(polyhedron) -> (double, Vector3d)

   .. warning::
      While the volume returned in this function is always correct, the centroid is only correct for convex polyhedra.  This should be generalized to work for all polyhedra in a future release.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  void clipPolyhedron(std::vector<Vertex3d<VA>>& poly, \
                                      const std::vector<Plane<VA>>& planes)

   Clip a polyhedron by a set of planes in place.  Examples are shown in :ref:`Clipping operations`.  The region of the polyhedron above the each plane (in the direction of the plane normal) is retained.

   After clipping, the ``Vertex3d::ID`` and ``Vertex3d::clips`` attribute of the vertices in the polyhedron are modified, such that ID holds a unique identifier for each remaining vertex, and clips holds the ID's of any planes used to create the vertex.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  void collapseDegenerates(std::vector<Vertex3d<VA>>& poly, \
                                           const double tol)

   Remove redundant vertices in the polyhedron in place, such that any edge of the input polyhedron with length less than ``tol`` is removed and their vertices combined.  It is possible entire faces of the polyhedron may be removed in this process, though only edge lengths are examined.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  std::vector<std::vector<int>> extractFaces(const std::vector<Vertex3d<VA>>& poly)

   Return an array of the vertex indices that represent the faces of the polyhedron.  The length of the returned array is the number of faces, while each element is the set of vertices (3 or more) in counter-clockwise order (viewed from the exterior of the polyhedron).

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  std::vector<std::set<int>> commonFaceClips(const std::vector<Vertex3d<VA>>& poly, \
                                                             const std::vector<std::vector<int>>& faces)

   Compute the unique plane ID's responsible for creating each face in the polyhedron.  Assumes the ``Vertex3d::clips`` attribute has been filled in by clipping the polyhedron.

.. cpp:function:: template<typename VA = internal::VectorAdapter<Vector3d>> \
                  std::vector<std::vector<int>> splitIntoTetrahedra(const std::vector<Vertex3d<VA>>& poly, \
                                                                    const double tol = 0.0)

   Return a tetrahedralization of the polyhedron.  The result is an array of quartets, with each quartet the indices of the vertices making up each tetrahedron.  The ``tol`` attribute is used to reject any tetrahedron with volumes less than ``tol``.

   .. warning::
      This method currently only works for convex polyhedra, and raises an assertion if called with a non-convex polyhedron.
