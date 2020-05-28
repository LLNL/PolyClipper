########################################
Polyhedron methods
########################################

These are the methods provided in the PolyClipper namespace (C++) and/or PolyClipper module (Python) for manipulating Polyhedra.  Recall that a Polyhedron is simply a ``std::vector<PolyClipper::Vertex3d>``, as described in :ref:`PolyClipper concepts`.  All methods here are included in the single header file ``polyclipper.hh``.

.. cpp:namespace:: PolyClipper

.. cpp:function:: void initializePolyhedron(Polyhedron& poly, const std::vector<PolyClipper::Vector3d>& positions, const std::vector<std::vector<int>>& neighbors)

   Initialize a PolyClipper::Polyhedron by constructing the necessary PolyClipper::Vertex3d objects described by the ``positions`` and ``neighbors``.  Note that the length of these arrays should be identical (i.e., the number of vertices in the resulting Polyhedron).  Each element of the ``neighbors`` array is an array of the neighbors for the vertex at the corresponding index in the ``positions`` array, where these neighbors are listed in counter-clockwise order around the vertex as viewed from the exterior of the Polyhedron.
   
   See the examples in :ref:`PolyClipper concepts`.

.. cpp:function:: std::string polyhedron2string(const Polyhedron& poly)

   Creates a human-readable string representation of the Polyhedron.

.. cpp:function:: void moments(double& zerothMoment, PolyClipper::Vector3d& firstMoment, const Polyhedron& polyhedron)

   Compute the zeroth (volume) and first (centroid) moment of a Polyhedron.

   .. Note::

      In Python this method has a different signature, as the moments are returned as a Python tuple:

      .. py:function:: moments(polyhedron) -> (double, Vector3d)

   .. Note::

      While the volume returned in this function is always correct, the centroid is only correct for convex polyhedra.  This is planned to be generalized to work for all polyhedra in a future release.

.. cpp:function:: void clipPolyhedron(Polyhedron& poly, const std::vector<Plane3d>& planes)

   Clip a Polyhedron by a set of planes in place.  Examples are shown in :ref:`Clipping operations`.  The region of the Polyhedron above the each plane (in the direction of the plane normal) is retained.

   After clipping, the ``Vertex3d::ID`` and ``Vertex3d::clips`` attribute of the vertices in the Polyhedron are modified, such that ID holds a unique identifier for each remaining vertex, and clips holds the ID's of any planes used to create the vertex.

.. cpp:function:: void collapseDegenerates(Polyhedron& poly, const double tol)

   Remove redundant vertices in the Polyhedron in place, such that any edge of the input Polyhedron with length less than ``tol`` is removed and their vertices combined.  It is possible entire faces of the Polyhedron may be removed in this process, though only edge lengths are examined.

.. cpp:function:: std::vector<std::vector<int>> extractFaces(const Polyhedron& poly)

   Return an array of the Vertex indices that represent the faces of the Polyhedron.  The length of the returned array is the number of faces, while each element is the set of vertices (3 or more) in counter-clockwise order (viewed from the exterior of the Polyhedron).

.. cpp:function:: std::vector<std::set<int>> commonFaceClips(const Polyhedron& poly, const std::vector<std::vector<int>>& faces)

   Compute the unique Plane ID's responsible for creating each face in the Polyhedron.  Assumes the ``Vertex3d::clips`` attribute has been filled in by clipping the Polyhedron.

.. cpp:function:: std::vector<std::vector<int>> splitIntoTetrahedra(const Polyhedron& poly, const double tol = 0.0)

   Return a tetrahedralization of the Polyhedron.  The result is an array of quartets, with each quartet the indices of the vertices making up each tetrahedron.  The ``tol`` attribute is used to reject any tetrahedra with volumes less than ``tol``.

   .. Note::

      This method currently only works for convex Polyhedra, and raises an assertion if called with a non-convex Polyhedron.
