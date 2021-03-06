########################################
PolyClipper data types
########################################

As mentioned in section :ref:`PolyClipper concepts`, PolyClipper includes the data types ``Vector2d``, ``Vector3d``, ``Plane2d``, ``Plane3d``, ``Vertex2d``, and ``Vertex3d``, all defined inside the C++ namespace ``PolyClipper``.  For 2D (polygon) operations include the header ``polyclipper2d.hh``, while for 3D (polyhedron) include ``polyclipper3d.hh``.
In this section we describe the C++ interface for each of these types.

.. note::
   Python interface notes

   The Python interface is identical to C++, with template parameters defaulted to PolyClippper's Vector2d/Vector3d implementations.  We also define a ``Polygon`` and ``Polyhedron`` class for in Python, which are ``std::vector<PolyClipper::Vertex2d<>>`` and ``std::vector<PolyClipper::Vertex3d<>>`` under the hood.

..
   Vector2d
   ----------

.. cpp:namespace:: PolyClipper

Vector classes
--------------------

Vector2d
~~~~~~~~

.. cpp:class:: Vector2d

  Vector2d represents a 2D vector coordinate in space, :math:`(x,y)`.  It supports a small number of simple vector manipulation operations:

  .. cpp:member:: double Vector2d::x

     :math:`x` coordinate

  .. cpp:member:: double Vector2d::y

     :math:`y` coordinate

  .. cpp:function:: Vector2d::Vector2d()

     Construct a ``Vector2d(0,0)``

  .. cpp:function:: Vector2d::Vector2d(double X, double Y)

     Construct a ``Vector2d(X,Y)``

  .. cpp:function:: double Vector2d::dot(const Vector2d& b) const

     Returns the dot product: :math:`\vec{a} \cdot \vec{b}`

  .. cpp:function:: double Vector2d::cross(const Vector2d& b) const

     Returns the :math:`z` component of the cross product: :math:`\vec{a} \times \vec{b}`

  .. cpp:function:: double Vector2d::magnitude2() const

     Returns the square of the magnitude of the vector: :math:`\vec{a} \cdot \vec{a}`

  .. cpp:function:: double Vector2d::magnitude() const

     Returns the magnitude of the vector: :math:`\sqrt{\vec{a} \cdot \vec{a}}`

  .. cpp:function:: Vector2d& Vector2d::operator*=(const double b)

     Inplace multiplication by a scalar: :math:`\vec{a} = b \vec{a}`

  .. cpp:function:: Vector2d& Vector2d::operator/=(const double b)

     Inplace division by a scalar: :math:`\vec{a} = \vec{a}/b`

  .. cpp:function:: Vector2d& Vector2d::operator+=(const Vector2d& b)

     Inplace addition with another vector: :math:`\vec{a} = \vec{a} + \vec{b}`

  .. cpp:function:: Vector2d& Vector2d::operator-=(const Vector2d& b)

     Inplace subtraction with another vector: :math:`\vec{a} = \vec{a} - \vec{b}`

  .. cpp:function:: Vector2d Vector2d::operator*(const double b) const

     Returns the result of multiplication by a scalar: :math:`b \vec{a}`

  .. cpp:function:: Vector2d Vector2d::operator/(const double b) const

     Returns the result of division by a scalar: :math:`\vec{a}/b`

  .. cpp:function:: Vector2d Vector2d::operator+(const Vector2d& b) const

     Returns the result of addition with another vector: :math:`\vec{a} + \vec{b}`

  .. cpp:function:: Vector2d Vector2d::operator-(const Vector2d& b) const

     Returns the result of subtraction with another vector: :math:`\vec{a} - \vec{b}`

  .. cpp:function:: Vector2d Vector2d::operator-() const

     Returns the negative of this vector: :math:`-\vec{a}`

  .. cpp:function:: Vector2d Vector2d::unitVector() const

     Returns the unit vector pointing in the direction of this one: :math:`\vec{a}/\sqrt{\vec{a} \cdot \vec{a}}`.

     If :math:`\vec{a} = (0,0)`, returns the unit vector in the :math:`x` direction: :math:`(1,0)`.
     
Vector3d
~~~~~~~~

.. cpp:class:: Vector3d

  Vector3d represents a 3D vector coordinate in space, :math:`(x,y,z)`.  It supports a small number of simple vector manipulation operations:

  .. cpp:member:: double Vector3d::x

     :math:`x` coordinate

  .. cpp:member:: double Vector3d::y

     :math:`y` coordinate

  .. cpp:member:: double Vector3d::z

     :math:`z` coordinate

  .. cpp:function:: Vector3d::Vector3d()

     Construct a ``Vector3d(0,0,0)``

  .. cpp:function:: Vector3d::Vector3d(double X, double Y, double Z)

     Construct a ``Vector3d(X,Y,Z)``

  .. cpp:function:: double Vector3d::dot(const Vector3d& b) const

     Returns the dot product: :math:`\vec{a} \cdot \vec{b}`

  .. cpp:function:: Vector3d Vector3d::cross(const Vector3d& b) const

     Returns the cross product: :math:`\vec{a} \times \vec{b}`

  .. cpp:function:: double Vector3d::magnitude2() const

     Returns the square of the magnitude of the vector: :math:`\vec{a} \cdot \vec{a}`

  .. cpp:function:: double Vector3d::magnitude() const

     Returns the magnitude of the vector: :math:`\sqrt{\vec{a} \cdot \vec{a}}`

  .. cpp:function:: Vector3d& Vector3d::operator*=(const double b)

     Inplace multiplication by a scalar: :math:`\vec{a} = b \vec{a}`

  .. cpp:function:: Vector3d& Vector3d::operator/=(const double b)

     Inplace division by a scalar: :math:`\vec{a} = \vec{a}/b`

  .. cpp:function:: Vector3d& Vector3d::operator+=(const Vector3d& b)

     Inplace addition with another vector: :math:`\vec{a} = \vec{a} + \vec{b}`

  .. cpp:function:: Vector3d& Vector3d::operator-=(const Vector3d& b)

     Inplace subtraction with another vector: :math:`\vec{a} = \vec{a} - \vec{b}`

  .. cpp:function:: Vector3d Vector3d::operator*(const double b) const

     Returns the result of multiplication by a scalar: :math:`b \vec{a}`

  .. cpp:function:: Vector3d Vector3d::operator/(const double b) const

     Returns the result of division by a scalar: :math:`\vec{a}/b`

  .. cpp:function:: Vector3d Vector3d::operator+(const Vector3d& b) const

     Returns the result of addition with another vector: :math:`\vec{a} + \vec{b}`

  .. cpp:function:: Vector3d Vector3d::operator-(const Vector3d& b) const

     Returns the result of subtraction with another vector: :math:`\vec{a} - \vec{b}`

  .. cpp:function:: Vector3d Vector3d::operator-() const

     Returns the negative of this vector: :math:`-\vec{a}`

  .. cpp:function:: Vector3d Vector3d::unitVector() const

     Returns the unit vector pointing in the direction of this one: :math:`\vec{a}/\sqrt{\vec{a} \cdot \vec{a}}`.

     If :math:`\vec{a} = (0,0,0)`, returns the unit vector in the :math:`x` direction: :math:`(1,0,0)`.
     
Plane classes
--------------------

.. cpp:class:: template<typename VA> Plane

    Plane represents a plane for clipping polytopes, and is templated on ``VA`` (a Vector adapter) type describing how to use the appropriate geometric Vector type.

    .. note::
      In the headers ``polyclipper2d.hh`` and ``polyclipper3d.hh`` we define the typedefs::

        using Plane2d = Plane<internal::VectorAdapter<Vector2d>>;
        using Plane3d = Plane<internal::VectorAdapter<Vector3d>>;

      (both in the namespace ``PolyClipper``) for convenience when working with PolyClipper's native Vectors.

    A plane is stored as a unit normal and closest signed distance from the plane to the origin: :math:`(\hat{n}, d)`.  The signed distance from the plane to any point :math:`\vec{p}` is

    .. math::
       d_s(\vec{p}) = (\vec{p} - \vec{p}_0) \cdot \hat{n} = d + \vec{p} \cdot \hat{n},

    where :math:`\vec{p}_0` is any point in the plane.  Note with this definition the :math:`d` parameter defining the plane is :math:`d = -\vec{p}_0 \cdot \hat{n}`.

  .. cpp:type:: Vector VA::VECTOR

  .. cpp:member:: double Plane::dist

     The minimum signed distance from the origin to the plane :math:`d`.

  .. cpp:member:: Vector Plane::normal

     The unit normal to the plane :math:`\hat{n}`.

  .. cpp:member:: int Plane::ID

     An optional integer identification number for the plane.  This is used by ``Vertex<VA>`` to record which plane(s) are responsible for creating the vertex.

  .. cpp:function:: Plane::Plane()

     Default constructor -- implies {:math:`\hat{n}, d`, ID} = {(1,0,0), 0.0, std::numeric_limits<int>::min()}

  .. cpp:function:: Plane::Plane(const double d, const Vector& nhat)

     Construct with {:math:`\hat{n}, d`, ID} = {nhat, d, std::numeric_limits<int>::min()}

  .. cpp:function:: Plane::Plane(const Vector& p, const Vector& nhat)

     Construct specifying the normal and a point in the plane, so {:math:`\hat{n}, d`, ID} = {nhat, :math:`-p\cdot\hat{n}`, std::numeric_limits<int>::min()}

  .. cpp:function:: Plane::Plane(const Vector& p, const Vector& nhat, const int id)

     Construct specifying the normal, a point in the plane, and ID, so {:math:`\hat{n}, d`, ID} = {nhat, :math:`-p\cdot\hat{n}`, id}

Vertex classes
--------------------

Vertex2d
~~~~~~~~

.. cpp:class:: template<typename VA = internal::VectorAdapter<Vector2d>> Vertex2d

  Vertex2d is used to encode Polygons in 2d.  A vertex includes a position and the connectivity to neighboring vertices in the Polygon.  In this 2d case, the connectivity is always 2 vertices, ordered such that going from the first neighbor, to this vertex, and on to the last neighbor goes around the Polygon in the counter-clockwise direction.  This is illustrated in the Polygon examples in :ref:`PolyClipper concepts`.

  .. note::
     Note that while ``Vertex2d`` is templated on a geometric Vector trait class, the template argument defaults to an implementation appropriate for use with ``PolyClipper::Vector2d``.

  .. cpp:type:: Vector VA::VECTOR

  .. cpp:member:: Vector Vertex2d::position

     The position of the vertex in :math:`(x,y)` coordinates.

  .. cpp:member:: std::pair<int, int> Vertex2d::neighbors

     The neighbor vertices this vertex is connected too.  These should be listed in counter-clockwise order going around the Polygon, so that ``neighbors.first`` is clockwise and ``neighbors.second`` is counter-clockwise from this vertex.

  .. cpp:member:: int Vertex2d::comp

     An internal state integer, for comparing this vertex to planes.  Used and overwritten during clipping operations.

  .. cpp:member:: int Vertex2d::ID

     An optional ID index for this vertex.  Used and overwritten during clipping operations.

  .. cpp:member:: std::set<int> Vertex2d::clips

     The set of Plane2d ID's responsible for creating this vertex during clipping operations.  Used and overwritten during clipping operations.

  .. cpp:function:: Vertex2d::Vertex2d()

     Default constructor, sets member data to {position, neighbors, comp, ID, clips} = {(0,0), (), 1, -1, {}}

  .. cpp:function:: Vertex2d::Vertex2d(const Vector& pos)

     Construct with just the position

  .. cpp:function:: Vertex2d::Vertex2d(const Vector& pos, const int c)

     Construct with {position, comp} = {pos, c}

Vertex3d
~~~~~~~~

.. cpp:class:: template<typename VA = internal::VectorAdapter<Vector3d>> Vertex3d

  Vertex3d is used to encode Polyhedra in 3d.  A vertex includes a position and the connectivity to neighboring vertices in the Polyhedron.  For Polyhedra, the neighbor connectivity should be 3 or more neighbors, listed counter-clockwise as viewed from the exterior side of the vertex (see the illustrations in :ref:`PolyClipper concepts` for examples).

  .. note::
     Note that while ``Vertex3d`` is templated on a geometric Vector trait class, the template argument defaults to an implementation appropriate for use with ``PolyClipper::Vector3d``.

  .. cpp:type:: Vector VA::VECTOR

  .. cpp:member:: Vector Vertex3d::position

     The position of the vertex in :math:`(x,y,z)` coordinates.

  .. cpp:member:: std::vector<int> Vertex3d::neighbors

     The neighbor vertices this vertex is connected too, listed in counter-clockwise order as viewed from the exterior of the Polyhedron.

  .. cpp:member:: int Vertex3d::comp

     An internal state integer, for comparing this vertex to planes.  Used and overwritten during clipping operations.

  .. cpp:member:: int Vertex3d::ID

     An optional ID index for this vertex.  Used and overwritten during clipping operations.

  .. cpp:member:: std::set<int> Vertex3d::clips

     The set of Plane3d ID's responsible for creating this vertex during clipping operations.  Used and overwritten during clipping operations.

  .. cpp:function:: Vertex3d::Vertex3d()

     Default constructor, sets member data to {position, neighbors, comp, ID, clips} = {(0,0,0), (), 1, -1, {}}

  .. cpp:function:: Vertex3d::Vertex3d(const Vector& pos)

     Construct with just the position

  .. cpp:function:: Vertex3d::Vertex3d(const Vector& pos, const int c)

     Construct with {position, comp} = {pos, c}
