########################################
PolyClipper data types
########################################

As mentioned in section :ref:`PolyClipper concepts`, PolyClipper includes the data types ``Vector2d``, ``Vector3d``, ``Plane2d``, ``Plane3d``, ``Vertex2d``, and ``Vertex3d``, all defined inside the C++ namespace ``PolyClipper``.  The entirity of PolyClipper's interface (including these types) are in the single header file ``polyclipper.hh``.  In this section we describe the C++ interface for each of these types (the Python interface is identical).

..
   Note this file also defines the types ``PolyClipper::Polygon`` and ``PolyClipper::Polyhedron``, but these are simply aliases for ``std::vector<PolyClipper::Vertex2d>`` and ``std::vector<PolyClipper::Vertex3d>``.  

..
   Vector2d
   ----------

.. cpp:namespace:: PolyClipper

Vector classes
--------------------
..
  ------------------------------------------------------------------------------
  Vector2d
  ------------------------------------------------------------------------------

.. cpp:class:: Vector2d

  Vector2d represents a 2D vector coordinate in space, :math:`(x,y)`.  It supports a small number of simple vector manipulation operations:

  .. cpp:member:: double x

     :math:`x` coordinate

  .. cpp:member:: double y

     :math:`y` coordinate

  .. cpp:function:: Vector2d()

     Construct a ``Vector2d(0,0)``

  .. cpp:function:: Vector2d(double X, double Y)

     Construct a ``Vector2d(X,Y)``

  .. cpp:function:: double dot(const Vector2d& b) const

     Returns the dot product: :math:`\vec{a} \cdot \vec{b}`

  .. cpp:function:: double cross(const Vector2d& b) const

     Returns the :math:`z` component of the cross product: :math:`\vec{a} \times \vec{b}`

  .. cpp:function:: double magnitude2() const

     Returns the square of the magnitude of the vector: :math:`\vec{a} \cdot \vec{a}`

  .. cpp:function:: double magnitude() const

     Returns the magnitude of the vector: :math:`\sqrt{\vec{a} \cdot \vec{a}}`

  .. cpp:function:: Vector2d& operator*=(const double b)

     Inplace multiplication by a scalar: :math:`\vec{a} = b \vec{a}`

  .. cpp:function:: Vector2d& operator/=(const double b)

     Inplace division by a scalar: :math:`\vec{a} = \vec{a}/b`

  .. cpp:function:: Vector2d& operator+=(const Vector2d& b)

     Inplace addition with another vector: :math:`\vec{a} = \vec{a} + \vec{b}`

  .. cpp:function:: Vector2d& operator-=(const Vector2d& b)

     Inplace subtraction with another vector: :math:`\vec{a} = \vec{a} - \vec{b}`

  .. cpp:function:: Vector2d operator*(const double b) const

     Returns the result of multiplication by a scalar: :math:`b \vec{a}`

  .. cpp:function:: Vector2d operator/(const double b) const

     Returns the result of division by a scalar: :math:`\vec{a}/b`

  .. cpp:function:: Vector2d operator+(const Vector2d& b) const

     Returns the result of addition with another vector: :math:`\vec{a} + \vec{b}`

  .. cpp:function:: Vector2d operator-(const Vector2d& b) const

     Returns the result of subtraction with another vector: :math:`\vec{a} - \vec{b}`

  .. cpp:function:: Vector2d operator-() const

     Returns the negative of this vector: :math:`-\vec{a}`

  .. cpp:function:: Vector2d unitVector() const

     Returns the unit vector pointing in the direction of this one: :math:`\vec{a}/\sqrt{\vec{a} \cdot \vec{a}}`.

     If :math:`\vec{a} = (0,0)`, returns the unit vector in the :math:`x` direction: :math:`(1,0)`.
     
..
  ------------------------------------------------------------------------------
  Vector3d
  ------------------------------------------------------------------------------

.. cpp:class:: Vector3d

  Vector3d represents a 3D vector coordinate in space, :math:`(x,y,z)`.  It supports a small number of simple vector manipulation operations:

  .. cpp:member:: double x

     :math:`x` coordinate

  .. cpp:member:: double y

     :math:`y` coordinate

  .. cpp:member:: double z

     :math:`z` coordinate

  .. cpp:function:: Vector3d()

     Construct a ``Vector3d(0,0,0)``

  .. cpp:function:: Vector3d(double X, double Y, double Z)

     Construct a ``Vector3d(X,Y,Z)``

  .. cpp:function:: double dot(const Vector3d& b) const

     Returns the dot product: :math:`\vec{a} \cdot \vec{b}`

  .. cpp:function:: Vector3d cross(const Vector3d& b) const

     Returns the cross product: :math:`\vec{a} \times \vec{b}`

  .. cpp:function:: double magnitude2() const

     Returns the square of the magnitude of the vector: :math:`\vec{a} \cdot \vec{a}`

  .. cpp:function:: double magnitude() const

     Returns the magnitude of the vector: :math:`\sqrt{\vec{a} \cdot \vec{a}}`

  .. cpp:function:: Vector3d& operator*=(const double b)

     Inplace multiplication by a scalar: :math:`\vec{a} = b \vec{a}`

  .. cpp:function:: Vector3d& operator/=(const double b)

     Inplace division by a scalar: :math:`\vec{a} = \vec{a}/b`

  .. cpp:function:: Vector3d& operator+=(const Vector3d& b)

     Inplace addition with another vector: :math:`\vec{a} = \vec{a} + \vec{b}`

  .. cpp:function:: Vector3d& operator-=(const Vector3d& b)

     Inplace subtraction with another vector: :math:`\vec{a} = \vec{a} - \vec{b}`

  .. cpp:function:: Vector3d operator*(const double b) const

     Returns the result of multiplication by a scalar: :math:`b \vec{a}`

  .. cpp:function:: Vector3d operator/(const double b) const

     Returns the result of division by a scalar: :math:`\vec{a}/b`

  .. cpp:function:: Vector3d operator+(const Vector3d& b) const

     Returns the result of addition with another vector: :math:`\vec{a} + \vec{b}`

  .. cpp:function:: Vector3d operator-(const Vector3d& b) const

     Returns the result of subtraction with another vector: :math:`\vec{a} - \vec{b}`

  .. cpp:function:: Vector3d operator-() const

     Returns the negative of this vector: :math:`-\vec{a}`

  .. cpp:function:: Vector3d unitVector() const

     Returns the unit vector pointing in the direction of this one: :math:`\vec{a}/\sqrt{\vec{a} \cdot \vec{a}}`.

     If :math:`\vec{a} = (0,0,0)`, returns the unit vector in the :math:`x` direction: :math:`(1,0,0)`.
     
Plane classes
--------------------
..
  ------------------------------------------------------------------------------
  Plane2d
  ------------------------------------------------------------------------------

.. cpp:class:: Plane2d

   Plane2d represents a plane in the :math:`(x,y)` coordinate system for clipping Polygons.  A plane is stored as a unit normal and closest signed distance from the plane to the origin: :math:`(\hat{n}, d)`.  The signed distance from the plane to any point :math:`\vec{p}` is

  .. math::
     d_s(\vec{p}) = (\vec{p} - \vec{p}_0) \cdot \hat{n} = d + \vec{p} \cdot \hat{n},

  where :math:`\vec{p}_0` is any point in the plane.  Note with this definition the :math:`d` parameter defining the plane is :math:`d = -\vec{p}_0 \cdot \hat{n}`.

  .. cpp:member:: Vector2d normal

     The unit normal to the plane :math:`\hat{n}`.

  .. cpp:member:: double dist

     The minimum signed distance from the origin to the plane :math:`d`.

  .. cpp:member:: int ID

     An optional integer identification number for the plane.  This is used by Vertex2d to record which plane(s) are responsible for creating the vertex.

  .. cpp:function:: Plane2d()

     Default constructor -- implies :math:`(\hat{n}, d) = (1,0,0), 0.0)`, ID = std::numeric_limits<int>::min().

  .. cpp:function:: Plane2d(const double d, const Vector2d& nhat)

     Construct with :math:`(\hat{n}, d)` = (nhat, d), ID = std::numeric_limits<int>::min().

  .. cpp:function:: Plane2d(const Vector2d& p, const Vector2d& nhat)

     Construct with :math:`(\hat{n}, d)` = (nhat, :math:`-p\cdot\hat{n}`), ID = std::numeric_limits<int>::min().

  .. cpp:function:: Plane2d(const Vector2d& p, const Vector2d& nhat, const int id)

     Construct with :math:`(\hat{n}, d)` = (nhat, :math:`-p\cdot \hat{n}`), ID = id.

Vertex classes
--------------------
..
  ------------------------------------------------------------------------------
  Vertex2d
  ------------------------------------------------------------------------------

.. cpp:class:: Vertex2d

  Vertex2d is used to encode Polygons in 2d.  A vertex includes a position and the connectivity to neighboring vertices in the Polygon.  In this 2d case, the connectivity is always 2 vertices, ordered such that going from the first neighbor, to this vertex, and on to the last neighbor goes around the Polygon in the counter-clockwise direction.
