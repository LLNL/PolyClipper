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

.. cpp:class:: Vector2d

  Vector2d represents a 2D vector coordinate in space, :math:`(x,y)`.  It supports a small number of simple vector manipulation operations:

  .. cpp:member:: double x

     The :math:`x` coordinate

  .. cpp:member:: double y

     The :math:`y` coordinate

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

  .. cpp:function:: operator-() const

     :return: negative of this vector: :math:`-\vec{a}`
     :rtype: Vector2d

  .. cpp:function:: unitVector() const

     :rtype: Vector2d

     Returns the unit vector pointing in the direction of this one: :math:`\vec{a}/\sqrt{\vec{a} \cdot \vec{a}}`.

     If :math:`\vec{a} = (0,0)`, returns the unit vector in the :math:`x` direction: :math:`(1,0)`.
     
