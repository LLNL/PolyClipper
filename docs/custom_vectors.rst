##################################################
Using PolyClipper with a user defined Vector type
##################################################

The examples shown in :ref:`Polyclipper concepts` all use PolyClipper's native Python interface, which employs PolyClipper's native (``Vector2d`` and ``Vector3d``) classes.  In C++ however the PolyClipper functions are all templated on a trait class describing how to invoke the assorted Vector operations (such as dot and cross products), allowing the C++ user to use PolyClipper methods directly with their own geometric Vector types.  The header file ``polyclipper_adapter.hh`` shows how the appropriate trait class for using the native PolyClipper Vector's is defined.  In this section we will describe how to write a trait class to adapt your own C++ Vector's, as well as show concrete examples adapting the STL class ``std::array`` for use as geometric Vectors in PolyClipper.

Writing a trait class for adapting a Vector with PolyClipper
------------------------------------------------------------

The following code is an example of all the methods we need to define in a trait class to adapt any C++ geometric Vector type for use with PolyClipper.  Depending on the user Vector type, the traited methods below may be simple wrappers around the corresponding methods of the users Vector type, or may entirely define the appropriate calculations.

.. code-block:: c++

    struct VectorAdapter {
      using VECTOR =                                             // Fill in your Vector type here
      static VECTOR Vector(double a, double b)                   // only 2D : construct and return a 2D Vector
      static VECTOR Vector(double a, double b, double c)         // only 3D : construct and return a 3D Vector
      static bool equal(const VECTOR& a, const VECTOR& b)        // test if a == b
      static double& x(VECTOR& a)                                // return the x component
      static double& y(VECTOR& a)                                // return the y component
      static double& z(VECTOR& a)                                // only 3D : return the z component
      static double  x(const VECTOR& a)                          // return the x component (const Vector)
      static double  y(const VECTOR& a)                          // return the y component (const Vector)
      static double  z(const VECTOR& a)                          // only 3D : return the z component (const Vector)
      static double  dot(const VECTOR& a, const VECTOR& b)       // return the dot product of a and b
      static double  crossmag(const VECTOR& a, const VECTOR& b)  // only 2D : return the z component of the cross product of a and b
      static VECTOR  cross(const VECTOR& a, const VECTOR& b)     // only 3D : return the cross product of a and b
      static double  magnitude2(const VECTOR& a)                 // return the square of the magnitude
      static double  magnitude(const VECTOR& a)                  // return the magnitude
      static VECTOR& imul(VECTOR& a, const double b)             // in place multiplication by a scalar
      static VECTOR& idiv(VECTOR& a, const double b)             // in place division by a scalar
      static VECTOR& iadd(VECTOR& a, const VECTOR& b)            // in place addition of Vectors a and b
      static VECTOR& isub(VECTOR& a, const VECTOR& b)            // in place subtraction of Vectors a and b
      static VECTOR  mul(const VECTOR& a, const double b)        // return the multiplication of Vector a by scalar b
      static VECTOR  div(const VECTOR& a, const double b)        // return the division of Vector a by scalar b
      static VECTOR  add(const VECTOR& a, const VECTOR& b)       // return the sum of a + b
      static VECTOR  sub(const VECTOR& a, const VECTOR& b)       // return the difference a - b
      static VECTOR  neg(const VECTOR& a)                        // return -a
      static VECTOR  unitVector(const VECTOR& a)                 // return a unit vector in the direction of a
      static std::string str(const VECTOR& a)                    // return a string representation of a
    };

.. note::
   These methods need only be provided for 2D vectors::

       VECTOR  VectorAdapter::Vector(double a, double b)
       double  VectorAdapter::crossmag(const VECTOR& a, const VECTOR& b)

.. note::
   These methods need only be provided for 3D vectors::

       VECTOR  VectorAdapter::Vector(double a, double b, double c)
       VECTOR  VectorAdapter::cross(const VECTOR& a, const VECTOR& b)
       double& VectorAdapter::z(VECTOR& a)
       double  VectorAdapter::z(const VECTOR& a)

Once you have your trait class for adapting your C++ Vector ready, you can call PolyClipper functions with this class a template parameter.  For instance, in 3D the method for clipping a polyhedron is templated with a default for PolyClipper's internal Vector type::

  template<typename VA = internal::VectorAdapter<Vector3d>>
  void clipPolyhedron(std::vector<Vertex3d<VA>>& poly,
                      const std::vector<Plane<VA>>& planes);

We can use our own templated Vector adapter (in this case called ``VectorAdapater``) by passing it like so::

  std::vector<PolyClipper::Vertex3d<VectorAdapter>> poly;
  std::vector<PolyClipper::Plane<VectorAdapter>> planes;
  // Fill in poly and planes at this point...
  clipPolyhedron(poly, planes);


A 2D example using ``std::array<double, 2>`` as a Vector
--------------------------------------------------------

In the test directory of PolyClipper (under ``test/test_array_vector``) is a complete example creating and clipping a polygon using ``std::array<double, 2>`` as the Vector type.  The source looks like the following::

  #include "polyclipper2d.hh"

  #include <array>
  #include <iostream>

  // Define a trait class for using a simple C style array of doubles as 2D Vector type for use in PolyClipper.
  struct ArrayAdapter2d {  
    using VECTOR = std::array<double, 2>;
    static VECTOR Vector(double a, double b)                   { return {a, b}; }                 // only 2D
    static bool equal(const VECTOR& a, const VECTOR& b)        { return (a[0] == b[0]) and (a[1] == b[1]); }
    static double& x(VECTOR& a)                                { return a[0]; }
    static double& y(VECTOR& a)                                { return a[1]; }
    static double  x(const VECTOR& a)                          { return a[0]; }
    static double  y(const VECTOR& a)                          { return a[1]; }
    static double  dot(const VECTOR& a, const VECTOR& b)       { return a[0]*b[0] + a[1]*b[1]; }
    static double  crossmag(const VECTOR& a, const VECTOR& b)  { return a[0]*b[1] - a[1]*b[0]; }   // only 2D
    static double  magnitude2(const VECTOR& a)                 { return a[0]*a[0] + a[1]*a[1]; }
    static double  magnitude(const VECTOR& a)                  { return std::sqrt(magnitude2(a)); }
    static VECTOR& imul(VECTOR& a, const double b)             { a[0] *= b; a[1] *= b; return a; }
    static VECTOR& idiv(VECTOR& a, const double b)             { a[0] /= b; a[1] /= b; return a; }
    static VECTOR& iadd(VECTOR& a, const VECTOR& b)            { a[0] += b[0]; a[1] += b[1]; return a; }
    static VECTOR& isub(VECTOR& a, const VECTOR& b)            { a[0] -= b[0]; a[1] -= b[1]; return a; }
    static VECTOR  mul(const VECTOR& a, const double b)        { return Vector(a[0] * b, a[1] * b); }
    static VECTOR  div(const VECTOR& a, const double b)        { return Vector(a[0] / b, a[1] / b); }
    static VECTOR  add(const VECTOR& a, const VECTOR& b)       { return Vector(a[0] + b[0], a[1] + b[1]); }
    static VECTOR  sub(const VECTOR& a, const VECTOR& b)       { return Vector(a[0] - b[0], a[1] - b[1]); }
    static VECTOR  neg(const VECTOR& a)                        { return Vector(-a[0], -a[1]); }
    static VECTOR  unitVector(const VECTOR& a)                 { auto mag = magnitude(a); return mag > 1.0e-15 ? div(a, mag) : Vector(1.0, 0.0); }
    static std::string str(const VECTOR& a)                    { std::ostringstream os; os << "(" << a[0] << " " << a[1] << ")"; return os.str(); }
  };

  int main() {

    using VA = ArrayAdapter2d;
    using Vector = VA::VECTOR;
    using Polygon = std::vector<PolyClipper::Vertex2d<VA>>;
    using Plane = PolyClipper::Plane<VA>;

    // Make a square
    //     3    2
    //     |----|
    //     |    |
    //     |----|
    //     0    1
    const std::vector<Vector> square_points = {VA::Vector(0.0,0.0), VA::Vector(10.0,0.0), VA::Vector(10.0,10.0), VA::Vector(0.0,10.0)};
    const std::vector<std::vector<int>> square_neighbors = {{3, 1}, {0, 2}, {1, 3}, {2, 0}};
    Polygon poly;
    double area;
    Vector cent;
    PolyClipper::initializePolygon(poly, square_points, square_neighbors);
    PolyClipper::moments(area, cent, poly);
    std::cout << "Initial polygon: " << polygon2string(poly) << std::endl
              << "Moments: " << area << " " << VA::str(cent) << std::endl << std::endl;

    // Clip by a couple of planes.
    const std::vector<Plane> planes = {Plane(VA::Vector(0.2, 0.2), VA::unitVector(VA::Vector(1.0,  1.0)), 10),
                                       Plane(VA::Vector(0.2, 0.8), VA::unitVector(VA::Vector(0.5, -1.0)), 20)};
    PolyClipper::clipPolygon(poly, planes);
    PolyClipper::moments(area, cent, poly);
    std::cout << "After clipping: " << polygon2string(poly) << std::endl
              << "Moments: " << area << " " << VA::str(cent) << std::endl;

    return 0;
  }
  
.. note::
   Note that the ``PolyClipper::Vector2d`` internal type defines its own versions of methods required by a Vector adapter trait class, so the adapter in that case (found in ``polyclipper_adapter.hh``) is only a thin wrapper around those methods.  In this example, however, the ``std::array`` does not have any notion of mathematical concepts like the Vector dot product, cross product, magnitude, etc., so we explicitly compute those in the appropriate methods of ``ArrayAdapter2d`` above.


A 3D example using ``std::array<double, 3>`` as a Vector
--------------------------------------------------------

Similarly to above we can also use ``std::array<double, 3>`` as a Vector type in PolyClipper with an appropriate trait class::

  #include "polyclipper3d.hh"

  #include <array>
  #include <iostream>

  // Define a trait class for using a simple C style array of doubles as 3D Vector type for use in PolyClipper.
  struct ArrayAdapter3d {  
    using VECTOR = std::array<double, 3>;
    static VECTOR Vector(double a, double b, double c)         { return {a, b, c}; }               // only 3D
    static bool equal(const VECTOR& a, const VECTOR& b)        { return (a[0] == b[0]) and (a[1] == b[1]) and (a[2] == b[2]); }
    static double& x(VECTOR& a)                                { return a[0]; }
    static double& y(VECTOR& a)                                { return a[1]; }
    static double& z(VECTOR& a)                                { return a[2]; }
    static double  x(const VECTOR& a)                          { return a[0]; }
    static double  y(const VECTOR& a)                          { return a[1]; }
    static double  z(const VECTOR& a)                          { return a[2]; }
    static double  dot(const VECTOR& a, const VECTOR& b)       { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
    static VECTOR  cross(const VECTOR& a, const VECTOR& b)     { return {a[1]*b[2] - a[2]*b[1],    // only 3D
                                                                         a[2]*b[0] - a[0]*b[2],
                                                                         a[0]*b[1] - a[1]*b[0]}; }
    static double  magnitude2(const VECTOR& a)                 { return a[0]*a[0] + a[1]*a[1] + a[2]*a[2]; }
    static double  magnitude(const VECTOR& a)                  { return std::sqrt(magnitude2(a)); }
    static VECTOR& imul(VECTOR& a, const double b)             { a[0] *= b; a[1] *= b; a[2] *= b; return a; }
    static VECTOR& idiv(VECTOR& a, const double b)             { a[0] /= b; a[1] /= b; a[2] /= b; return a; }
    static VECTOR& iadd(VECTOR& a, const VECTOR& b)            { a[0] += b[0]; a[1] += b[1]; a[2] += b[2]; return a; }
    static VECTOR& isub(VECTOR& a, const VECTOR& b)            { a[0] -= b[0]; a[1] -= b[1]; a[2] -= b[2]; return a; }
    static VECTOR  mul(const VECTOR& a, const double b)        { return Vector(a[0] * b, a[1] * b, a[2] * b); }
    static VECTOR  div(const VECTOR& a, const double b)        { return Vector(a[0] / b, a[1] / b, a[2] / b); }
    static VECTOR  add(const VECTOR& a, const VECTOR& b)       { return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]); }
    static VECTOR  sub(const VECTOR& a, const VECTOR& b)       { return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]); }
    static VECTOR  neg(const VECTOR& a)                        { return Vector(-a[0], -a[1], -a[2]); }
    static VECTOR  unitVector(const VECTOR& a)                 { auto mag = magnitude(a); return mag > 1.0e-15 ? div(a, mag) : Vector(1.0, 0.0, 0.0); }
    static std::string str(const VECTOR& a)                    { std::ostringstream os; os << "(" << a[0] << " " << a[1] << " " << a[2] << ")"; return os.str(); }
  };

  int main() {

    using VA = ArrayAdapter3d;
    using Vector = VA::VECTOR;
    using Polyhedron = std::vector<PolyClipper::Vertex3d<VA>>;
    using Plane = PolyClipper::Plane<VA>;

    // Make a cube                   |y     
    //                               |      
    //                               |____x 
    //   3/-----/2                  /       
    //   /     /|                  /z       
    // 7|-----|6|
    //  |     | |
    //  |  0  | /1
    //  |_____|/
    //  4     5
    //
    const std::vector<Vector> cube_points = {VA::Vector(0,0,0),  VA::Vector(10,0,0),  VA::Vector(10,10,0),  VA::Vector(0,10,0),
                                             VA::Vector(0,0,10), VA::Vector(10,0,10), VA::Vector(10,10,10), VA::Vector(0,10,10)};
    const std::vector<std::vector<int>> cube_neighbors = {{1, 4, 3},
                                                          {5, 0, 2},
                                                          {3, 6, 1},
                                                          {7, 2, 0},
                                                          {5, 7, 0},
                                                          {1, 6, 4},
                                                          {5, 2, 7},
                                                          {4, 6, 3}};
    const std::vector<std::vector<int>> cube_facets = {{4, 5, 6, 7},
                                                       {1, 2, 6, 5},
                                                       {0, 3, 2, 1},
                                                       {4, 7, 3, 0},
                                                       {6, 2, 3, 7},
                                                       {1, 5, 4, 0}};
    Polyhedron poly;
    double vol;
    Vector cent;
    PolyClipper::initializePolyhedron(poly, cube_points, cube_neighbors);
    PolyClipper::moments(vol, cent, poly);
    std::cout << "Initial polyhedron: " << polyhedron2string(poly) << std::endl
              << "Moments: " << vol << " " << VA::str(cent) << std::endl << std::endl;

    // Clip by a couple of planes.
    const std::vector<Plane> planes = {Plane(VA::Vector(0.2, 0.2, 0.2), VA::unitVector(VA::Vector(1.0,  1.0,  1.0)), 10),
                                       Plane(VA::Vector(0.2, 0.8, 0.8), VA::unitVector(VA::Vector(0.5, -1.0, -1.0)), 20)};
    PolyClipper::clipPolyhedron(poly, planes);
    PolyClipper::moments(vol, cent, poly);
    std::cout << "After clipping: " << polyhedron2string(poly) << std::endl
              << "Moments: " << vol << " " << VA::str(cent) << std::endl;

    return 0;
  }


