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

  
