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
  static std::array<double, 3> get_triple(const VECTOR& a)   { return a; }
  static void set_triple(VECTOR& a,
                         const std::array<double, 3> vals)   { a = vals; }
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

  
