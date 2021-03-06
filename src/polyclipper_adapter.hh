//------------------------------------------------------------------------------
// Trait class to use PolyClipper's native Vector types.
// This is an example of the traits necessary to use your own Vector types
// directly with PolyClipper.
//------------------------------------------------------------------------------
#ifndef __PolyClipper_adapter__
#define __PolyClipper_adapter__

#include "polyclipper_utilities.hh"

#include <string>
#include <sstream>
#include <array>

namespace PolyClipper {
namespace internal {

template<typename VectorType>
struct VectorAdapter {
  using VECTOR = VectorType;
  static VECTOR  Vector(double a, double b)                  { return VECTOR(a, b); }    // only 2D
  static VECTOR  Vector(double a, double b, double c)        { return VECTOR(a, b, c); } // only 3D
  static bool    equal(const VECTOR& a, const VECTOR& b)     { return a == b; }
  static double& x(VECTOR& a)                                { return a.x; }
  static double& y(VECTOR& a)                                { return a.y; }
  static double& z(VECTOR& a)                                { return a.z; }             // only 3D
  static double  x(const VECTOR& a)                          { return a.x; }
  static double  y(const VECTOR& a)                          { return a.y; }
  static double  z(const VECTOR& a)                          { return a.z; }             // only 3D
  static double  dot(const VECTOR& a, const VECTOR& b)       { return a.dot(b); }
  static double  crossmag(const VECTOR& a, const VECTOR& b)  { return a.crossmag(b); }   // only 2D
  static VECTOR  cross(const VECTOR& a, const VECTOR& b)     { return a.cross(b); }      // only 3D
  static double  magnitude2(const VECTOR& a)                 { return a.magnitude2(); }
  static double  magnitude(const VECTOR& a)                  { return a.magnitude(); }
  static VECTOR& imul(VECTOR& a, const double b)             { a *= b; return a; }
  static VECTOR& idiv(VECTOR& a, const double b)             { a /= b; return a; }
  static VECTOR& iadd(VECTOR& a, const VECTOR& b)            { a += b; return a; }
  static VECTOR& isub(VECTOR& a, const VECTOR& b)            { a -= b; return a; }
  static VECTOR  mul(const VECTOR& a, const double b)        { return a * b; }
  static VECTOR  div(const VECTOR& a, const double b)        { return a / b; }
  static VECTOR  add(const VECTOR& a, const VECTOR& b)       { return a + b; }
  static VECTOR  sub(const VECTOR& a, const VECTOR& b)       { return a - b; }
  static VECTOR  neg(const VECTOR& a)                        { return -a; }
  static VECTOR  unitVector(const VECTOR& a)                 { return a.unitVector(); }
  static std::string str(const VECTOR& a)                    { std::ostringstream os; os << a; return os.str(); }  // Pretty human readable representation
  static std::array<double, 3> get_triple(const VECTOR& a)   { return a.triple(); }
  static void set_triple(VECTOR& a,
                         const std::array<double, 3>& vals)  { a.set_triple(vals); }
};

}
}

#endif
