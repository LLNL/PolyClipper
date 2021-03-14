//------------------------------------------------------------------------------
// 3D Vector
//------------------------------------------------------------------------------
#ifndef __PolyClipper_Vector3d__
#define __PolyClipper_Vector3d__

#include <cmath>
#include <vector>
#include <array>

namespace PolyClipper {

struct Vector3d {
  double x, y, z;
  Vector3d(): x(0.0), y(0.0), z(0.0) {}
  Vector3d(double X, double Y, double Z): x(X), y(Y), z(Z) {}
  bool      operator==(const Vector3d& rhs) const { return x == rhs.x and y == rhs.y and z == rhs.z; }
  double    dot(const Vector3d& rhs) const        { return x*rhs.x + y*rhs.y + z*rhs.z; }
  Vector3d  cross(const Vector3d& rhs) const      { return Vector3d(y*rhs.z - z*rhs.y,
                                                                    z*rhs.x - x*rhs.z,
                                                                    x*rhs.y - y*rhs.x); }
  double    magnitude2() const                    { return x*x + y*y + z*z; }
  double    magnitude() const                     { return std::sqrt(x*x + y*y + z*z); }
  Vector3d& operator*=(const double rhs)          { x *= rhs; y *= rhs; z *= rhs; return *this; }
  Vector3d& operator/=(const double rhs)          { x /= rhs; y /= rhs; z /= rhs; return *this; }
  Vector3d  operator+=(const Vector3d rhs)        { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
  Vector3d  operator-=(const Vector3d rhs)        { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }
  Vector3d  operator*(const double rhs) const     { return Vector3d(rhs*x, rhs*y, rhs*z); }
  Vector3d  operator/(const double rhs) const     { return Vector3d(x/rhs, y/rhs, z/rhs); }
  Vector3d  operator+(const Vector3d rhs) const   { return Vector3d(x + rhs.x, y + rhs.y, z + rhs.z); }
  Vector3d  operator-(const Vector3d rhs) const   { return Vector3d(x - rhs.x, y - rhs.y, z - rhs.z); }
  Vector3d  operator-() const                     { return Vector3d(-x, -y, -z); }
  Vector3d  unitVector() const {
    const auto mag = this->magnitude();
    return (mag > 0.0 ? Vector3d(x/mag, y/mag, z/mag) : Vector3d(1.0, 0.0, 0.0));
  }
  std::array<double, 3> triple() const               { return {x, y, z}; }
  void set_triple(const std::array<double, 3>& vals) { x = vals[0]; y = vals[1]; z = vals[2]; }
};

inline Vector3d operator*(const double lhs, const Vector3d& rhs) { return rhs*lhs; }

}

#endif

