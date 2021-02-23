//------------------------------------------------------------------------------
// 2D Vector
//------------------------------------------------------------------------------
#ifndef __PolyClipper_Vector2d__
#define __PolyClipper_Vector2d__

#include <cmath>

namespace PolyClipper {

struct Vector2d {
  double x, y;
  Vector2d(): x(0.0), y(0.0) {}
  Vector2d(double X, double Y, double Z = 0.0)    : x(X), y(Y) {}  // Z dummied out for interface consistency with 3D
  bool      operator==(const Vector2d& rhs) const { return x == rhs.x and y == rhs.y; }
  double    dot(const Vector2d& rhs) const        { return x*rhs.x + y*rhs.y; }
  double    crossmag(const Vector2d& rhs) const   { return x*rhs.y - y*rhs.x; }
  double    magnitude2() const                    { return x*x + y*y; }
  double    magnitude() const                     { return std::sqrt(x*x + y*y); }
  Vector2d& operator*=(const double rhs)          { x *= rhs; y *= rhs; return *this; }
  Vector2d& operator/=(const double rhs)          { x /= rhs; y /= rhs; return *this; }
  Vector2d  operator+=(const Vector2d rhs)        { x += rhs.x; y += rhs.y; return *this; }
  Vector2d  operator-=(const Vector2d rhs)        { x -= rhs.x; y -= rhs.y; return *this; }
  Vector2d  operator*(const double rhs) const     { return Vector2d(rhs*x, rhs*y); }
  Vector2d  operator/(const double rhs) const     { return Vector2d(x/rhs, y/rhs); }
  Vector2d  operator+(const Vector2d rhs) const   { return Vector2d(x + rhs.x, y + rhs.y); }
  Vector2d  operator-(const Vector2d rhs) const   { return Vector2d(x - rhs.x, y - rhs.y); }
  Vector2d  operator-() const                     { return Vector2d(-x, -y); }
  Vector2d  unitVector() const {
    const auto mag = this->magnitude();
    return (mag > 0.0 ? Vector2d(x/mag, y/mag) : Vector2d(1.0, 0.0));
  }
};

inline Vector2d operator*(const double lhs, const Vector2d& rhs) { return rhs*lhs; }

}

#endif
