//---------------------------------PolyClipper--------------------------------//
// Clip a faceted volume (polygon or polyhedron) by a set of planes in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane is retained.  Above here is taken to mean the direction the plane
// normal points in.  This can also be defined in terms of a signed distance.
//
// We can define a plane by either a unit normal and point in the plane (n, p0),
// or the plane normal and distance from the plane to the origin (n, d0).  The
// signed distance of a point (p) from the plane is then given by
//
//   d_s = (p - p0).dot(n) = d0 + p.dot(n)
//
// The the volume with d_s > 0 is above the plane, while d_s < 0 is below.
//
// The algorithms herein are based on R3D as outlined in 
// Powell, D., & Abel, T. (2015). An exact general remeshing scheme applied to 
// physically conservative voxelization. Journal of Computational Physics, 297, 340–356.
//
// Created by J. Michael Owen, Tue Nov 28 10:00:51 PST 2017
//----------------------------------------------------------------------------//
#ifndef __PolyClipper_hh__
#define __PolyClipper_hh__

#include <cmath>
#include <string>
#include <limits>
#include <vector>
#include <set>

namespace PolyClipper {

//------------------------------------------------------------------------------
// Trait class to use PolyClipper's native Vector types.
// This is an example of the traits necessary to use your own Vector types
// directly with PolyClipper.
//------------------------------------------------------------------------------
namespace internal {
template<typename VectorType>
struct VectorAdapter {
  using VECTOR = VectorType;
  static VECTOR Vector(double a, double b)                 { return VECTOR(a, b); }    // only 2D
  static VECTOR Vector(double a, double b, double c)       { return VECTOR(a, b, c); } // only 3D
  static bool equal(const VECTOR& a, const VECTOR& b)      { return a == b; }
  static double& x(VECTOR& a)                              { return a.x; }
  static double& y(VECTOR& a)                              { return a.y; }
  static double& z(VECTOR& a)                              { return a.z; }
  static double dot(const VECTOR& a, const VECTOR& b)      { return a.dot(b); }
  static double crossmag(const VECTOR& a, const VECTOR& b) { return a.crossmag(b); }   // only 2D
  static VECTOR cross(const VECTOR& a, const VECTOR& b)    { return a.cross(b); }      // only 3D
  static double magnitude2(const VECTOR& a)                { return a.magnitude2(); }
  static double magnitude(const VECTOR& a)                 { return a.magnitude(); }
  static void imul(const VECTOR& a, const double b)        { a *= b; }
  static void idiv(const VECTOR& a, const double b)        { a /= b; }
  static void iadd(const VECTOR& a, const VECTOR& b)       { a += b; }
  static void isub(const VECTOR& a, const VECTOR& b)       { a -= b; }
  static VECTOR mul(const VECTOR& a, const double b)       { return a * b; }
  static VECTOR div(const VECTOR& a, const double b)       { return a / b; }
  static VECTOR add(const VECTOR& a, const VECTOR& b)      { return a + b; }
  static VECTOR sub(const VECTOR& a, const VECTOR& b)      { return a - b; }
  static VECTOR neg(const VECTOR& a)                       { return -a; }
  static VECTOR unitVector(const VECTOR& a)                { return a.unitVector(); }
};
}

//------------------------------------------------------------------------------
// 2D Vector
//------------------------------------------------------------------------------
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

//------------------------------------------------------------------------------
// 3D Vector
//------------------------------------------------------------------------------
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
};

inline Vector3d operator*(const double lhs, const Vector3d& rhs) { return rhs*lhs; }

//------------------------------------------------------------------------------
// Plane.
//------------------------------------------------------------------------------
template<typename VectorType,
         typename VA = internal::VectorAdapter<VectorType>>
struct Plane {
  using Vector = VectorType;
  double dist;                       // Signed distance to the origin
  Vector normal;                     // Unit normal
  int ID;                            // ID for the plane, used to label vertices
  Plane()                                                  : dist(0.0), normal(VA::VECTOR(1.0, 0.0, 0.0)), ID(std::numeric_limits<int>::min()) {}
  Plane(const double d, const Vector& nhat)                : dist(d), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane(const Vector& p, const Vector& nhat)               : dist(VA::dot(-p, nhat)), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane(const Vector& p, const Vector& nhat, const int id) : dist(VA::dot(-p, nhat)), normal(nhat), ID(id) {}
  Plane(const Plane2d& rhs)                                : dist(rhs.dist), normal(rhs.normal), ID(rhs.ID) {}
  Plane& operator=(const Plane& rhs)                       { dist = rhs.dist; normal = rhs.normal; ID = rhs.ID; return *this; }
  bool operator==(const Plane& rhs) const                  { return (dist == rhs.dist and VA::equal(normal, rhs.normal)); }
  bool operator!=(const Plane& rhs) const                  { return not (*this == rhs); }
  bool operator< (const Plane& rhs) const                  { return (dist < rhs.dist); }
  bool operator> (const Plane& rhs) const                  { return (dist > rhs.dist); }
};

//------------------------------------------------------------------------------
// The 2D vertex struct, which we use to encode polygons.
//------------------------------------------------------------------------------
template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
struct Vertex2d {
  using Vector = VectorType;
  Vector position;
  std::pair<int, int> neighbors;
  int comp;
  mutable int ID;                            // convenient, but sneaky
  mutable std::set<int> clips;               // the planes (if any) that created this point
  Vertex2d()                                 : position(VA::VECTOR(0.0,0.0)), neighbors(), comp(1), ID(-1), clips() {}
  Vertex2d(const Vector& pos)                : position(pos),                 neighbors(), comp(1), ID(-1), clips() {}
  Vertex2d(const Vector& pos, const int c)   : position(pos),                 neighbors(), comp(c), ID(-1), clips() {}
  Vertex2d(const Vertex2d& rhs)              : position(rhs.position), neighbors(rhs.neighbors), comp(rhs.comp), ID(rhs.ID), clips(rhs.clips) {}
  Vertex2d& operator=(const Vertex2d& rhs)   { position = rhs.position; neighbors = rhs.neighbors; comp = rhs.comp; ID = rhs.ID; clips = rhs.clips; return *this; }
  bool operator==(const Vertex2d& rhs) const {
    return (VA::equal(position, rhs.position) and
            neighbors == rhs.neighbors and
            comp      == rhs.comp and
            ID        == rhs.ID);
  }
};

//------------------------------------------------------------------------------
// The 3D vertex struct, which we use to encode polyhedra.
//------------------------------------------------------------------------------
template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
struct Vertex3d {
  using Vector = VectorType;
  Vector position;
  std::vector<int> neighbors;
  int comp;
  mutable int ID;                            // convenient, but sneaky
  mutable std::set<int> clips;               // the planes (if any) that created this point
  Vertex3d()                                 : position(VA::VECTOR(0.0, 0.0, 0.0)), neighbors(), comp(1), ID(-1), clips() {}
  Vertex3d(const Vector& pos)                : position(pos),                       neighbors(), comp(1), ID(-1), clips() {}
  Vertex3d(const Vector& pos, const int c)   : position(pos),                       neighbors(), comp(c), ID(-1), clips() {}
  Vertex3d(const Vertex3d& rhs)              : position(rhs.position), neighbors(rhs.neighbors), comp(rhs.comp), ID(rhs.ID), clips(rhs.clips) {}
  Vertex3d& operator=(const Vertex3d& rhs)   { position = rhs.position; neighbors = rhs.neighbors; comp = rhs.comp; ID = rhs.ID; clips = rhs.clips; return *this; }
  bool operator==(const Vertex3d& rhs) const {
    return (VA::equal(position, rhs.position) and
            neighbors == rhs.neighbors and
            comp      == rhs.comp and
            ID        == rhs.ID);
  }
};

//------------------------------------------------------------------------------
// 2D (polygon) methods.
//------------------------------------------------------------------------------
template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
void initializePolygon(std::vector<Vertex2d<VectorType, VA>& poly,
                       const std::vector<VectorType>& positions,
                       const std::vector<std::vector<int>>& neighbors);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
std::string polygon2string(const std::vector<Vertex2d<VectorType, VA>>& poly);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
void moments(double& zerothMoment, VectorType& firstMoment,
             const std::vector<Vertex2d<VectorType, VA>>& polygon);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
void clipPolygon(std::vector<Vertex2d<VectorType, VA>>& poly,
                 const std::vector<Plane<VertexType, VA>>& planes);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
void collapseDegenerates(std::vector<Vertex2d<VectorType, VA>>& poly,
                         const double tol);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
std::vector<std::vector<int>> extractFaces(const std::vector<Vertex2d<VectorType, VA>>& poly);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
std::vector<std::set<int>> commonFaceClips(const std::vector<Vertex2d<VectorType, VA>>& poly,
                                           const std::vector<std::vector<int>>& faces);

template<typename VectorType = Vector2d,
         typename VA = internal::VectorAdapter<Vector2d>>
std::vector<std::vector<int>> splitIntoTriangles(const std::vector<Vertex2d<VectorType, VA>>& poly,
                                                 const double tol = 0.0);

// //------------------------------------------------------------------------------
// // 3D (polyhedron) methods.
// //------------------------------------------------------------------------------
// typedef std::vector<Vertex3d> Polyhedron;

// void initializePolyhedron(Polyhedron& poly,
//                           const std::vector<PolyClipper::Vector3d>& positions,
//                           const std::vector<std::vector<int>>& neighbors);

// std::string polyhedron2string(const Polyhedron& poly);

// void moments(double& zerothMoment, PolyClipper::Vector3d& firstMoment,
//               const Polyhedron& polyhedron);

// void clipPolyhedron(Polyhedron& poly,
//                      const std::vector<Plane3d>& planes);

// void collapseDegenerates(Polyhedron& poly,
//                            const double tol);

// std::vector<std::vector<int>> extractFaces(const Polyhedron& poly);

// std::vector<std::set<int>> commonFaceClips(const Polyhedron& poly,
//                                            const std::vector<std::vector<int>>& faces);

// std::vector<std::vector<int>> splitIntoTetrahedra(const Polyhedron& poly, 
//                                                   const double tol = 0.0);


}

#endif

