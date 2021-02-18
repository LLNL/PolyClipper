//------------------------------------------------------------------------------
// Plane.
//------------------------------------------------------------------------------
#ifndef __PolyClipper_Plane__
#define __PolyClipper_Plane__

namespace PolyClipper {

template<typename VectorType,
         typename VA = internal::VectorAdapter<VectorType>>
struct Plane {
  using Vector = VectorType;
  double dist;                       // Signed distance to the origin
  Vector normal;                     // Unit normal
  int ID;                            // ID for the plane, used to label vertices
  Plane()                                                  : dist(0.0), normal(VA::Vector(1.0, 0.0, 0.0)), ID(std::numeric_limits<int>::min()) {}
  Plane(const double d, const Vector& nhat)                : dist(d), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane(const Vector& p, const Vector& nhat)               : dist(VA::dot(-p, nhat)), normal(nhat), ID(std::numeric_limits<int>::min()) {}
  Plane(const Vector& p, const Vector& nhat, const int id) : dist(VA::dot(-p, nhat)), normal(nhat), ID(id) {}
  Plane(const Plane& rhs)                                  : dist(rhs.dist), normal(rhs.normal), ID(rhs.ID) {}
  Plane& operator=(const Plane& rhs)                       { dist = rhs.dist; normal = rhs.normal; ID = rhs.ID; return *this; }
  bool operator==(const Plane& rhs) const                  { return (dist == rhs.dist and VA::equal(normal, rhs.normal)); }
  bool operator!=(const Plane& rhs) const                  { return not (*this == rhs); }
  bool operator< (const Plane& rhs) const                  { return (dist < rhs.dist); }
  bool operator> (const Plane& rhs) const                  { return (dist > rhs.dist); }
};

}

#endif
