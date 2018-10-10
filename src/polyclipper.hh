//---------------------------------PolyClipper--------------------------------//
// Clip a faceted volume (polygon or polyhedron) by a set of planes in place.
//
// We use the convention that any portion of the faceted volume "below" the 
// plane is clipped, i.e., only the portion of the faceted volume "above" the 
// plane
//    plane.compare(point) >= 0
// is retained.
//
// The algorithms herein are based on R3D as outlined in 
// Powell, D., & Abel, T. (2015). An exact general remeshing scheme applied to 
// physically conservative voxelization. Journal of Computational Physics, 297, 340–356.
//
// Created by J. Michael Owen, Tue Nov 28 10:00:51 PST 2017
//----------------------------------------------------------------------------//
#ifndef __PolyClipper_hh__
#define __PolyClipper_hh__

#include "Geometry/Dimension.hh"
#include "Geometry/GeomPlane.hh"

#include <string>
#include <vector>

namespace PolyClipper {

//------------------------------------------------------------------------------
// 2D Vector
//------------------------------------------------------------------------------
struct Vector2d {
  double x, y;
  Vector(double X, double Y): x(X), y(Y) {}
  bool operator==(const Vector2d& rhs) { return x == rhs.x and y == rhs.y; }
};

//------------------------------------------------------------------------------
// 3D Vector
//------------------------------------------------------------------------------
struct Vector3d {
  double x, y, z;
  Vector(double X, double Y, double Z): x(X), y(Y), z(Z) {}
  bool operator==(const Vector3d& rhs) { return x == rhs.x and y == rhs.y and z == rhs.z; }
};

//------------------------------------------------------------------------------
// 2D plane.
//------------------------------------------------------------------------------
struct Plane2d {
  typedef PolyClipper::Vector2d Vector;
  double dist;                       // Signed distance to the origin
  Vector normal;                     // Unit normal
  Plane2d(): dist(0.0), normal(1,0) {}
  Plane2d(const double d, const Vector& nhat): dist(d), normal(nhat) {}
  Plane2d(const Vector& p, const Vector& nhat): dist(-p.dot(nhat)), normal(nhat) {}
  Plane2d& operator=(const Plane2d& rhs) { dist = rhs.dist; normal = rhs.normal; return *this; }
  bool operator==(const Plane2d& rhs) { return (dist == rhs.dist and normal == rhs.normal); }
};

//------------------------------------------------------------------------------
// 3D plane.
//------------------------------------------------------------------------------
struct Plane3d {
  typedef PolyClipper::Vector3d Vector;
  Vector normal;                     // Unit normal
  double dist;                       // Signed distance to the origin
  Plane3d(): dist(0.0), normal(1,0,0) {}
  Plane3d(const double d, const Vector& nhat): dist(d), normal(nhat) {}
  Plane3d(const Vector& p, const Vector& nhat): dist(-p.dot(nhat)), normal(nhat) {}
  Plane3d& operator=(const Plane3d& rhs) { dist = rhs.dist; normal = rhs.normal; return *this; }
  bool operator==(const Plane3d& rhs) { return (dist == rhs.dist and normal == rhs.normal); }
};

//------------------------------------------------------------------------------
// The 2D vertex struct, which we use to encode polygons.
//------------------------------------------------------------------------------
struct Vertex2d {
  typedef PolyClipper::Vector2d Vector;
  Vector position;
  std::pair<int, int> neighbors;
  int comp;
  mutable int ID;    // convenient, but sneaky
  Vertex2d():                               position(),    neighbors(), comp(1), ID(-1) {}
  Vertex2d(const Vector& pos):              position(pos), neighbors(), comp(1), ID(-1) {}
  Vertex2d(const Vector& pos, const int c): position(pos), neighbors(), comp(c), ID(-1) {}
  bool operator==(const Vertex2d& rhs) const {
    return (position  == rhs.position and
            neighbors == rhs.neighbors and
            comp      == rhs.comp and
            ID        == rhs.ID);
  }
};

//------------------------------------------------------------------------------
// The 3D vertex struct, which we use to encode polyhedra.
//------------------------------------------------------------------------------
struct Vertex3d {
  typedef PolyClipper::Vector3d Vector;
  Vector position;
  std::vector<int> neighbors;
  int comp;
  mutable int ID;    // convenient, but sneaky
  Vertex3d():                               position(),    neighbors(), comp(1), ID(-1) {}
  Vertex3d(const Vector& pos):              position(pos), neighbors(), comp(1), ID(-1) {}
  Vertex3d(const Vector& pos, const int c): position(pos), neighbors(), comp(c), ID(-1) {}
  bool operator==(const Vertex3d& rhs) const {
    return (position  == rhs.position and
            neighbors == rhs.neighbors and
            comp      == rhs.comp and
            ID        == rhs.ID);
  }
};

//------------------------------------------------------------------------------
// 2D (polygon) methods.
//------------------------------------------------------------------------------
typedef std::vector<Vertex2d> Polygon;

void initializePolygon(Polygon& poly,
                       const std::vector<PolyClipper::Vector2d>& positions,
                       const std::vector<std::vector<int>>& neighbors);

std::string polygon2string(const Polygon& poly);

void moments(double& zerothMoment, PolyClipper::Vector2d& firstMoment,
             const Polygon& polygon);

void clipPolygon(Polygon& poly,
                 const std::vector<Plane2d>& planes);

void collapseDegenerates(Polygon& poly,
                         const double tol);

  std::vector<std::vector<int>> splitIntoTriangles(const Polygon& poly,
                                                   const double tol = 0.0);

//------------------------------------------------------------------------------
// 3D (polyhedron) methods.
//------------------------------------------------------------------------------
typedef std::vector<Vertex3d> Polyhedron;

std::vector<std::vector<int>> extractFaces(const Polyhedron& poly);

void initializePolyhedron(Polyhedron& poly,
                          const std::vector<PolyClipper::Vector3d>& positions,
                          const std::vector<std::vector<int>>& neighbors);

std::string polyhedron2string(const Polyhedron& poly);

void moments(double& zerothMoment, PolyClipper::Vector3d& firstMoment,
             const Polyhedron& polyhedron);

void clipPolyhedron(Polyhedron& poly,
                    const std::vector<Plane3d>& planes);

void collapseDegenerates(Polyhedron& poly,
                         const double tol);

std::vector<std::vector<int>> splitIntoTetrahedra(const Polyhedron& poly, 
                                                  const double tol = 0.0);
}

#endif

