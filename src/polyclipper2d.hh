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
#ifndef __PolyClipper2d__
#define __PolyClipper2d__

#include "polyclipper_adapter.hh"
#include "polyclipper_plane.hh"
#include "polyclipper_vector2d.hh"
#include "polyclipper_utilities.hh"

#include <cmath>
#include <string>
#include <limits>
#include <vector>
#include <set>

namespace PolyClipper {

// A convenient alias for the 2D plane with PolyClippers internal Vector type
using Plane2d = Plane<internal::VectorAdapter<Vector2d>>;

//------------------------------------------------------------------------------
// The 2D vertex struct, which we use to encode polygons, i.e., 
// polygons are specified as std::vector<Vertex2d>.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
struct Vertex2d {
  using Vector = typename VA::VECTOR;
  Vector position;
  std::pair<int, int> neighbors;
  int comp;
  mutable int ID;                            // convenient, but sneaky
  mutable std::set<int> clips;               // the planes (if any) that created this point
  Vertex2d()                                 : position(VA::Vector(0.0,0.0)), neighbors(), comp(1), ID(-1), clips() {}
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
// Initialize a polygon given the vertex coordinates and connectivity.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
void initializePolygon(std::vector<Vertex2d<VA>>& poly,
                       const std::vector<typename VA::VECTOR>& positions,
                       const std::vector<std::vector<int>>& neighbors);

//------------------------------------------------------------------------------
// Return a nicely formatted string representing the polygon.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
std::string polygon2string(const std::vector<Vertex2d<VA>>& poly);

//------------------------------------------------------------------------------
// Compute the zeroth and first moment of a Polygon.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
void moments(double& zerothMoment, typename VA::VECTOR& firstMoment,
             const std::vector<Vertex2d<VA>>& polygon);

//------------------------------------------------------------------------------
// Clip a polygon by planes.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
void clipPolygon(std::vector<Vertex2d<VA>>& poly,
                 const std::vector<Plane<VA>>& planes);

//------------------------------------------------------------------------------
// Collapse degenerate vertices.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
void collapseDegenerates(std::vector<Vertex2d<VA>>& poly,
                         const double tol);

//------------------------------------------------------------------------------
// Return the vertices ordered in faces.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
std::vector<std::vector<int>> extractFaces(const std::vector<Vertex2d<VA>>& poly);

//------------------------------------------------------------------------------
// Compute the set of clips common to each face.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
std::vector<std::set<int>> commonFaceClips(const std::vector<Vertex2d<VA>>& poly,
                                           const std::vector<std::vector<int>>& faces);

//------------------------------------------------------------------------------
// Split a polygon into a set of triangles.
//------------------------------------------------------------------------------
template<typename VA = internal::VectorAdapter<Vector2d>>
std::vector<std::vector<int>> splitIntoTriangles(const std::vector<Vertex2d<VA>>& poly,
                                                 const double tol = 0.0);

}

#include "polyclipper2dImpl.hh"

#endif

