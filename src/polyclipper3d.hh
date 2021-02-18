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
#ifndef __PolyClipper3d_hh__
#define __PolyClipper3d_hh__

#include "polyclipper_adapter.hh"
#include "polyclipper_plane.hh"
#include "polyclipper_vector3d.hh"
#include "polyclipper_utilities.hh"

#include <cmath>
#include <string>
#include <limits>
#include <vector>
#include <set>

namespace PolyClipper {

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
  Vertex3d()                                 : position(VA::Vector(0.0, 0.0, 0.0)), neighbors(), comp(1), ID(-1), clips() {}
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
// 3D (polyhedron) methods.
//------------------------------------------------------------------------------
template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
void initializePolyhedron(std::vector<Vertex3d<VectorType, VA>>& poly,
                          const std::vector<VectorType>& positions,
                          const std::vector<std::vector<int>>& neighbors);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
std::string polyhedron2string(const std::vector<Vertex3d<VectorType, VA>>& poly);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
void moments(double& zerothMoment, VectorType& firstMoment,
             const std::vector<Vertex3d<VectorType, VA>>& polyhedron);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
void clipPolyhedron(std::vector<Vertex3d<VectorType, VA>>& poly,
                    const std::vector<Plane<VectorType, VA>>& planes);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
void collapseDegenerates(std::vector<Vertex3d<VectorType, VA>>& poly,
                         const double tol);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
std::vector<std::vector<int>> extractFaces(const std::vector<Vertex3d<VectorType, VA>>& poly);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
std::vector<std::set<int>> commonFaceClips(const std::vector<Vertex3d<VectorType, VA>>& poly,
                                           const std::vector<std::vector<int>>& faces);

template<typename VectorType = Vector3d,
         typename VA = internal::VectorAdapter<Vector3d>>
std::vector<std::vector<int>> splitIntoTetrahedra(const std::vector<Vertex3d<VectorType, VA>>& poly, 
                                                  const double tol = 0.0);


}

#include "polyclipper3dImpl.hh"

#endif

