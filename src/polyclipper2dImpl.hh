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
#include <list>
#include <map>
#include <set>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <iterator>
#include <algorithm>

using std::vector;
using std::list;
using std::map;
using std::set;
using std::ostream_iterator;
using std::cerr;
using std::endl;

namespace PolyClipper {

namespace internal {

//------------------------------------------------------------------------------
// Compare a plane and a box (defined by it's min/max coordinates).
//   -1 ==> box below plane
//    0 ==> plane cuts through box
//    1 ==> box above plane
//------------------------------------------------------------------------------
template<typename VA>
inline
int compare(const Plane<VA>& plane,
            const double xmin,
            const double ymin,
            const double xmax,
            const double ymax) {
  using Vector = typename VA::VECTOR;
  const auto c1 = internal::compare<VA>(plane, VA::Vector(xmin, ymin));
  const auto c2 = internal::compare<VA>(plane, VA::Vector(xmax, ymin));
  const auto c3 = internal::compare<VA>(plane, VA::Vector(xmax, ymax));
  const auto c4 = internal::compare<VA>(plane, VA::Vector(xmin, ymax));
  const auto cmin = std::min(c1, std::min(c2, std::min(c3, c4)));
  const auto cmax = std::max(c1, std::max(c2, std::max(c3, c4)));
  if (cmin >= 0) {
    return  1;
  } else if (cmax <= 0) {
    return -1;
  } else {
    return  0;
  }
}

//------------------------------------------------------------------------------
// Check if two line segments intersect.
//------------------------------------------------------------------------------
template<typename VA>
inline
bool
segmentsIntersect(const typename VA::VECTOR& a,
                  const typename VA::VECTOR& b,
                  const typename VA::VECTOR& c,
                  const typename VA::VECTOR& d) {
  using PlaneType = Plane<VA>;

  // The plane in the (c,c) orientation.
  PlaneType cdplane;
  cdplane.normal = VA::Vector(-(VA::y(c) - VA::y(d)), VA::x(c) - VA::x(d)).unitVector();
  cdplane.dist = VA::dot(-c, cdplane.normal);

  // Does the (a,b) segment straddle the plane?
  if (internal::compare<VA>(cdplane, a)*internal::compare<VA>(cdplane, b) == 1) return false;

  // Is the point where (a,b) intersects the plane between (c,d)?
  const auto g = internal::segmentPlaneIntersection(a, b, cdplane);
  return VA::dot(VA::sub(c, g), VA::sub(d, g)) <= 0;
}

//------------------------------------------------------------------------------
// Check if a line segment intersects the polygon.
//------------------------------------------------------------------------------
template<typename VA>
inline
bool
intersect(const typename VA::VECTOR& a,            // line-segment begin
          const typename VA::VECTOR& b,            // line-segment end
          const std::vector<Vertex2d<VA>>& poly) { // Polygon
  auto result = false;
  const auto n = poly.size();
  auto i = 0;
  while (i < n and (not result)) {
    result = internal::segmentsIntersect(a, b, poly[i].position, poly[(i+1)%n].position);
  }
  return result;
}

}              // internal namespace methods

//------------------------------------------------------------------------------
// Initialize a polygon given the vertex coordinates and connectivity.
//------------------------------------------------------------------------------
template<typename VA>
void
initializePolygon(std::vector<Vertex2d<VA>>& poly,
                  const vector<typename VA::VECTOR>& positions,
                  const vector<vector<int>>& neighbors) {

  // Pre-conditions
  const auto n = positions.size();
  PCASSERT(neighbors.size() == n);
  poly.resize(n);
  for (auto i = 0; i < n; ++i) {
    PCASSERT(neighbors[i].size() == 2);
    poly[i].position = positions[i];
    poly[i].neighbors = {neighbors[i][0], neighbors[i][1]};
  }
}

//------------------------------------------------------------------------------
// Return a nicely formatted string representing the polygon.
//------------------------------------------------------------------------------
template<typename VA>
std::string
polygon2string(const std::vector<Vertex2d<VA>>& poly) {
  using Vertex = Vertex2d<VA>;

  // Numbers of vertices.
  const auto nverts = poly.size();
  const auto nactive = std::count_if(poly.begin(), poly.end(),
                                     [](const Vertex& x) { return x.comp >= 0; });
  set<int> usedVertices;

  // Dump the raw vertex info.
  std::ostringstream s;
  s << "{\n";
  for (auto i = 0; i < nverts; ++i) {
    s << "  " << i << " " << VA::str(poly[i].position) << " comp=" << poly[i].comp
      << " [" << poly[i].neighbors.first << " " << poly[i].neighbors.second << "]"
      << " clips[";
    std::copy(poly[i].clips.begin(), poly[i].clips.end(), ostream_iterator<int>(s, " "));
    s << "]\n";
  }
  s << "}\n";

  // // Go until we hit all the active vertices.
  // s << "[";
  // while (usedVertices.size() < nactive) {
  //   s << "[";

  //   // Look for the first active unused vertex.
  //   auto vstart = 0;
  //   while (vstart < nverts and
  //          (poly[vstart].comp < 0 or usedVertices.find(vstart) != usedVertices.end())) vstart++;
  //   PCASSERT(vstart < nverts);
  //   auto vnext = vstart;

  //   // Read out this loop.
  //   auto force = true;
  //   while (force or vnext != vstart) {
  //     s << " " << poly[vnext].position;
  //     force = false;
  //     usedVertices.insert(vnext);
  //     vnext = poly[vnext].neighbors.second;
  //   }
  //   s << "]";
  // }
  // s << "]";
  return s.str();
}

//------------------------------------------------------------------------------
// Compute the zeroth and first moment of a Polygon.
//------------------------------------------------------------------------------
template<typename VA>
void moments(double& zerothMoment, typename VA::VECTOR& firstMoment,
             const std::vector<Vertex2d<VA>>& polygon) {

  // Useful types.
  using Vector = typename VA::VECTOR;
  const double nearlyZero = 1.0e-15;

  // Clear the result for accumulation.
  zerothMoment = 0.0;
  firstMoment = VA::Vector(0.0, 0.0);

  // Walk the polygon, and add up our results triangle by triangle.
  if (polygon.size() > 2) {             // Require at least a triangle
    const auto nverts = polygon.size();
    const auto v0 = polygon[0];
    for (const auto v1: polygon) {
      const auto v2 = polygon[v1.neighbors.second];
      const auto triA = VA::crossmag(VA::sub(v1.position, v0.position), VA::sub(v2.position, v0.position));
      zerothMoment += triA;
      VA::iadd(firstMoment, VA::mul(VA::sub(VA::add(v1.position, v2.position), VA::mul(v0.position, 2.0)), triA));
    }
    VA::idiv(firstMoment, 3.0*std::max(nearlyZero, zerothMoment));
    VA::iadd(firstMoment, v0.position);
    zerothMoment *= 0.5;
  }
}

//------------------------------------------------------------------------------
// Clip a polygon by planes.
//------------------------------------------------------------------------------
template<typename VA>
void clipPolygon(std::vector<Vertex2d<VA>>& polygon,
                 const std::vector<Plane<VA>>& planes) {

  // Useful types.
  using Vector = typename VA::VECTOR;
  using Vertex = Vertex2d<VA>;
  const double nearlyZero = 1.0e-15;

  // Prepare to dump the input state if we hit an exception
  std::string initial_state;
#ifndef NDEBUG
  {
    internal::serialize(polygon, initial_state);
    internal::serialize<VA>(planes, initial_state);
  }
#endif

  // Check the input.
  double V0;
  Vector C0;
  moments(V0, C0, polygon);
  if (V0 < nearlyZero) polygon.clear();
  // cerr << "Initial polygon: " << polygon2string(polygon) << " " << V0 << endl;

  // Find the bounding box of the polygon.
  auto xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
  auto ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
  for (auto& v: polygon) {
    xmin = std::min(xmin, VA::x(v.position));
    xmax = std::max(xmax, VA::x(v.position));
    ymin = std::min(ymin, VA::y(v.position));
    ymax = std::max(ymax, VA::y(v.position));
  }

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not polygon.empty()) {
    const auto& plane = planes[kplane++];
    // cerr << "Clip plane: " << plane.dist << " " << VA::str(plane.normal) << endl;

    // First check against the bounding box.
    auto boxcomp = internal::compare(plane, xmin, ymin, xmax, ymax);
    auto above = boxcomp ==  1;
    auto below = boxcomp == -1;
    PCALWAYSASSERT(not (above and below));

    // Check the current set of vertices against this plane.
    if (not (above or below)) {
      above = true;
      below = true;
      for (auto& v: polygon) {
        v.comp = internal::compare<VA>(plane, v.position);
        if (v.comp == 1) {
          below = false;
        } else if (v.comp == -1) {
          above = false;
        }
      }
      PCASSERT2(not (above and below), internal::dumpSerializedState(initial_state));
    }

    // Did we get a simple case?
    if (below) {
      // The polygon is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes -- we're done.
      polygon.clear();

    } else if (not above) {

      // This plane passes through the polygon.
      // Insert any new vertices.
      vector<int> hangingVertices;
      int vprev, vnext, vnew;
      const auto nverts0 = polygon.size();
      for (auto v = 0; v < nverts0; ++v) {
        std::tie(vprev, vnext) = polygon[v].neighbors;

        // BCS XXX note that we are about to attempt to index into polygon.
        // If vnext = -1, for example, this will be a memory error with the sanitizer.
        PCASSERT2(vnext >= 0, internal::dumpSerializedState(initial_state));
        // if (vnext < 0) throw std::logic_error("vnext < 0");
        // if (vprev < 0) throw std::logic_error("vprev < 0");
        // const auto nverts = polygon.size();
        // if (vnext >= nverts) throw std::logic_error("vnext >= nverts");
        // if (vprev >= nverts) throw std::logic_error("vprev >= nverts");

        if ((polygon[v].comp)*(polygon[vnext].comp) == -1) {
          // This pair straddles the plane and creates a new vertex.
          vnew = polygon.size();
          polygon.push_back(Vertex(internal::segmentPlaneIntersection(polygon[v].position,
                                                                      polygon[vnext].position,
                                                                      plane),
                                   2));         // 2 indicates new vertex
          polygon[vnew].neighbors = {v, vnext};
          polygon[vnew].clips.insert(plane.ID);
          // Patch up clip info for existing clips
          if (polygon[v].comp == -1) {
            for (const auto cp: polygon[v].clips) {
              if (polygon[vnext].clips.find(cp) != polygon[vnext].clips.end()) polygon[vnew].clips.insert(cp);
            }
          } else {
            for (const auto cp: polygon[vnext].clips) {
              if (polygon[v].clips.find(cp) != polygon[v].clips.end()) polygon[vnew].clips.insert(cp);
            }
          }
          polygon[v].neighbors.second = vnew;
          polygon[vnext].neighbors.first = vnew;
          hangingVertices.push_back(vnew);
          // cerr << " --> Inserting new vertex @ " << VA::str(polygon.back().position) << endl;

        } else if (polygon[v].comp == 0 and 
                   (polygon[vprev].comp == -1 xor polygon[vnext].comp == -1)) {
          // This vertex is exactly in-plane, but has exactly one neighbor edge that will be entirely clipped.
          // No new vertex, but vptr will be hanging.
          hangingVertices.push_back(v);
          polygon[v].clips.insert(plane.ID);
          // cerr << " --> Hanging vertex @ " << VA::str(polygon[v].position) << endl;

        }
      }
      // cerr << "After insertion: " << polygon2string(polygon) << endl;

      // Look for any topology links to clipped nodes we need to patch.
      const auto nverts = polygon.size();
      size_t i, j, k;
      bool clipped;
      for (i = 0u; i < nverts; ++i) {
        if (polygon[i].comp == 0 or polygon[i].comp == 2) {

          // If my previous neighbor was clipped, walk forward until we encounter
          // another clipped point, and the last unclipped becomes our previous.
          j = polygon[i].neighbors.first;
          if (polygon[j].comp == -1) {
            j = polygon[i].neighbors.second;
            while (polygon[j].comp != -1 and j != i) j = polygon[j].neighbors.second;
            j = polygon[j].neighbors.first;
            PCASSERT2(polygon[j].comp != -1 and j != i,
                      internal::dumpSerializedState(initial_state));
            polygon[i].neighbors.first = j;
            polygon[j].neighbors.second = i;
          }

          // Same thing in reverse for our next neighbor.
          j = polygon[i].neighbors.second;
          if (polygon[j].comp == -1) {
            j = polygon[i].neighbors.first;
            while (polygon[j].comp != -1 and j != i) j = polygon[j].neighbors.first;
            j = polygon[j].neighbors.second;
            PCASSERT2(polygon[j].comp != -1 and j != i,
                      internal::dumpSerializedState(initial_state));
            polygon[i].neighbors.second = j;
            polygon[j].neighbors.first = i;
          }
        }
      }
      // cerr << "After clipping topology: " << polygon2string(polygon) << endl;

      // // For each hanging vertex, link to the neighbors that survive the clipping.
      // // If there are more than two hanging vertices, we've clipped a non-convex face and need to check
      // // how to hook up each section, possibly resulting in new faces.
      // PCASSERT2(hangingVertices.size() % 2 == 0, internal::dumpSerializedState(initial_state));
      // // if (hangingVertices.size() % 2 != 0) throw std::logic_error("hangingVertices mod 2 is not zero");
      // if (true) { //(hangingVertices.size() > 2) {

      //   // Yep, more than one new edge here.
      //   const auto direction = VA::Vector(VA::y(plane.normal), -VA::x(plane.normal));
      //   std::sort(hangingVertices.begin(), hangingVertices.end(), 
      //             [&](const int a, const int b) { return (VA::dot(polygon[a].position, direction) < VA::dot(polygon[b].position, direction)); });

      //   // Now the ordered pairs of these new vertices form the new edges.
      //   int v1, v2;
      //   for (auto k = 0; k < hangingVertices.size(); k += 2) {
      //     v1 = hangingVertices[k];
      //     v2 = hangingVertices[k + 1];
      //     polygon[v1].neighbors.second = v2;
      //     polygon[v2].neighbors.first  = v1;
      //   }

      // } else {

      //   // Just hook across the vertices and we're done.
      //   for (auto v: hangingVertices) {
      //     std::tie(vprev, vnext) = polygon[v].neighbors;
      //     PCASSERT2(polygon[v].comp == 0 or polygon[v].comp == 2, internal::dumpSerializedState(initial_state));
      //     PCASSERT2(polygon[vprev].comp == -1 xor polygon[vnext].comp == -1, internal::dumpSerializedState(initial_state));

      //     if (polygon[vprev].comp == -1) {
      //       // We have to search backwards.
      //       while (polygon[vprev].comp == -1) {
      //         vprev = polygon[vprev].neighbors.first;
      //       }
      //       PCASSERT2(vprev != v, internal::dumpSerializedState(initial_state));
      //       polygon[v].neighbors.first = vprev;

      //     } else {
      //       // We have to search forward.
      //       while (polygon[vnext].comp == -1) {
      //         vnext = polygon[vnext].neighbors.second;
      //       }
      //       PCASSERT2(vnext != v, internal::dumpSerializedState(initial_state));
      //       polygon[v].neighbors.second = vnext;

      //     }
      //   }

      // }

      // Remove the clipped vertices, compressing the polygon.

      // First, number the active vertices sequentially.
      xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
      ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
      size_t nkill = 0u;
      i = 0u;
      for (auto& v: polygon) {
        if (v.comp < 0) {
          nkill++;
        } else {
          v.ID = i++;
          xmin = std::min(xmin, VA::x(v.position));
          xmax = std::max(xmax, VA::x(v.position));
          ymin = std::min(ymin, VA::y(v.position));
          ymax = std::max(ymax, VA::y(v.position));
        }
      }

      // Find the vertices to remove, and renumber the neighbors.
      if (nkill > 0u) {
        vector<int> verts2kill;
        for (auto k = 0; k < polygon.size(); ++k) {
          if (polygon[k].comp < 0) {
            verts2kill.push_back(k);
          } else {
            polygon[k].neighbors.first = polygon[polygon[k].neighbors.first].ID;
            polygon[k].neighbors.second = polygon[polygon[k].neighbors.second].ID;
          }
        }
        internal::removeElements(polygon, verts2kill);
      }

      // cerr << "After compression: " << polygon2string(polygon) << endl;

#ifndef NDEBUG
      if (polygon.size() >= 3) {
        for (const auto& v: polygon) {
          PCASSERT2(v.neighbors.first >= 0 and v.neighbors.second >= 0,
                    v.neighbors.first << " " << v.neighbors.second << "\n" << internal::dumpSerializedState(initial_state));
        }
      }
#endif

      // Is the polygon gone?
      if (polygon.size() < 3) {
        polygon.clear();
      } else {
        double V1;
        Vector C1;
        moments(V1, C1, polygon);
      if (V1 < nearlyZero or
          V1/V0 < 100.0*nearlyZero) polygon.clear();
      }
    }
  }
}

//------------------------------------------------------------------------------
// Collapse degenerate vertices.
//------------------------------------------------------------------------------
template<typename VA>
void collapseDegenerates(std::vector<Vertex2d<VA>>& polygon,
                         const double tol) {

  using VertexType = Vertex2d<VA>;

  const auto tol2 = tol*tol;
  auto n = polygon.size();
  if (n > 0) {

    // Set the initial ID's the vertices.
    for (auto i = 0; i < n; ++i) polygon[i].ID = i;

    // Walk the polygon removing degenerate edges until we make a sweep without
    // removing any.
    auto done = false;
    auto active = false;
    while (not done) {
      done = true;
      for (auto i = 0; i < n; ++i) {
        if (polygon[i].ID >= 0) {
          auto j = polygon[i].neighbors.second;
          PCASSERT(j != i);
          PCASSERT(polygon[j].ID >= 0);
          // if (j==i) throw std::logic_error("got vertex and neighbor identical in collapseDegenerates");
          // if (polygon[j].ID < 0) throw std::logic_error("polygon[j].ID is negative in collapseDegenerates");
          if (VA::magnitude2(VA::sub(polygon[i].position, polygon[j].position)) < tol2) {
            done = false;
            active = true;
            polygon[j].ID = -1;
            int numAttempts = 0;
            while (polygon[j].ID < 0) {
              polygon[i].clips.insert(polygon[j].clips.begin(), polygon[j].clips.end());
              j = polygon[j].neighbors.second;
              ++numAttempts;
              // With the clang-xlf compiler, it's possible to get an infinite loop
              // here if the neighbor patter of the input polygon is bad. It would
              // be better to verify that the neighbor pattern is valid before attempting
              // to collapse degenerates. This is a workaround.
              if (numAttempts > 100) {
                 /*
                 for (int k=0; k<n; ++k) {
                     std::cout << "  orig k = " << k << ", ID = " << copyPoly[k].ID
                       << ", n first = " << copyPoly[k].neighbors.first
                       << ", n second = " << copyPoly[k].neighbors.second
                       << ", pos = ("
                       << std::scientific << std::setprecision(15) << copyPoly[k].position.x
                       << ", "
                       << std::scientific << std::setprecision(15) << copyPoly[k].position.y
                       << ") " << std::endl;
                 }
                 */

                 throw std::logic_error("too many attempts in collapseDegenerates");
              }
            }
            polygon[i].neighbors.second = j;
            polygon[j].neighbors.first  = i;
          }
        }
      }
    }

    if (active) {

      // Renumber the nodes assuming we're going to clear out the degenerates.
      auto offset = 0;
      for (auto i = 0; i < n; ++i) {
        if (polygon[i].ID == -1) {
          --offset;
        } else {
          polygon[i].ID += offset;
        }
      }
      for (auto& v: polygon) {
        v.neighbors.first = polygon[v.neighbors.first].ID;
        v.neighbors.second = polygon[v.neighbors.second].ID;
      }

      // Erase the inactive vertices.
      polygon.erase(remove_if(polygon.begin(), polygon.end(), [](const VertexType& v) { return v.ID < 0; }), polygon.end());
      if (polygon.size() < 3) polygon.clear();
    }
  }

  // Post-conditions.
#ifndef NDEBUG
  {
    const auto n = polygon.size();
    for (auto i = 0; i < n; ++i) {
      PCASSERT(polygon[i].ID == i);
      PCASSERT(polygon[i].neighbors.first < n);
      PCASSERT(polygon[i].neighbors.second < n);
    }
  }
#endif
}

//------------------------------------------------------------------------------
// Return the vertices ordered in faces.
//------------------------------------------------------------------------------
template<typename VA>
vector<vector<int>>
extractFaces(const std::vector<Vertex2d<VA>>& poly) {

  using VertexType = Vertex2d<VA>;
  
  // Numbers of vertices.
  const auto nverts = poly.size();
  const auto nactive = count_if(poly.begin(), poly.end(),
                                [](const VertexType& x) { return x.comp >= 0; });

  // Allocate the result arrays.
  vector<vector<int>> faceVertices(nactive, vector<int>(2));

  // Go until we hit all the active vertices.
  auto k = 0;
  set<int> usedVertices;
  while (usedVertices.size() < nactive) {

    // Look for the first active unused vertex.
    auto vstart = 0;
    while (vstart < nverts and
           (poly[vstart].comp < 0 or usedVertices.find(vstart) != usedVertices.end())) vstart++;
    PCASSERT(vstart < nverts);
    auto vnext = vstart;

    // Read out this loop.
    auto force = true;
    while (force or vnext != vstart) {
      PCASSERT(k < nactive);
      faceVertices[k][0] = vnext;
      vnext = poly[vnext].neighbors.second;
      faceVertices[k][1] = vnext;
      ++k;
      force = false;
      usedVertices.insert(vnext);
    }
    faceVertices[k-1][1] = vstart;
  }

  PCASSERT(k == nactive);
  return faceVertices;
}

//------------------------------------------------------------------------------
// Compute the set of clips common to each face.
//------------------------------------------------------------------------------
template<typename VA>
vector<set<int>>
commonFaceClips(const std::vector<Vertex2d<VA>>& poly,
                const vector<vector<int>>& faceVertices) {

  const auto nfaces = faceVertices.size();
  vector<set<int>> faceClips(nfaces);
  for (auto k = 0; k < nfaces; ++k) {
    PCASSERT(faceVertices[k].size() == 2);
    std::set_intersection(poly[faceVertices[k][0]].clips.begin(), poly[faceVertices[k][0]].clips.end(),
                          poly[faceVertices[k][1]].clips.begin(), poly[faceVertices[k][1]].clips.end(),
                          std::inserter(faceClips[k], faceClips[k].begin()));
  }
  return faceClips;
}

//------------------------------------------------------------------------------
// Split a polygon into a set of triangles.
//------------------------------------------------------------------------------
template<typename VA>
vector<vector<int>> splitIntoTriangles(const std::vector<Vertex2d<VA>>& poly,
                                       const double tol) {

  // Prepare the result, which will be triples of indices in the input polygon vertices.
  vector<vector<int>> result;

  // Check if we're convex.
  const auto n0 = poly.size();
  bool convex = true;
  auto i = 0;
  while (convex and i < n0) {
    convex = VA::crossmag(VA::sub(poly[poly[i].neighbors.second].position, poly[i].position), VA::sub(poly[poly[i].neighbors.first].position, poly[i].position)) >= 0.0;
    ++i;
  }

  // If the polygon is convex we can just make a fan of triangles from the first point.
  if (convex) {
    const auto& v0 = poly[0].position;
    double a;
    for (auto i = 2; i < n0; ++i) {
      const auto& v1 = poly[i-1].position;
      const auto& v2 = poly[i].position;
      a = VA::crossmag(VA::sub(v1, v0), VA::sub(v2, v0));  // really should be 0.5*
      if (a > tol) result.push_back({0, i - 1, i});
    }
    return result;
  }

  // PolyClipper::splitIntoTriangles ERROR: non-convex polygons not supported yet.
  PCASSERT(false);
}

}
