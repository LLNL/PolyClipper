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

#include "polyclipper.hh"
#include "polyclipper_utilities.hh"

#include <list>
#include <map>
#include <set>
#include <iostream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <assert.h>
using std::vector;
using std::list;
using std::map;
using std::set;
using std::ostream_iterator;
using std::cerr;
using std::endl;

namespace PolyClipper {

namespace {    // anonymous methods

//------------------------------------------------------------------------------
// Compare a plane and point.
//------------------------------------------------------------------------------
inline
int compare(const Plane2d& plane,
            const PolyClipper::Vector2d& point) {
  const auto sgndist = plane.dist + plane.normal.dot(point);
  if (std::abs(sgndist) < 1.0e-10) return 0;
  return sgn0(sgndist);
}

//------------------------------------------------------------------------------
// Compare a plane and a box (defined by it's min/max coordinates).
//   -1 ==> box below plane
//    0 ==> plane cuts through box
//    1 ==> box above plane
//------------------------------------------------------------------------------
inline
int compare(const Plane2d& plane,
            const double xmin,
            const double ymin,
            const double xmax,
            const double ymax) {
  typedef PolyClipper::Vector2d Vector;
  const auto c1 = compare(plane, Vector(xmin, ymin));
  const auto c2 = compare(plane, Vector(xmax, ymin));
  const auto c3 = compare(plane, Vector(xmax, ymax));
  const auto c4 = compare(plane, Vector(xmin, ymax));
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
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
inline
PolyClipper::Vector2d
segmentPlaneIntersection(const PolyClipper::Vector2d& a,       // line-segment begin
                         const PolyClipper::Vector2d& b,       // line-segment end
                         const Plane2d& plane) {               // plane
  const auto asgndist = plane.dist + plane.normal.dot(a);
  const auto bsgndist = plane.dist + plane.normal.dot(b);
  assert(asgndist != bsgndist);
  return (a*bsgndist - b*asgndist)/(bsgndist - asgndist);
}

//------------------------------------------------------------------------------
// Check if two line segments intersect.
//------------------------------------------------------------------------------
inline
bool
segmentsIntersect(const PolyClipper::Vector2d& a,
                  const PolyClipper::Vector2d& b,
                  const PolyClipper::Vector2d& c,
                  const PolyClipper::Vector2d& d) {

  // The plane in the (c,c) orientation.
  Plane2d cdplane;
  cdplane.normal = PolyClipper::Vector2d(-(c.y - d.y), c.x - d.x).unitVector();
  cdplane.dist = -c.dot(cdplane.normal);

  // Does the (a,b) segment straddle the plane?
  if (compare(cdplane, a)*compare(cdplane, b) == 1) return false;

  // Is the point where (a,b) intersects the plane between (c,d)?
  const auto g = segmentPlaneIntersection(a, b, cdplane);
  return (c - g).dot(d - g) <= 0;
}

//------------------------------------------------------------------------------
// Check if a line segment intersects the polygon.
//------------------------------------------------------------------------------
inline
bool
intersect(const PolyClipper::Vector2d& a,       // line-segment begin
          const PolyClipper::Vector2d& b,       // line-segment end
          const Polygon& poly) {                  // Polygon
  auto result = false;
  const auto n = poly.size();
  auto i = 0;
  while (i < n and (not result)) {
    result = segmentsIntersect(a, b, poly[i].position, poly[(i+1)%n].position);
  }
  return result;
}

}              // anonymous methods

//------------------------------------------------------------------------------
// Initialize a polygon given the vertex coordinates and connectivity.
//------------------------------------------------------------------------------
void
initializePolygon(Polygon& poly,
                  const vector<PolyClipper::Vector2d>& positions,
                  const vector<vector<int>>& neighbors) {

  // Pre-conditions
  const auto n = positions.size();
  assert (neighbors.size() == n);
  poly.resize(n);
  for (auto i = 0; i < n; ++i) {
    assert (neighbors[i].size() == 2);
    poly[i].position = positions[i];
    poly[i].neighbors = {neighbors[i][0], neighbors[i][1]};
  }
}

//------------------------------------------------------------------------------
// Return a nicely formatted string representing the polygon.
//------------------------------------------------------------------------------
std::string
polygon2string(const Polygon& poly) {
  std::ostringstream s;

  // Numbers of vertices.
  const auto nverts = poly.size();
  const auto nactive = count_if(poly.begin(), poly.end(),
                                [](const Vertex2d& x) { return x.comp >= 0; });
  set<int> usedVertices;

  // Dump the raw vertex info.
  s << "{\n";
  for (auto i = 0; i < nverts; ++i) {
    s << "  " << i << " " << poly[i].position
      << " [" << poly[i].neighbors.first << " " << poly[i].neighbors.second << "]"
      << " clips[";
    copy(poly[i].clips.begin(), poly[i].clips.end(), ostream_iterator<int>(s, " "));
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
  //   assert (vstart < nverts);
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
void moments(double& zerothMoment, PolyClipper::Vector2d& firstMoment,
             const Polygon& polygon) {

  // Useful types.
  typedef PolyClipper::Vector2d Vector;

  // Clear the result for accumulation.
  zerothMoment = 0.0;
  firstMoment = {0.0, 0.0};

  // Walk the polygon, and add up our results triangle by triangle.
  if (not polygon.empty()) {
    const auto nverts = polygon.size();
    const auto v0 = polygon[0];
    for (const auto v1: polygon) {
      const auto v2 = polygon[v1.neighbors.second];
      const auto triA = (v1.position - v0.position).cross(v2.position - v0.position);
      zerothMoment += triA;
      firstMoment += triA * (v1.position + v2.position - 2.0*v0.position);
    }
    assert (zerothMoment != 0.0);
    firstMoment = firstMoment/(3.0*zerothMoment) + v0.position;
    zerothMoment *= 0.5;
  }
}

//------------------------------------------------------------------------------
// Clip a polygon by planes.
//------------------------------------------------------------------------------
void clipPolygon(Polygon& polygon,
                 const std::vector<Plane2d>& planes) {

  // Useful types.
  typedef PolyClipper::Vector2d Vector;

  // cerr << "Initial polygon: " << polygon2string(polygon) << endl;

  // Find the bounding box of the polygon.
  auto xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
  auto ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
  for (auto& v: polygon) {
    xmin = std::min(xmin, v.position.x);
    xmax = std::max(xmax, v.position.x);
    ymin = std::min(ymin, v.position.y);
    ymax = std::max(ymax, v.position.y);
  }

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not polygon.empty()) {
    const auto& plane = planes[kplane++];
    // cerr << "Clip plane: " << plane.dist << " " << plane.normal << endl;

    // First check against the bounding box.
    auto boxcomp = compare(plane, xmin, ymin, xmax, ymax);
    auto above = boxcomp ==  1;
    auto below = boxcomp == -1;
    assert (not (above and below));

    // Check the current set of vertices against this plane.
    if (not (above or below)) {
      for (auto& v: polygon) {
        v.comp = compare(plane, v.position);
        if (v.comp == 1) {
          below = false;
        } else if (v.comp == -1) {
          above = false;
        }
      }
      assert (not (above and below));
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

        if ((polygon[v].comp)*(polygon[vnext].comp) == -1) {
          // This pair straddles the plane and creates a new vertex.
          vnew = polygon.size();
          polygon.push_back(Vertex2d(segmentPlaneIntersection(polygon[v].position,
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
          // cerr << " --> Inserting new vertex @ " << polygon.back().position << endl;

        } else if (polygon[v].comp == 0 and 
                   (polygon[vprev].comp == -1 xor polygon[vnext].comp == -1)) {
          // This vertex is exactly in-plane, but has exactly one neighbor edge that will be entirely clipped.
          // No new vertex, but vptr will be hanging.
          hangingVertices.push_back(v);
          polygon[v].clips.insert(plane.ID);
          // cerr << " --> Hanging vertex @ " << polygon[v].position << endl;

        }
      }

      // cerr << "After insertion: " << polygon2string(polygon) << endl;

      // For each hanging vertex, link to the neighbors that survive the clipping.
      // If there are more than two hanging vertices, we've clipped a non-convex face and need to check
      // how to hook up each section, possibly resulting in new faces.
      assert (hangingVertices.size() % 2 == 0);
      if (true) { //(hangingVertices.size() > 2) {

        // Yep, more than one new edge here.
        const Vector direction(plane.normal.y, -(plane.normal.x));
        sort(hangingVertices.begin(), hangingVertices.end(), 
             [&](const int a, const int b) { return (polygon[a].position).dot(direction) < (polygon[b].position).dot(direction); });

        // Now the ordered pairs of these new vertices form the new edges.
        int v1, v2;
        for (auto k = 0; k < hangingVertices.size(); k += 2) {
          v1 = hangingVertices[k];
          v2 = hangingVertices[k + 1];
          polygon[v1].neighbors.second = v2;
          polygon[v2].neighbors.first  = v1;
        }

      } else {

        // Just hook across the vertices and we're done.
        for (auto v: hangingVertices) {
          std::tie(vprev, vnext) = polygon[v].neighbors;
          assert (polygon[v].comp == 0 or polygon[v].comp == 2);
          assert (polygon[vprev].comp == -1 xor polygon[vnext].comp == -1);

          if (polygon[vprev].comp == -1) {
            // We have to search backwards.
            while (polygon[vprev].comp == -1) {
              vprev = polygon[vprev].neighbors.first;
            }
            assert (vprev != v);
            polygon[v].neighbors.first = vprev;

          } else {
            // We have to search forward.
            while (polygon[vnext].comp == -1) {
              vnext = polygon[vnext].neighbors.second;
            }
            assert (vnext != v);
            polygon[v].neighbors.second = vnext;

          }
        }

      }

      // Remove the clipped vertices, compressing the polygon.

      // First, number the active vertices sequentially.
      xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
      ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
      auto i = 0;
      auto nkill = 0;
      for (auto& v: polygon) {
        if (v.comp < 0) {
          nkill++;
        } else {
          v.ID = i++;
          xmin = std::min(xmin, v.position.x);
          xmax = std::max(xmax, v.position.x);
          ymin = std::min(ymin, v.position.y);
          ymax = std::max(ymax, v.position.y);
        }
      }

      // Find the vertices to remove, and renumber the neighbors.
      if (nkill > 0) {
        vector<int> verts2kill;
        for (auto k = 0; k < polygon.size(); ++k) {
          if (polygon[k].comp < 0) {
            verts2kill.push_back(k);
          } else {
            polygon[k].neighbors.first = polygon[polygon[k].neighbors.first].ID;
            polygon[k].neighbors.second = polygon[polygon[k].neighbors.second].ID;
          }
        }
        removeElements(polygon, verts2kill);
      }

      // cerr << "After compression: " << polygon2string(polygon) << endl;

      // Is the polygon gone?
      if (polygon.size() < 3) polygon.clear();
    }
  }
}

//------------------------------------------------------------------------------
// Collapse degenerate vertices.
//------------------------------------------------------------------------------
void collapseDegenerates(Polygon& polygon,
                         const double tol) {

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
          assert (polygon[j].ID >= 0);
          if ((polygon[i].position - polygon[j].position).magnitude2() < tol2) {
            done = false;
            active = true;
            polygon[j].ID = -1;
            while (polygon[j].ID < 0) {
              polygon[i].clips.insert(polygon[j].clips.begin(), polygon[j].clips.end());
              j = polygon[j].neighbors.second;
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
      polygon.erase(remove_if(polygon.begin(), polygon.end(), [](const Vertex2d& v) { return v.ID < 0; }), polygon.end());
      if (polygon.size() < 3) polygon.clear();
    }
  }

  // Post-conditions.
#ifndef NDEBUG
  {
    const auto n = polygon.size();
    for (auto i = 0; i < n; ++i) {
      assert(polygon[i].ID == i);
      assert(polygon[i].neighbors.first < n);
      assert(polygon[i].neighbors.second < n);
    }
  }
#endif
}

//------------------------------------------------------------------------------
// Return the vertices ordered in faces.
//------------------------------------------------------------------------------
vector<vector<int>>
extractFaces(const Polygon& poly) {

  // Numbers of vertices.
  const auto nverts = poly.size();
  const auto nactive = count_if(poly.begin(), poly.end(),
                                [](const Vertex2d& x) { return x.comp >= 0; });

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
    assert(vstart < nverts);
    auto vnext = vstart;

    // Read out this loop.
    auto force = true;
    while (force or vnext != vstart) {
      assert(k < nactive);
      faceVertices[k][0] = vnext;
      vnext = poly[vnext].neighbors.second;
      faceVertices[k][1] = vnext;
      ++k;
      force = false;
      usedVertices.insert(vnext);
    }
    faceVertices[k-1][1] = vstart;
  }

  assert(k == nactive);
  return faceVertices;
}

//------------------------------------------------------------------------------
// Compute the set of clips common to each face.
//------------------------------------------------------------------------------
vector<set<int>>
commonFaceClips(const Polygon& poly,
                const vector<vector<int>>& faceVertices) {

  const auto nfaces = faceVertices.size();
  vector<set<int>> faceClips(nfaces);
  for (auto k = 0; k < nfaces; ++k) {
    assert(faceVertices[k].size() == 2);
    std::set_intersection(poly[faceVertices[k][0]].clips.begin(), poly[faceVertices[k][0]].clips.end(),
                          poly[faceVertices[k][1]].clips.begin(), poly[faceVertices[k][1]].clips.end(),
                          std::inserter(faceClips[k], faceClips[k].begin()));
  }
  return faceClips;
}

//------------------------------------------------------------------------------
// Split a polygon into a set of triangles.
//------------------------------------------------------------------------------
vector<vector<int>> splitIntoTriangles(const Polygon& poly,
                                       const double tol) {

  // Prepare the result, which will be triples of indices in the input polygon vertices.
  vector<vector<int>> result;

  // Check if we're convex.
  const auto n0 = poly.size();
  bool convex = true;
  auto i = 0;
  while (convex and i < n0) {
    convex = ((poly[poly[i].neighbors.second].position - poly[i].position).cross((poly[poly[i].neighbors.first].position - poly[i].position)) >= 0.0);
    ++i;
  }

  // If the polygon is convex we can just make a fan of triangles from the first point.
  if (convex) {
    const auto& v0 = poly[0].position;
    double a;
    for (auto i = 2; i < n0; ++i) {
      const auto& v1 = poly[i-1].position;
      const auto& v2 = poly[i].position;
      a = 0.5*(v1 - v0).cross(v2 - v0);
      if (a > tol) result.push_back({0, i - 1, i});
    }
    return result;
  }

  // PolyClipper::splitIntoTriangles ERROR: non-convex polygons not supported yet.
  assert (false);
}

}
