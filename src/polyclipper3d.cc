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
using std::pair;
using std::make_pair;
using std::min;
using std::max;
using std::abs;
using std::ostream_iterator;
using std::cerr;
using std::endl;

namespace PolyClipper {

namespace {    // anonymous methods

//------------------------------------------------------------------------------
// Compare a plane and point.
//------------------------------------------------------------------------------
inline
int compare(const Plane3d& plane,
            const PolyClipper::Vector3d& point) {
  return sgn(plane.dist + plane.normal.dot(point));
}

//------------------------------------------------------------------------------
// Compare a plane and a box (defined by it's min/max coordinates).
//   -1 ==> box below plane
//    0 ==> plane cuts through box
//    1 ==> box above plane
//------------------------------------------------------------------------------
inline
int compare(const Plane3d& plane,
            const double xmin,
            const double ymin,
            const double zmin,
            const double xmax,
            const double ymax,
            const double zmax) {
  typedef PolyClipper::Vector3d Vector;
  const auto c1 = compare(plane, Vector(xmin, ymin, zmin));
  const auto c2 = compare(plane, Vector(xmax, ymin, zmin));
  const auto c3 = compare(plane, Vector(xmax, ymax, zmin));
  const auto c4 = compare(plane, Vector(xmin, ymax, zmin));
  const auto c5 = compare(plane, Vector(xmin, ymin, zmax));
  const auto c6 = compare(plane, Vector(xmax, ymin, zmax));
  const auto c7 = compare(plane, Vector(xmax, ymax, zmax));
  const auto c8 = compare(plane, Vector(xmin, ymax, zmax));
  const auto cmin = min(c1, min(c2, min(c3, min(c4, min(c5, min(c6, min(c7, c8)))))));
  const auto cmax = max(c1, max(c2, max(c3, max(c4, max(c5, max(c6, max(c7, c8)))))));
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
PolyClipper::Vector3d
segmentPlaneIntersection(const PolyClipper::Vector3d& a,       // line-segment begin
                         const PolyClipper::Vector3d& b,       // line-segment end
                         const Plane3d& plane) {                 // plane
  const auto asgndist = plane.dist + plane.normal.dot(a);
  const auto bsgndist = plane.dist + plane.normal.dot(b);
  assert (asgndist != bsgndist);
  return (a*bsgndist - b*asgndist)/(bsgndist - asgndist);
}

//------------------------------------------------------------------------------
// Find the next neighbor in CCW order in the neighbor set of a vertex.
// This should be the previous vertex from our entry value.
//------------------------------------------------------------------------------
inline
int
nextInFaceLoop(const Vertex3d& v, const int vprev) {
  const auto itr = find(v.neighbors.begin(), v.neighbors.end(), vprev);
  assert (itr != v.neighbors.end());
  if (itr == v.neighbors.begin()) {
    return v.neighbors.back();
  } else {
    return *(itr - 1);
  }
}

//------------------------------------------------------------------------------
// Perform an in-place intersection of two sets, leaving the result in the
// first argument.
//------------------------------------------------------------------------------
inline
void
inPlaceSetIntersection(set<int>& a, const set<int>& b) {
  auto it1 = a.begin();
  auto it2 = b.begin();
  while ((it1 != a.end()) and (it2 != b.end())) {
    if (*it1 < *it2) {
      a.erase(it1++);
    } else if (*it2 < *it1) {
      ++it2;
    } else { // *it1 == *it2
      ++it1;
      ++it2;
    }
  }

  // Anything left in a from here on did not appear in b, so we remove it.
  a.erase(it1, a.end());
}

}              // anonymous methods

//------------------------------------------------------------------------------
// Initialize a polyhedron given the vertex coordinates and connectivity.
//------------------------------------------------------------------------------
void
initializePolyhedron(Polyhedron& poly,
                     const vector<PolyClipper::Vector3d>& positions,
                     const vector<vector<int>>& neighbors) {

  // Pre-conditions
  const auto n = positions.size();
  assert (neighbors.size() == n);
  poly.resize(n);
  for (auto i = 0; i < n; ++i) {
    assert(neighbors[i].size() >= 3);
    poly[i].position = positions[i];
    poly[i].neighbors = neighbors[i];
  }
}

//------------------------------------------------------------------------------
// Return a nicely formatted string representing the polyhedron.
//------------------------------------------------------------------------------
std::string
polyhedron2string(const Polyhedron& poly) {

  std::ostringstream s;
  const auto nverts = poly.size();
  for (auto i = 0; i < nverts; ++i) {
    s << i << " ID=" << poly[i].ID << " comp=" << poly[i].comp << " @ " << poly[i].position
      << " neighbors=[";
    copy(poly[i].neighbors.begin(), poly[i].neighbors.end(), ostream_iterator<int>(s, " "));
    s << "] clips[";
    copy(poly[i].clips.begin(), poly[i].clips.end(), ostream_iterator<int>(s, " "));
    s << "]\n";
  }

  return s.str();
}

// std::string
// polyhedron2string(const Polyhedron& poly) {
//   ostringstream s;
//   s << "[";

//   // Get the vertices in face ordering.
//   const vector<vector<const Vertex3d*>> faces;

//   // Now output the face vertex coordinates.
//   for (const auto& face: faces) {
//     s << "[";
//     for (const auto vptr: face) {
//       s << " " << vptr->position;
//     }
//     s << "]\n ";
//   }
//   s << "]";

//   return s.str();
// }

//------------------------------------------------------------------------------
// Compute the zeroth and first moment of a Polyhedron.
//------------------------------------------------------------------------------
void moments(double& zerothMoment, PolyClipper::Vector3d& firstMoment,
             const Polyhedron& polyhedron) {

  // Useful types.
  typedef PolyClipper::Vector3d Vector;

  // Clear the result for accumulation.
  zerothMoment = 0.0;
  firstMoment = {0.0, 0.0, 0.0};

  if (not polyhedron.empty()) {

    // Pick a central point to work from.
    // auto nactive = 0;
    // Vector3d origin;
    // for (const auto v: polyhedron) {
    //   if (v.comp >= 0.0) {
    //     ++nactive;
    //     origin += v.position;
    //   }
    // }
    // origin /= nactive;
    const auto origin = polyhedron[0].position;

    // Walk the facets
    const auto facets = extractFaces(polyhedron);
    for (const auto& facet: facets) {
      const auto n = facet.size();
      const auto p0 = polyhedron[facet[0]].position - origin;
      for (auto k = 1; k < n - 1; ++k) {
        const auto i = facet[k];
        const auto j = facet[(k + 1) % n];
        const auto p1 = polyhedron[i].position - origin;
        const auto p2 = polyhedron[j].position - origin;
        const auto dV = p0.dot(p1.cross(p2));
        zerothMoment += dV;                               // 6x
        firstMoment += dV*(p0 + p1 + p2);                 // 24x
      }
    }
    zerothMoment /= 6.0;
    firstMoment *= safeInv(24.0*zerothMoment);
    firstMoment += origin;
  }
}

//------------------------------------------------------------------------------
// Clip a polyhedron by planes.
//------------------------------------------------------------------------------
void clipPolyhedron(Polyhedron& polyhedron,
                    const std::vector<Plane3d>& planes) {

  // Pre-declare variables.  Normally I prefer local declaration, but this
  // seems to slightly help performance.
  bool above, below;
  int nverts0, nverts, nneigh, i, j, k, jn, inew, iprev, inext, itmp;
  vector<int>::iterator nitr;

  // cerr << "Initial:\n" << polyhedron2string(polyhedron) << endl;

  // Find the bounding box of the polyhedron.
  auto xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
  auto ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
  auto zmin = std::numeric_limits<double>::max(), zmax = std::numeric_limits<double>::lowest();
  for (auto& v: polyhedron) {
    xmin = std::min(xmin, v.position.x);
    xmax = std::max(xmax, v.position.x);
    ymin = std::min(ymin, v.position.y);
    ymax = std::max(ymax, v.position.y);
    zmin = std::min(zmin, v.position.z);
    zmax = std::max(zmax, v.position.z);
  }

  // Loop over the planes.
  auto kplane = 0;
  const auto nplanes = planes.size();
  while (kplane < nplanes and not polyhedron.empty()) {
    const auto& plane = planes[kplane++];
    // cerr << "Clip plane: " << plane.dist << " " << plane.normal << endl;

    // First check against the bounding box.
    auto boxcomp = compare(plane, xmin, ymin, zmin, xmax, ymax, zmax);
    auto above = boxcomp ==  1;
    auto below = boxcomp == -1;
    assert (not (above and below));

    // Check the current set of vertices against this plane.
    // Also keep track of any vertices that landed exactly in-plane.
    if (not (above or below)) {
      for (auto& v: polyhedron) {
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
      // The polyhedron is entirely below the clip plane, and is therefore entirely removed.
      // No need to check any more clipping planes -- we're done.
      polyhedron.clear();

    } else if (not above) {

      // This plane passes through the polyhedron.
      // Insert any new vertices.
      nverts0 = polyhedron.size();
      for (i = 0; i < nverts0; ++i) {   // Only check vertices before we start adding new ones.
        if (polyhedron[i].comp == 1) {

          // This vertex survives clipping -- check the neighbors for any new vertices we need to insert.
          nneigh = polyhedron[i].neighbors.size();
          assert (nneigh >= 3);
          for (auto j = 0; j < nneigh; ++j) {
            jn = polyhedron[i].neighbors[j];
            assert (jn < nverts0);
            if (polyhedron[jn].comp == -1) {

              // This edge straddles the clip plane, so insert a new vertex.
              inew = polyhedron.size();
              polyhedron.push_back(Vertex3d(segmentPlaneIntersection(polyhedron[i].position,
                                                                     polyhedron[jn].position,
                                                                     plane),
                                            2));         // 2 indicates new vertex
              assert (polyhedron.size() == inew + 1);
              polyhedron[inew].neighbors = vector<int>({jn, i});
              polyhedron[inew].clips.insert(plane.ID);

              // Patch up clip info -- gotta scan for common elements in the neighbors of the clipped guy.
              std::set<int> common_clips;
              std::set_intersection(polyhedron[i].clips.begin(), polyhedron[i].clips.end(),
                                    polyhedron[jn].clips.begin(), polyhedron[jn].clips.end(),
                                    std::inserter(common_clips, common_clips.begin()));
              polyhedron[inew].clips.insert(common_clips.begin(), common_clips.end());

              nitr = find(polyhedron[jn].neighbors.begin(), polyhedron[jn].neighbors.end(), i);
              assert (nitr != polyhedron[jn].neighbors.end());
              *nitr = inew;
              polyhedron[i].neighbors[j] = inew;
              // cerr << " --> Inserting new vertex @ " << polyhedron.back().position << endl;

            }
          }
        } else if (polyhedron[i].comp == 0) {
          // This vertex is exactly in plane, so just add this plane as a clip.
          polyhedron[i].clips.insert(plane.ID);
        }
      }
      nverts = polyhedron.size();
      // cerr << "After insertion:\n" << polyhedron2string(polyhedron) << endl;

      // // Next handle reconnecting any vertices that were exactly in-plane.
      // for (auto vptr: planeVertices) {
      //   assert (vptr->comp == 0);
      //   const auto nneigh = vptr->neighbors.size();
      //   assert (nneigh >= 3);
      //   for (auto j = 0; j < nneigh; ++j) {
      //     auto nptr = vptr->neighbors[j];
      //     if (nptr->comp == -1) {

      //       // Yep, this edge is clipped, so look for where the in-plane vertex should hook up.
      //       auto vprev = vptr;
      //       auto vnext = nptr;
      //       auto tmp = vnext;
      //       auto k = 0;
      //       while (vnext->comp == -1 and k++ < nverts) {
      //         tmp = vnext;
      //         vnext = nextInFaceLoop(vnext, vprev);
      //         vprev = tmp;
      //       }
      //       assert (vprev->comp == -1);
      //       assert (vnext->comp != -1);
      //       assert (vnext != vptr);
      //       vptr->neighbors[j] = vnext;

      //       const auto barf = ((vnext->position - Vector( 3.46945e-18, 0.0657523, -0.0728318 )).magnitude() < 1.0e-3);
      //       if (barf) cerr << "Deg: " << vptr->position << " " << nneigh << endl;

      //       // Figure out which pointer on the new neighbor should point back at vptr.
      //       auto itr = find(vnext->neighbors.begin(), vnext->neighbors.end(), vprev);
      //       assert (itr != vnext->neighbors.end());
      //       *itr = vptr;
      //     }
      //   }
      // }

      // For each new vertex, link to the neighbors that survive the clipping.
      {
        for (i = nverts0; i < nverts; ++i) {
          assert (polyhedron[i].comp == 2);
          nneigh = polyhedron[i].neighbors.size();

          // Look for any neighbors of the vertex that are clipped.
          for (j = 0; j < nneigh; ++j) {
            jn = polyhedron[i].neighbors[j];
            if (polyhedron[jn].comp == -1) {

              // This neighbor is clipped, so look for the first unclipped vertex along this face loop.
              iprev = i;
              inext = jn;
              itmp = inext;
              // cerr << vitr->ID << ": ( " << vprev->ID << " " << vnext->ID << ")";
              k = 0;
              while (polyhedron[inext].comp == -1 and k++ < nverts) {
                itmp = inext;
                inext = nextInFaceLoop(polyhedron[inext], iprev);
                iprev = itmp;
                // cerr << " (" << vprev->ID << " " << vnext->ID << ")";
              }
              // cerr << endl;
              assert (polyhedron[inext].comp != -1);
              polyhedron[i].neighbors[j] = inext;
              polyhedron[inext].neighbors.insert(polyhedron[inext].neighbors.begin(), i);
              // newEdges.push_back(make_pair(vnext, vitr));
            }
          }
        }
      }
      // const auto nNewEdges = newEdges.size();
      // assert (nNewEdges >= 3);

      // // Collapse any degenerate vertices back onto the originals.
      // {
      //   Vertex3d *vdeg, *vptr;
      //   for (auto& vpair: degenerateVertices) {
      //     tie(vdeg, vptr) = vpair;
      //     vptr->neighbors.reserve(vptr->neighbors.size() + vdeg->neighbors.size());
      //     vptr->neighbors.insert(vptr->neighbors.end(), vdeg->neighbors.begin(), vdeg->neighbors.end());
      //     vdeg->comp = -1;
      //     for (auto& v: polyhedron) { // Get rid of this!
      //       replace(v.neighbors.begin(), v.neighbors.end(), vdeg, vptr);
      //     }
      //   }
      // }

      // Remove the clipped vertices and collapse degenerates, compressing the polyhedron.
      i = 0;
      xmin = std::numeric_limits<double>::max(), xmax = std::numeric_limits<double>::lowest();
      ymin = std::numeric_limits<double>::max(), ymax = std::numeric_limits<double>::lowest();
      zmin = std::numeric_limits<double>::max(), zmax = std::numeric_limits<double>::lowest();
      for (auto& v: polyhedron) {
        if (v.comp >= 0) {
          v.ID = i++;
          xmin = std::min(xmin, v.position.x);
          xmax = std::max(xmax, v.position.x);
          ymin = std::min(ymin, v.position.y);
          ymax = std::max(ymax, v.position.y);
          zmin = std::min(zmin, v.position.z);
          zmax = std::max(zmax, v.position.z);
        }
      }

      // Renumber the neighbor links.
      for (i = 0; i < nverts; ++i) {
        if (polyhedron[i].comp >= 0) {
          assert (polyhedron[i].neighbors.size() >= 3);
          for (j = 0; j < polyhedron[i].neighbors.size(); ++j) {
            polyhedron[i].neighbors[j] = polyhedron[polyhedron[i].neighbors[j]].ID;
          }
        }
      }
      polyhedron.erase(std::remove_if(polyhedron.begin(), polyhedron.end(), [](Vertex3d& v) { return v.comp < 0; }), polyhedron.end());
      // cerr << "After compression:\n" << polyhedron2string(polyhedron) << endl;

      // Is the polyhedron gone?
      if (polyhedron.size() < 4) polyhedron.clear();
    }
  }
}

//------------------------------------------------------------------------------
// Collapse degenerate vertices.
//------------------------------------------------------------------------------
void collapseDegenerates(Polyhedron& polyhedron,
                         const double tol) {

  const auto tol2 = tol*tol;
  auto n = polyhedron.size();
  if (n > 0) {

    // Set the initial ID's the vertices.
    for (auto i = 0; i < n; ++i) polyhedron[i].ID = i;

    // cerr << "Initial: " << endl << polyhedron2string(polyhedron) << endl;

    // Walk the polyhedron removing degenerate edges until we make a sweep without
    // removing any.  Don't worry about ordering of the neighbors yet.
    auto active = false;
    for (auto i = 0; i < n; ++i) {
      if (polyhedron[i].ID >= 0) {
        auto idone = false;
        while (not idone) {
          idone = true;
          for (auto jneigh = 0; jneigh < polyhedron[i].neighbors.size(); ++jneigh) {
            const auto j = polyhedron[i].neighbors[jneigh];
            assert (polyhedron[j].ID >= 0);
            if ((polyhedron[i].position - polyhedron[j].position).magnitude2() < tol2) {
              // cerr << " --> collapasing " << j << " to " << i;
              active = true;
              idone = false;
              polyhedron[j].ID = -1;
              polyhedron[i].clips.insert(polyhedron[j].clips.begin(), polyhedron[j].clips.end());

              // Merge the neighbors of j->i.
              auto jitr = polyhedron[i].neighbors.begin() + jneigh;
              auto kitr = find(polyhedron[j].neighbors.begin(), polyhedron[j].neighbors.end(), i);
              assert (kitr != polyhedron[j].neighbors.end());
              jitr = polyhedron[i].neighbors.insert(jitr, polyhedron[j].neighbors.begin(), kitr);
              jitr = polyhedron[i].neighbors.insert(jitr, kitr+1, polyhedron[j].neighbors.end());  // jitr now points at the first of the newly inserted neighbors

              // Make sure i & j are removed from the neighbor set of i.
              polyhedron[i].neighbors.erase(remove_if(polyhedron[i].neighbors.begin(), polyhedron[i].neighbors.end(),
                                                      [&](const int val) { return val == i or val == j; }), 
                                            polyhedron[i].neighbors.end());

              // Remove any adjacent repeats.
              for (auto kitr = polyhedron[i].neighbors.begin(); kitr < polyhedron[i].neighbors.end() - 1; ++kitr) {
                if (*kitr == *(kitr + 1)) kitr = polyhedron[i].neighbors.erase(kitr);
              }
              if (polyhedron[i].neighbors.front() == polyhedron[i].neighbors.back()) polyhedron[i].neighbors.pop_back();

              // {
              //   cerr << " : new neighbors of " << i << " [";
              //   std::copy(polyhedron[i].neighbors.begin(), polyhedron[i].neighbors.end(), ostream_iterator<int>(cerr, " "));
              //   cerr << "]" << endl;
              // }

              // Make all the neighbors of j point back at i instead of j.
              for (auto k: polyhedron[j].neighbors) {
                 // i is a neighbor to j, and j has already been removed from list
                 // also, i can not be a neighbor to itself
                 if (k != i) {
                    auto itr = find(polyhedron[k].neighbors.begin(), polyhedron[k].neighbors.end(), j);
                    // assert (itr != polyhedron[k].neighbors.end());
                    if (itr != polyhedron[k].neighbors.end()) *itr = i;
                 }
              }
            }
          }
        }
      }
    }

    // cerr << "After relinking: " << endl << polyhedron2string(polyhedron) << endl;

    if (active) {

      // Renumber the nodes assuming we're going to clear out the degenerates.
      auto offset = 0;
      for (auto i = 0; i < n; ++i) {
        if (polyhedron[i].ID == -1) {
          --offset;
        } else {
          polyhedron[i].ID += offset;
        }
      }
      for (auto& v: polyhedron) {
        for (auto itr = v.neighbors.begin(); itr < v.neighbors.end(); ++itr) {
          *itr = polyhedron[*itr].ID;
        }
        v.neighbors.erase(remove_if(v.neighbors.begin(), v.neighbors.end(), [](const int x) { return x < 0; }), v.neighbors.end());

        // Remove any adjacent repeats.
        for (auto kitr = v.neighbors.begin(); kitr < v.neighbors.end() - 1; ++kitr) {
          if (*kitr == *(kitr + 1)) kitr = v.neighbors.erase(kitr);
        }
        if (v.neighbors.front() == v.neighbors.back()) v.neighbors.pop_back();
      }

      // Erase the inactive vertices.
      polyhedron.erase(remove_if(polyhedron.begin(), polyhedron.end(), [](const Vertex3d& v) { return v.ID < 0; }), polyhedron.end());
      if (polyhedron.size() < 4) polyhedron.clear();
    }
  }

  // cerr << "Final: " << endl << polyhedron2string(polyhedron) << endl;

  // Post-conditions.
#ifndef NDEBUG
  {
    const auto n = polyhedron.size();
    for (auto i = 0; i < n; ++i) {
      assert (polyhedron[i].ID == i);
      for (auto j: polyhedron[i].neighbors) assert (j >= 0 and j < n);
    }
  }
#endif

}

//------------------------------------------------------------------------------
// Return the vertices ordered in faces.
// Implicitly uses the convention that neighbors for each vertex are arranged
// counter-clockwise viewed from the exterior.
//------------------------------------------------------------------------------
vector<vector<int>>
extractFaces(const Polyhedron& poly) {

  typedef pair<int, int> Edge;
  typedef vector<int> Face;

  // Prepare the result.
  vector<vector<int>> faceVertices;

  // {
  //   int i = 0;
  //   for (const auto& v: poly) v.ID = i++;
  //   for (const auto& v: poly) {
  //     cerr << " **> " << v.ID << " " << v.position << " :";
  //     for (const auto j: v.neighbors) cerr << " " << poly[j].ID;
  //     cerr << endl;
  //   }
  // }

  // Walk each vertex in the polyhedron.
  set<Edge> edgesWalked;
  const auto nverts = poly.size();
  for (auto i = 0; i < nverts; ++i) {
    const auto& v = poly[i];
    if (v.comp >= 0) {

      // {
      //   cerr << " --> " << v.ID << " " << v.position << " :";
      //   for (const auto j: v.neighbors) cerr << " " << poly[j].ID;
      //   cerr << endl;
      // }

      // Check every (outgoing) edge attached to this vertex.
      for (const auto ni: v.neighbors) {
        assert (poly[ni].comp >= 0);

        // Has this edge been walked yet?
        if (edgesWalked.find(make_pair(ni, i)) == edgesWalked.end()) {
          Face face(1, ni);
          auto vstart = ni;
          auto vnext = i;
          auto vprev = ni;
          // cerr << "Face: " << ni;

          // Follow around the face represented by this edge until we get back
          // to our starting vertex.
          while (vnext != vstart) {
            // cerr << " " << vnext;
            face.push_back(vnext);
            assert (edgesWalked.find(make_pair(vprev, vnext)) == edgesWalked.end());
            edgesWalked.insert(make_pair(vprev, vnext));
            auto itr = find(poly[vnext].neighbors.begin(), poly[vnext].neighbors.end(), vprev);
            assert (itr != poly[vnext].neighbors.end());
            vprev = vnext;
            if (itr == poly[vnext].neighbors.begin()) {
              vnext = poly[vnext].neighbors.back();
            } else {
              vnext = *(itr - 1);
            }
          }
          // cerr << endl;
          edgesWalked.insert(make_pair(vprev, vnext));   // Final edge connecting last->first vertex
          assert (face.size() >= 3);
          faceVertices.push_back(face);
        }
      }
    }
  }

  // Post-conditions.
#ifndef NDEBUG
  {
    // Every pair should have been walked twice, once in each direction.
    for (auto i = 0; i < nverts; ++i) {
      for (const auto ni: poly[i].neighbors) {
        assert (edgesWalked.find(make_pair(i, ni)) != edgesWalked.end());
        assert (edgesWalked.find(make_pair(ni, i)) != edgesWalked.end());
      }
    }
  }
#endif
  return faceVertices;
}

//------------------------------------------------------------------------------
// Compute the set of clips common to each face.
//------------------------------------------------------------------------------
vector<set<int>>
commonFaceClips(const Polyhedron& poly,
                const vector<vector<int>>& faceVertices) {
  const auto nfaces = faceVertices.size();
  vector<set<int>> faceClips(nfaces);
  for (auto k = 0; k < nfaces; ++k) {
    const auto n = faceVertices[k].size();
    assert(n >= 3);
    faceClips[k] = poly[faceVertices[k][0]].clips;
    for (auto j = 1; j < n; ++j) {
      inPlaceSetIntersection(faceClips[k], poly[faceVertices[k][j]].clips);
    }
  }
  return faceClips;
}

//------------------------------------------------------------------------------
// Split a polyhedron into a set of tetrahedra.
//------------------------------------------------------------------------------
vector<vector<int>> splitIntoTetrahedra(const Polyhedron& poly,
                                        const double tol) {

  // Prepare the result, which will be quadruples of indices in the input polyhedron vertices.
  vector<vector<int>> result;

  // Check if we're convex.
  const auto n0 = poly.size();
  bool convex = true;
  auto i = 0;
  while (convex and i < n0) {
    const auto nv = poly[i].neighbors.size();
    assert (nv >= 3);
    auto j = 0;
    while (convex and j < nv - 2) {
      const auto& v0 = poly[i].position;
      const auto& v1 = poly[poly[i].neighbors[j  ]].position;
      const auto& v2 = poly[poly[i].neighbors[j+1]].position;
      const auto& v3 = poly[poly[i].neighbors[j+2]].position;
      convex = ((v1 - v0).dot((v2 - v0).cross(v3 - v0)) <= 1.0e-10);  // Note we have to flip signs cause of neighbor ordering
      ++j;
    }
    ++i;
  }

  // If the polyhedron is convex we can just make a fan of tetrahedra from the first point.
  if (convex) {

    // Get the faces.
    const auto faces = extractFaces(poly);

    // Create tetrahedra from each face back to the starting point, except for faces which contain
    // the starting point since those would be zero volume.
    double vol;
    const auto& v0 = poly[0].position;
    for (const auto& face: faces) {
      // cout << " --> [";
      // copy(face.begin(), face.end(), ostream_iterator<int>(cout, " "));
      // cout << "]" << endl;
      if (find(face.begin(), face.end(), 0) == face.end()) {
        const auto nf = face.size();
        assert (nf >= 3);
        for (auto i = 2; i < nf; ++i) {
          const auto& v1 = poly[face[0  ]].position;
          const auto& v2 = poly[face[i-1]].position;
          const auto& v3 = poly[face[i  ]].position;
          vol = (v1 - v0).dot((v2 - v0).cross(v3 - v0))/3.0;
          if (vol > tol) result.push_back({0, face[0], face[i-1], face[i]});
        }
      }
    }
    return result;
  }

  // PolyClipper::splitIntoTetrahedra ERROR: non-convex polyhedra not supported yet:\n" + polyhedron2string(poly)
  assert(false);
}

}
