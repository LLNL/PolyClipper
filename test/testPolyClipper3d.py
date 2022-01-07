#ATS:test(SELF, label="Polyhedron clipping tests")

import unittest
from math import *
import time

from PolyClipper import *
from PolyClipperTestUtilities import *

# Create a global random number generator.
import random
rangen = random.Random()

#-------------------------------------------------------------------------------
# Make a cube                   |y     
#                               |      
#                               |____x 
#   3/-----/2                  /       
#   /     /|                  /z       
# 7|-----|6|
#  |     | |
#  |  0  | /1
#  |_____|/
#  4     5
#
cube_points = [Vector3d(*coords) for coords in 
               [(0,0,0),  (10,0,0),  (10,10,0),  (0,10,0),
                (0,0,10), (10,0,10), (10,10,10), (0,10,10)]]

cube_neighbors = [[1, 4, 3],
                  [5, 0, 2],
                  [3, 6, 1],
                  [7, 2, 0],
                  [5, 7, 0],
                  [1, 6, 4],
                  [5, 2, 7],
                  [4, 6, 3]]

cube_facets = [[4, 5, 6, 7],
               [1, 2, 6, 5],
               [0, 3, 2, 1],
               [4, 7, 3, 0],
               [6, 2, 3, 7],
               [1, 5, 4, 0]]

#-------------------------------------------------------------------------------
# Make a diamond (2 pyramids (tets) stuck together on the base at z=0.

diamond_points = [Vector3d(*coords) for coords in
                  [(0,0,-1), (0,1,0), (-1,-1,0), (1,-1,0), (0,0,1)]]

diamond_neighbors = [[1, 3, 2],
                     [0, 2, 4, 3],
                     [0, 3, 4, 1],
                     [0, 1, 4, 2],
                     [1, 2, 3]]

diamond_facets = [[0, 2, 1],
                  [0, 1, 3],
                  [0, 3, 2],
                  [1, 4, 3],
                  [2, 4, 1],
                  [2, 3, 4]]

#-------------------------------------------------------------------------------
# Make a non-convex notched thingy.                            |y     
#                                                              |      
#                                                              |____x 
#       6            5       3          2                     /       
#       /------------/      /----------/                     /z       
#      /            / \    /          /|
#     /            /   \ 4/          / |
#    /            /     \/          /  |
#    |------------\      /---------|9  |
#    |13        12 \    / 10       |   |
#    |              \  /           |   |
#    |               \/            |   |
#    |               11            |   |
#    |   0                         |  / 1
#    |                             | /
#    |------------------------------/
#    7                             8

notched_points = [Vector3d(*coords)
                  for coords in [(0,0,0), (4,0,0), (4,2,0), (3,2,0), (2,1,0), (1,2,0), (0,2,0),
                                 (0,0,1), (4,0,1), (4,2,1), (3,2,1), (2,1,1), (1,2,1), (0,2,1)]]
notched_neighbors = [[7, 6, 1],   # 0
                     [0, 2, 8],   # 1
                     [1, 3, 9],   # 2
                     [4, 10, 2],  # 3
                     [5, 11, 3],  # 4
                     [6, 12, 4],  # 5
                     [13, 5, 0],  # 6
                     [8, 13, 0],  # 7
                     [1, 9, 7],   # 8
                     [2, 10, 8],  # 9
                     [9, 3, 11],  # 10
                     [10, 4, 12], # 11
                     [11, 5, 13], # 12
                     [7, 12, 6]]  # 13
notched_facets = [[6, 5, 4, 3, 2, 1, 0],
                  [7, 8, 9, 10, 11, 12, 13],
                  [1, 2, 9, 8],
                  [2, 3, 10, 9],
                  [3, 4, 11, 10],
                  [4, 5, 12, 11],
                  [5, 6, 13, 12],
                  [7, 13, 6, 0],
                  [0, 1, 8, 7]]

#-------------------------------------------------------------------------------
# A degenerate pyramid.  Just reuse the cube, but collapse one face.
degenerate_cube_points1 = [Vector3d(*coords)
                           for coords in [(0,0,0), (1,0,0), (1,1,0), (0,1,0),
                                          (0,0,1), (0,0,1), (0,0,1), (0,0,1)]]

# Another one collapsing a different vertex.
degenerate_cube_points2 = [Vector3d(*coords)
                           for coords in [(0,0,0),  (10,0,0),  (10,10,0),  (10,10,0),
                                          (0,0,10), (10,0,10), (10,10,0), (10,10,0)]]

#-------------------------------------------------------------------------------
# Compute the answer for zeroth and first moments
#-------------------------------------------------------------------------------
def moments_answer(poly):
    assert len(poly) >= 4
    origin = poly[0].position
    m0, m1 = 0.0, Vector3d(0, 0, 0)
    facets = extractFaces(poly)
    for facet in facets:
        n = len(facet)
        p0 = poly[facet[0]].position - origin
        for k in range(1,n-1):
            i = facet[k]
            j = facet[(k + 1) % n]
            p1 = poly[i].position - origin
            p2 = poly[j].position - origin
            dV = p0.dot(p1.cross(p2))
            m0 += dV                                    # 6x
            m1 += dV*(p0 + p1 + p2)                     # 24x
    m0 /= 6.0
    m1 /= 24.0*m0
    m1 += origin
    return m0, m1

#-------------------------------------------------------------------------------
# Compute the average outward normal to a vertex
#-------------------------------------------------------------------------------
def vertexNormal(vertex, poly):
    result = Vector3d(0,0,0)
    for j in vertex.neighbors:
        result += poly[j].position - vertex.position
    return -result.unitVector()

#-------------------------------------------------------------------------------
# Find the maximum length across a polyhedron
#-------------------------------------------------------------------------------
def max_chord(poly):
    xmin = min([v.position.x for v in poly])
    xmax = max([v.position.x for v in poly])
    ymin = min([v.position.y for v in poly])
    ymax = max([v.position.y for v in poly])
    zmin = min([v.position.z for v in poly])
    zmax = max([v.position.z for v in poly])
    return (Vector3d(xmax, ymax, zmax) - Vector3d(xmin, ymin, zmin)).magnitude()

#-------------------------------------------------------------------------------
# Test harness
#-------------------------------------------------------------------------------
class TestPolyhedronClipping(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.convexPolyData = [(cube_points, cube_neighbors, cube_facets),
                               (diamond_points, diamond_neighbors, diamond_facets),
                               (degenerate_cube_points1, cube_neighbors, cube_facets),
                               (degenerate_cube_points2, cube_neighbors, cube_facets)]
        self.nonconvexPolyData = [(notched_points, notched_neighbors, notched_facets)]
        self.degeneratePolyData = [(degenerate_cube_points1, cube_neighbors, cube_facets),
                                   (degenerate_cube_points2, cube_neighbors, cube_facets)]
        self.polyData = self.convexPolyData + self.nonconvexPolyData
        self.ntests = 1000
        return

    #---------------------------------------------------------------------------
    # initializePolyhedron
    #---------------------------------------------------------------------------
    def test_initializePolyhedron(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            m0, m1 = moments(poly)
            M0, M1 = moments_answer(poly)
            self.assertTrue(fuzzyEqual(m0, M0),
                            "Volume comparison failure: %g != %g" % (m0, M0))
            self.assertTrue(fuzzyEqual((m1 - M1).magnitude(), 0.0),
                            "Centroid comparison failure: %s != %s" % (m1, M1))

    #---------------------------------------------------------------------------
    # moments (shift -- moments should be invariant to shifting frame)
    #---------------------------------------------------------------------------
    def test_shiftMoments(self):
        for points, neighbors, facets in self.polyData:
            poly0 = Polyhedron()
            initializePolyhedron(poly0, points, neighbors)
            vol0, centroid0 = moments(poly0)

            # Now shift by some random amount
            delta = Vector3d(rangen.uniform(-1.0, 1.0),
                             rangen.uniform(-1.0, 1.0),
                             rangen.uniform(-1.0, 1.0)) * rangen.uniform(0.1, 10.0)
            poly1 = Polyhedron()
            points1 = [p + delta for p in points]
            initializePolyhedron(poly1, points1, neighbors)
            vol1, centroid1 = moments(poly1)

            self.assertTrue(fuzzyEqual(vol0, vol1),
                            "Volume comparison failure: %g != %g" % (vol0, vol1))
            self.assertTrue(fuzzyEqual((centroid1 - (centroid0 + delta)).magnitude(), 0.0),
                            "Centroid comparison failure: %s != %s" % (centroid0, centroid1 + delta))

    #---------------------------------------------------------------------------
    # collapseDegenerates
    #---------------------------------------------------------------------------
    def test_collapseDegenerates(self):
        for points, neighbors, facets in self.degeneratePolyData:
            poly0 = Polyhedron()
            initializePolyhedron(poly0, points, neighbors)
            assert len(poly0) == len(points)
            #writePolyOBJ(poly0, "poly0.obj")
            poly1 = Polyhedron(poly0)
            #writePolyOBJ(poly1, "poly1_initial.obj")
            collapseDegenerates(poly1, 1.0e-10)
            #writePolyOBJ(poly1, "poly1_collapse.obj")
            assert len(poly1) == 5
            vol0, centroid0 = moments(poly0)
            vol1, centroid1 = moments(poly1)
            #writePolyOBJ(poly1, "poly1_collapse.obj")
            self.assertTrue(fuzzyEqual(vol1, vol0),
                            "Volume comparison failure: %g != %g" % (vol1, vol0))
            self.assertTrue(fuzzyEqual((centroid1 - centroid0).magnitude(), 0.0),
                            "Centroid comparison failure: %s != %s" % (centroid1, centroid0))

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polyhedron.
    #---------------------------------------------------------------------------
    def testClipInternalOnePlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            for i in range(self.ntests):
                planes1, planes2 = [], []
                p0 = Vector3d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                phat = Vector3d(rangen.uniform(-1.0, 1.0), 
                                rangen.uniform(-1.0, 1.0), 
                                rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(Plane3d(p0,  phat))
                planes2.append(Plane3d(p0, -phat))
                chunk1 = Polyhedron(poly)
                chunk2 = Polyhedron(poly)
                clipPolyhedron(chunk1, planes1)
                clipPolyhedron(chunk2, planes2)
                v0, c0 = moments(poly)
                v1, c1 = moments(chunk1)
                v2, c2 = moments(chunk2)
                success = fuzzyEqual(v1 + v2, v0)
                if not success:
                    print("Failed on pass ", i)
                    print("Plane: ", p0, phat)
                    print("Poly:\n", list(poly))
                    print("Chunk 1:\n ", list(chunk1))
                    print("Chunk 2:\n ", list(chunk2))
                    print(moments(chunk1))
                    print(moments(chunk2))
                    print("Vol check: %g + %g = %g" % (v1, v2, v0))
                    writePolyOBJ(poly, "poly.obj")
                    writePolyOBJ(chunk1, "chunk1.obj")
                    writePolyOBJ(chunk2, "chunk2.obj")
                self.assertTrue(success, "Plane clipping summing to wrong volumes: %s + %s = %s != %s" % (v1, v2, v1 + v2, v0))
        return

    #---------------------------------------------------------------------------
    # Clip with a single plane along z = 0
    #---------------------------------------------------------------------------
    def testClipZeroPlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            planes1 = [Plane3d(Vector3d(0,0,0), Vector3d(0,0,1))]
            planes2 = [Plane3d(Vector3d(0,0,0), Vector3d(0,0,-1))]
            chunk1 = Polyhedron(poly)
            chunk2 = Polyhedron(poly)
            clipPolyhedron(chunk1, planes1)
            clipPolyhedron(chunk2, planes2)
            v0, c0 = moments(poly)
            v1, c1 = moments(chunk1)
            v2, c2 = moments(chunk2)
            success = fuzzyEqual(v1 + v2, v0)
            if not success:
                print("Poly:\n", list(poly))
                print("Chunk 1:\n ", list(chunk1))
                print("Chunk 2:\n ", list(chunk2))
                print(moments(chunk1))
                print(moments(chunk2))
                print("Vol check: %g + %g = %g" % (v1, v2, v0))
                writePolyOBJ(poly, "poly.obj")
                writePolyOBJ(chunk1, "chunk1.obj")
                writePolyOBJ(chunk2, "chunk2.obj")
            self.assertTrue(success, "Plane clipping summing to wrong volumes: %s + %s = %s != %s" % (v1, v2, v1 + v2, v0))
        return

    #---------------------------------------------------------------------------
    # Clip with the same plane repeatedly.
    #---------------------------------------------------------------------------
    def testRedundantClip(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            for i in range(self.ntests):
                planes1, planes2 = [], []
                p0 = Vector3d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                phat = Vector3d(rangen.uniform(-1.0, 1.0), 
                                rangen.uniform(-1.0, 1.0),
                                rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(Plane3d(p0,  phat))
                planes2.append(Plane3d(p0,  phat))
                planes2.append(Plane3d(p0,  phat))
                chunk1 = Polyhedron(poly)
                chunk2 = Polyhedron(poly)
                clipPolyhedron(chunk1, planes1)
                clipPolyhedron(chunk2, planes2)
                v0, c0 = moments(poly)
                v1, c1 = moments(chunk1)
                v2, c2 = moments(chunk2)
                success = fuzzyEqual(v1, v2)
                if not success:
                    print("Failed on pass ", i)
                    print("Plane: ", p0, phat)
                    print("Poly:\n", list(poly))
                    print("Chunk 1:\n ", list(chunk1))
                    print("Chunk 2:\n ", list(chunk2))
                    print(moments(chunk1))
                    print(moments(chunk2))
                    print("Vol check: %g = %g" % (v1, v2))
                self.assertTrue(success,
                                "Redundant plane clipping wrong volumes: %s != %s" % (v1, v2))
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polyhedron -- null test.
    #---------------------------------------------------------------------------
    def testNullClipOnePlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            v0, c0 = moments(poly)
            for i in range(self.ntests):
                r = rangen.uniform(2.0, 100.0) * max_chord(poly)
                theta = rangen.uniform(0.0, 2.0*pi)
                phi = rangen.uniform(0.0, pi)
                phat = Vector3d(cos(theta)*sin(phi),
                                sin(theta)*sin(phi),
                                cos(phi))
                p0 = c0 + r*phat
                planes = [Plane3d(p0, -phat)]
                chunk = Polyhedron(poly)
                clipPolyhedron(chunk, planes)
                v1, c1 = moments(chunk)
                self.assertTrue(v1 == v0,
                                "Null plane clipping failure: %s != %s" % (v1, v0))
        return

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polyhedron and rejecting the whole thing.
    #---------------------------------------------------------------------------
    def testFullClipOnePlane(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            v0, c0 = moments(poly)
            for i in range(self.ntests):
                r = rangen.uniform(2.0, 100.0) * max_chord(poly)
                theta = rangen.uniform(0.0, 2.0*pi)
                phi = rangen.uniform(0.0, pi)
                phat = Vector3d(cos(theta)*sin(phi),
                                sin(theta)*sin(phi),
                                cos(phi))
                p0 = c0 + r*phat
                planes = [Plane3d(p0, phat)]
                chunk = Polyhedron(poly)
                clipPolyhedron(chunk, planes)
                self.assertTrue(len(chunk) == 0,
                                "Full plane clipping failure: %s" % chunk)
        return

    #---------------------------------------------------------------------------
    # Clip with planes that exactly hit at least one vertex.  Should reject
    # entire polyhedron.
    #---------------------------------------------------------------------------
    def testFullClipOnePlaneDegenerate(self):
        for points, neighbors in [(cube_points, cube_neighbors),
                                  (diamond_points, diamond_neighbors),
                                  (degenerate_cube_points1, cube_neighbors),
                                  (degenerate_cube_points2, cube_neighbors)]:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            for v in poly:
                phat = vertexNormal(v, poly)
                planes = [Plane3d(v.position, phat)]
                chunk = Polyhedron(poly)
                clipPolyhedron(chunk, planes)
                self.assertTrue(len(chunk) == 0,
                                "Full plane clipping failure: %s" % chunk)

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polyhedron.
    #---------------------------------------------------------------------------
    def testClipInternalTwoPlanes(self):
        for points, neighbors, facets in self.polyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            v0, c0 = moments(poly)
            for i in range(self.ntests):
                p0 = Vector3d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                norm1 = Vector3d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0),
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector3d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0),
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = [Plane3d(p0,  norm1),
                           Plane3d(p0,  norm2)]
                planes2 = [Plane3d(p0,  norm1),
                           Plane3d(p0, -norm2)]
                planes3 = [Plane3d(p0, -norm1),
                           Plane3d(p0,  norm2)]
                planes4 = [Plane3d(p0, -norm1),
                           Plane3d(p0, -norm2)]
                chunk1 = Polyhedron(poly)
                chunk2 = Polyhedron(poly)
                chunk3 = Polyhedron(poly)
                chunk4 = Polyhedron(poly)
                clipPolyhedron(chunk1, planes1)
                clipPolyhedron(chunk2, planes2)
                clipPolyhedron(chunk3, planes3)
                clipPolyhedron(chunk4, planes4)
                v1, c1 = moments(chunk1)
                v2, c2 = moments(chunk2)
                v3, c3 = moments(chunk3)
                v4, c4 = moments(chunk4)
                self.assertTrue(fuzzyEqual(v1 + v2 + v3 + v4, v0),
                                "Two plane clipping summing to wrong volumes: %s + %s + %s + %s = %s != %s" % (v1, v2, v3, v4, v1 + v2 + v3 + v4, v0))
        return

    #---------------------------------------------------------------------------
    # Split a (convex) polyhedron into tetrahedra.
    #---------------------------------------------------------------------------
    def testSplitIntoTetrahedra(self):
        for points, neighbors, facets in self.convexPolyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            tets = splitIntoTetrahedra(poly)
            vol0, centroid0 = moments(poly)
            volTets = 0.0
            centroidTets = Vector3d()
            for inds in tets:
                assert len(inds) == 4
                v0 = poly[inds[0]].position
                v1 = poly[inds[1]].position
                v2 = poly[inds[2]].position
                v3 = poly[inds[3]].position
                V = (v1 - v0).dot((v2 - v0).cross(v3 - v0))
                assert V >= 0.0
                volTets += V
                centroidTets += V*(v0 + v1 + v2 + v3)
            volTets /= 6.0
            centroidTets /= 24.0*volTets
            assert abs(volTets - vol0) < 1.0e-20
            assert (centroidTets - centroid0).magnitude() < 1.0e-20

    #---------------------------------------------------------------------------
    # extractFaces
    #---------------------------------------------------------------------------
    def testExtractFaces(self):
        def minLeadingPermutation(xlist):
            imin = xlist.index(min(xlist))
            return xlist[imin:] + xlist[:imin]
        for points, neighbors, answer in self.convexPolyData:
            answer = [minLeadingPermutation(x) for x in answer]
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            faces = extractFaces(poly)
            faces = [minLeadingPermutation(x) for x in faces]
            assert len(faces) == len(answer)
            for face in faces:
                assert face in answer

    #---------------------------------------------------------------------------
    # commonFaceClips (convex)
    #---------------------------------------------------------------------------
    def testCommonFaceClipsConvex(self):
        for points, neighbors,facets in self.convexPolyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            for i in range(self.ntests):
                p0 = Vector3d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                norm1 = Vector3d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0),
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector3d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0),
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = [Plane3d(p0,  norm1, 10),
                           Plane3d(p0,  norm2, 20)]
                chunk1 = Polyhedron(poly)
                clipPolyhedron(chunk1, planes1)
                faces1 = extractFaces(chunk1)
                clips1 = commonFaceClips(chunk1, faces1)
                assert len(clips1) == len(faces1)
                for clip in clips1:
                    assert len(clip) in (0, 1)
                    for iclip in clip:
                        assert iclip in (10, 20)

    #---------------------------------------------------------------------------
    # commonFaceClips (non-convex)
    #---------------------------------------------------------------------------
    def testCommonFaceClipsNonConvex(self):
        for points, neighbors,facets in self.nonconvexPolyData:
            poly = Polyhedron()
            initializePolyhedron(poly, points, neighbors)
            for i in range(self.ntests):
                p0 = Vector3d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                norm1 = Vector3d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0),
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector3d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0),
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = [Plane3d(p0,  norm1, 10),
                           Plane3d(p0,  norm2, 20)]
                chunk1 = Polyhedron(poly)
                clipPolyhedron(chunk1, planes1)
                faces1 = extractFaces(chunk1)
                clips1 = commonFaceClips(chunk1, faces1)
                assert len(clips1) == len(faces1)
                for clip in clips1:
                    assert len(clip) in (0, 1, 2)
                    for iclip in clip:
                        assert iclip in (10, 20)

if __name__ == "__main__":
    unittest.main()
