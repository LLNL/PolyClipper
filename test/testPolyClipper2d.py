#ATS:test(SELF, label="PolyClipper 2D (polygon) tests")

import unittest
from math import *
import time

from PolyClipper import *
from PolyClipperTestUtilities import *

# Create a global random number generator.
import random
rangen = random.Random()

#-------------------------------------------------------------------------------
# Make a square
#   3    2
#   |----|
#   |    |
#   |----|
#   0    1
square_points = []
for coords in [(0,0), (10,0), (10,10), (0,10)]:
    square_points.append(Vector2d(*coords))

#-------------------------------------------------------------------------------
# Make a non-convex notched thingy.
#    6           5      3          2
#    |------------\    /-----------|
#    |             \  /            |
#    |              \/             |
#    |               4             |
#    |                             |
#    |------------------------------
#    0                             1
notched_points = []
for coords in [(0,0), (4,0), (4,2), (3,2), (2,1), (1,2), (0,2)]:
    notched_points.append(Vector2d(*coords))

#-------------------------------------------------------------------------------
# A degenerate square
#   4    3
#  5|----|
#   |    |
#   |----|2
#   0    1
degenerate_square_points = []
for coords in [(0,0), (1,0), (1,0), (1,1), (0,1), (0,1)]:
    degenerate_square_points.append(Vector2d(*coords))

#-------------------------------------------------------------------------------
# Compute the answer for zeroth and first moments
#-------------------------------------------------------------------------------
def moments_answer(poly):
    n = len(poly)
    assert n >= 3
    v1 = poly[0]
    m0, m1 = 0.0, Vector2d(0, 0)
    for v2 in poly:
        v3 = poly[v2.neighbors[1]]
        Ai = (v2.position - v1.position).cross(v3.position - v1.position)  # 2x
        m0 += Ai
        m1 += (v1.position + v2.position + v3.position)*Ai                 # 6x
    return 0.5*m0, (v1.position + m1/(3.0*m0))

#-------------------------------------------------------------------------------
# Compute the vertex neighbors assuming an ordered ring of the given size.
#-------------------------------------------------------------------------------
def vertexNeighbors(points):
    n = len(points)
    neighbors = []
    for i in xrange(n):
        neighbors.append([(i - 1) % n,
                          (i + 1) % n])
    return neighbors

#-------------------------------------------------------------------------------
# Compute the facets assuming an ordered ring of the given size.
#-------------------------------------------------------------------------------
def facets(points):
    n = len(points)
    facets = []
    for i in xrange(n):
        facets.append([i, (i + 1) % n])
    return facets

#-------------------------------------------------------------------------------
# Find the maximum length across a polygon
#-------------------------------------------------------------------------------
def max_chord(poly):
    xmin = min([v.position.x for v in poly])
    xmax = max([v.position.x for v in poly])
    ymin = min([v.position.y for v in poly])
    ymax = max([v.position.y for v in poly])
    return (Vector2d(xmax, ymax) - Vector2d(xmin, ymin)).magnitude()

#-------------------------------------------------------------------------------
# Test harness
#-------------------------------------------------------------------------------
class TestPolyClipper2d(unittest.TestCase):

    #---------------------------------------------------------------------------
    # setUp
    #---------------------------------------------------------------------------
    def setUp(self):
        self.convexPointSets = [square_points, degenerate_square_points]
        self.nonconvexPointSets = [notched_points]
        self.pointSets = self.convexPointSets + self.nonconvexPointSets
        self.ntests = 5 # 10000

    #---------------------------------------------------------------------------
    # initializePolygon  (Also tests moments)
    #---------------------------------------------------------------------------
    def test_initializePolygon(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            vol, centroid = moments(poly)
            vol0, centroid0 = moments_answer(poly)
            self.failUnless(vol == vol0,
                            "Volume comparison failure: %g != %g" % (vol, vol0))
            self.failUnless(centroid == centroid0,
                            "Centroid comparison failure: %s != %s" % (centroid, centroid0))

    #---------------------------------------------------------------------------
    # collapseDegenerates
    #---------------------------------------------------------------------------
    def test_collapseDegenerates(self):
        poly0 = Polygon()
        initializePolygon(poly0, degenerate_square_points, vertexNeighbors(degenerate_square_points))
        assert len(poly0) == len(degenerate_square_points)
        poly1 = Polygon(poly0)
        collapseDegenerates(poly1, 1.0e-10)
        assert len(poly1) == 4
        vol0, centroid0 = moments(poly0)
        vol1, centroid1 = moments(poly1)
        assert vol1 == vol0
        assert centroid1 == centroid0

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polygon.
    #---------------------------------------------------------------------------
    def testClipInternalOnePlane(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            for i in xrange(self.ntests):
                planes1, planes2 = [], []
                p0 = Vector2d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                phat = Vector2d(rangen.uniform(-1.0, 1.0), 
                                rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(Plane2d(p0,  phat))
                planes2.append(Plane2d(p0, -phat))
                chunk1 = Polygon(poly)
                chunk2 = Polygon(poly)
                clipPolygon(chunk1, planes1)
                clipPolygon(chunk2, planes2)
                v0, c0 = moments(poly)
                v1, c1 = moments(chunk1)
                v2, c2 = moments(chunk2)
                success = fuzzyEqual(v1 + v2, v0)
                if not success:
                    print "Failed on pass ", i
                    print "Plane: ", p0, phat
                    print "Poly:\n", list(poly)
                    print "Chunk 1:\n ", list(chunk1)
                    print "Chunk 2:\n ", list(chunk2)
                    print moments(chunk1)
                    print moments(chunk2)
                    print "Vol check: %g + %g = %g" % (v1, v2, v0)
                self.failUnless(success, "Plane clipping summing to wrong volumes: %s + %s = %s != %s" % (v1, v2, v1 + v2, v0))

    #---------------------------------------------------------------------------
    # Clip with the same plane repeatedly.
    #---------------------------------------------------------------------------
    def testRedundantClip(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            for i in xrange(self.ntests):
                planes1, planes2 = [], []
                p0 = Vector2d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                phat = Vector2d(rangen.uniform(-1.0, 1.0), 
                                rangen.uniform(-1.0, 1.0)).unitVector()
                planes1.append(Plane2d(p0,  phat))
                planes2.append(Plane2d(p0,  phat))
                planes2.append(Plane2d(p0,  phat))
                chunk1 = Polygon(poly)
                chunk2 = Polygon(poly)
                clipPolygon(chunk1, planes1)
                clipPolygon(chunk2, planes2)
                v0, c0 = moments(poly)
                v1, c1 = moments(chunk1)
                v2, c2 = moments(chunk2)
                success = fuzzyEqual(v1, v2)
                if not success:
                    print "Failed on pass ", i
                    print "Plane: ", p0, phat
                    print "Poly:\n", list(poly)
                    print "Chunk 1:\n ", list(chunk1)
                    print "Chunk 2:\n ", list(chunk2)
                    print moments(chunk1)
                    print moments(chunk2)
                    print "Vol check: %g = %g" % (v1, v2)
                self.failUnless(success,
                                "Redundant plane clipping wrong volumes: %s != %s" % (v1, v2))

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polygon -- null test.
    #---------------------------------------------------------------------------
    def testNullClipOnePlane(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            v0, c0 = moments(poly)
            for i in xrange(self.ntests):
                r = rangen.uniform(2.0, 100.0) * max_chord(poly)
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector2d(cos(theta), sin(theta))
                p0 = c0 + r*phat
                planes = [Plane2d(p0, -phat)]
                chunk = Polygon(poly)
                clipPolygon(chunk, planes)
                v1, c1 = moments(chunk)
                self.failUnless(v1 == v0,
                                "Null plane clipping failure: %s != %s" % (v1, v0))

    #---------------------------------------------------------------------------
    # Clip with planes passing outside the polygon and rejecting the whole thing.
    #---------------------------------------------------------------------------
    def testFullClipOnePlane(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            v0, c0 = moments(poly)
            for i in xrange(self.ntests):
                r = rangen.uniform(2.0, 100.0) * max_chord(poly)
                theta = rangen.uniform(0.0, 2.0*pi)
                phat = Vector2d(cos(theta), sin(theta))
                p0 = c0 + r*phat
                planes = [Plane2d(p0, phat)]
                chunk = Polygon(poly)
                clipPolygon(chunk, planes)
                v1, c1 = moments(chunk)
                self.failUnless(v1 == 0.0,
                                "Full plane clipping failure: %s != %s" % (v1, 0.0))

    #---------------------------------------------------------------------------
    # Clip with planes passing through the polygon.
    #---------------------------------------------------------------------------
    def testClipInternalTwoPlanes(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            v0, c0 = moments(poly)
            for i in xrange(self.ntests):
                p0 = Vector2d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                norm1 = Vector2d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector2d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = [Plane2d(p0,  norm1),
                           Plane2d(p0,  norm2)]
                planes2 = [Plane2d(p0,  norm1),
                           Plane2d(p0, -norm2)]
                planes3 = [Plane2d(p0, -norm1),
                           Plane2d(p0,  norm2)]
                planes4 = [Plane2d(p0, -norm1),
                           Plane2d(p0, -norm2)]
                chunk1 = Polygon(poly)
                chunk2 = Polygon(poly)
                chunk3 = Polygon(poly)
                chunk4 = Polygon(poly)
                clipPolygon(chunk1, planes1)
                clipPolygon(chunk2, planes2)
                clipPolygon(chunk3, planes3)
                clipPolygon(chunk4, planes4)
                v1, c1 = moments(chunk1)
                v2, c2 = moments(chunk2)
                v3, c3 = moments(chunk3)
                v4, c4 = moments(chunk4)
                self.failUnless(fuzzyEqual(v1 + v2 + v3 + v4, v0),
                                "Two plane clipping summing to wrong volumes: %s + %s + %s + %s = %s != %s" % (v1, v2, v3, v4, v1 + v2 + v3 + v4, v0))

    # #---------------------------------------------------------------------------
    # # Split a (convex) polygon into triangles.
    # #---------------------------------------------------------------------------
    # def testSplitIntoTriangles(self):
    #     for points in self.convexPointSets:
    #         PCpoly = Polygon()
    #         initializePolygon(PCpoly, points, vertexNeighbors(points))
    #         tris = splitIntoTriangles(PCpoly)
    #         vol0, centroid0 = moments(PCpoly)
    #         volTris = 0.0
    #         centroidTris = Vector2d()
    #         for inds in tris:
    #             assert len(inds) == 3
    #             a = ((PCpoly[inds[1]].position - PCpoly[inds[0]].position).cross(PCpoly[inds[2]].position - PCpoly[inds[0]].position).z)
    #             assert a >= 0.0
    #             volTris += a
    #             centroidTris += a*(PCpoly[inds[0]].position + PCpoly[inds[1]].position + PCpoly[inds[2]].position)
    #         volTris *= 0.5
    #         centroidTris /= 6.0*volTris
    #         assert abs(volTris - vol0) < 1.0e-20
    #         assert (centroidTris - centroid0).magnitude() < 1.0e-20

if __name__ == "__main__":
    unittest.main()
