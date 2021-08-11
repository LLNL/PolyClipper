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
square_points = [Vector2d(*coords) 
                 for coords in [(0,0), (10,0), (10,10), (0,10)]]

#-------------------------------------------------------------------------------
# Diamond
#      2
#     / \
#    /   \
#   3     1
#    \   /
#     \ /
#      0
diamond_points = [Vector2d(*coords)
                  for coords in [(0,-1), (1,0), (0,1), (-1,0)]]

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
notched_points = [Vector2d(*coords)
                  for coords in [(0,0), (4,0), (4,2), (3,2), (2,1), (1,2), (0,2)]]

#-------------------------------------------------------------------------------
# A degenerate square
#   4    3
#  5|----|
#   |    |
#   |----|2
#   0    1
degenerate_square_points = [Vector2d(*coords)
                            for coords in [(0,0), (1,0), (1,0), (1,1), (0,1), (0,1)]]

#-------------------------------------------------------------------------------
# Compute the answer for zeroth and first moments
#-------------------------------------------------------------------------------
def moments_answer(poly):
    n = len(poly)
    assert n >= 3
    v1 = poly[0]
    m0, m1 = 0.0, Vector2d(0, 0)
    for v2 in poly[1:-1]:
        v3 = poly[v2.neighbors[1]]
        Ai = (v2.position - v1.position).crossmag(v3.position - v1.position)  # 2x
        m0 += Ai
        m1 += (v1.position + v2.position + v3.position)*Ai                    # 6x
    return 0.5*m0, m1/(3.0*m0)

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
# Compute the average outward normal to a vertex
#-------------------------------------------------------------------------------
def vertexNormal(vertex, poly):
    return -(poly[vertex.neighbors[0]].position - vertex.position +
             poly[vertex.neighbors[1]].position - vertex.position).unitVector()

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
        self.convexPointSets = [square_points, diamond_points, degenerate_square_points]
        self.nonconvexPointSets = [notched_points]
        self.pointSets = self.convexPointSets + self.nonconvexPointSets
        self.ntests = 1000

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
    # moments (shift -- moments should be invariant to shifting frame)
    #---------------------------------------------------------------------------
    def test_shiftMoments(self):
        for points in self.pointSets:
            poly0 = Polygon()
            initializePolygon(poly0, points, vertexNeighbors(points))
            vol0, centroid0 = moments(poly0)

            # Now shift by some random amount
            delta = Vector2d(rangen.uniform(-1.0, 1.0),
                             rangen.uniform(-1.0, 1.0)) * rangen.uniform(0.1, 10.0)
            poly1 = Polygon()
            points1 = [p + delta for p in points]
            initializePolygon(poly1, points1, vertexNeighbors(points1))
            vol1, centroid1 = moments(poly1)

            self.failUnless(fuzzyEqual(vol0, vol1),
                            "Volume comparison failure: %g != %g" % (vol0, vol1))
            self.failUnless(fuzzyEqual((centroid1 - (centroid0 + delta)).magnitude(), 0.0),
                            "Centroid comparison failure: %s != %s" % (centroid0, centroid1 + delta))

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
    # Clip with a single plane along y = 0
    #---------------------------------------------------------------------------
    def testClipZeroPlane(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            planes1 = [Plane2d(Vector2d(0,0), Vector2d(0,1))]
            planes2 = [Plane2d(Vector2d(0,0), Vector2d(0,-1))]
            chunk1 = Polygon(poly)
            chunk2 = Polygon(poly)
            clipPolygon(chunk1, planes1)
            clipPolygon(chunk2, planes2)
            v0, c0 = moments(poly)
            v1, c1 = moments(chunk1)
            v2, c2 = moments(chunk2)
            success = fuzzyEqual(v1 + v2, v0)
            if not success:
                print "Plane: ", p0, phat
                print "Poly:\n", list(poly)
                print "Chunk 1:\n ", list(chunk1)
                print "Chunk 2:\n ", list(chunk2)
                print moments(chunk1)
                print moments(chunk2)
                print "Vol check: %g + %g = %g" % (v1, v2, v0)
                writePolyOBJ(poly, "poly.obj")
                writePolyOBJ(chunk1, "chunk1.obj")
                writePolyOBJ(chunk2, "chunk2.obj")
            self.failUnless(success, "Plane clipping summing to wrong volumes: %s + %s = %s != %s" % (v1, v2, v1 + v2, v0))

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
                    writePolyOBJ(poly, "poly.obj")
                    writePolyOBJ(chunk1, "chunk1.obj")
                    writePolyOBJ(chunk2, "chunk2.obj")
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
                self.failUnless(len(chunk) == 0,
                                "Full plane clipping failure: %s" % chunk)

    #---------------------------------------------------------------------------
    # Clip with planes that exactly hit at least one vertex.  Should reject
    # entire polygon.
    #---------------------------------------------------------------------------
    def testFullClipOnePlaneDegenerate(self):
        for points in (square_points, diamond_points, degenerate_square_points):
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            for v in poly:
                phat = vertexNormal(v, poly)
                planes = [Plane2d(v.position, phat)]
                chunk = Polygon(poly)
                clipPolygon(chunk, planes)
                self.failUnless(len(chunk) == 0,
                                "Full plane clipping failure: %s" % chunk)

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


    #---------------------------------------------------------------------------
    # Intersect two triangles with a small intersection
    #---------------------------------------------------------------------------
    def testIntersectTriangles(self):
        # Make the clipping polygon
        clipperTriVertices = [
           Vector2d( 1.028506494059886700E+00,  2.999999990617289169E-01),
           Vector2d( 1.028506493432178370E+00,  1.999999996863228713E-01),
           Vector2d( 1.067486385958200890E+00,  2.499999991021088219E-01),
        ]
        clipperTri = Polygon()
        initializePolygon(clipperTri, clipperTriVertices, [[2,1], [0,2], [1,0], ])

        # Make the polygon that will be clipped
        targetTriVertices = [
           Vector2d( 1.011431130628088759E+00,  2.999999993207053683E-01),
           Vector2d( 1.006141848526197480E+00,  2.999999991474828209E-01),
           Vector2d( 1.106466278308764428E+00,  2.999999978793384536E-01),
        ]
        targetTri = Polygon()
        initializePolygon(targetTri, targetTriVertices, [[2,1], [0,2], [1,0], ])

        # Now turn the clipper triangle into list of clip planes
        clipPlanes = []
        clipperTriArea, clipperTriCentroid = moments(clipperTri)
        for f in extractFaces(clipperTri):
            v0 = clipperTri[f[0]].position
            v1 = clipperTri[f[1]].position
            edgeVector  = v1 - v0
            halfEdgePos = 0.5*(v1 + v0)
            unitVector = Vector2d( edgeVector.y, -edgeVector.x ).unitVector()

            # Make sure the clipperTri face normal vectors point inward
            centroidDir = clipperTriCentroid - halfEdgePos
            if unitVector.dot(centroidDir) < 0.0:
               unitVector *= -1.0

            clipPlanes.append( Plane2d(v0, unitVector) )

        # Test passes if clipping is successful
        clipPolygon(targetTri, clipPlanes)

    #---------------------------------------------------------------------------
    # Create a pathological polygon by clipping, and then call collapseDegenerates
    #---------------------------------------------------------------------------
    def testIntersectAndCollapse(self):

        # Make the clipping polygon
        clipperTriVertices = [
           Vector2d( 5.597115204949861672E-02, -9.699999999948993867E-01),
           Vector2d( 6.000000004947565052E-02, -9.699999999914215021E-01),
           Vector2d( 5.519423042922926015E-02, -9.659999999884136823E-01),
        ]
        clipperTri = Polygon()
        initializePolygon(clipperTri, clipperTriVertices, [[2,1], [0,2], [1,0], ])

        # Make the polygon that will be clipped
        targetTriVertices = [
           Vector2d( 5.000000000142665324E-02, -9.700000000000538192E-01),
           Vector2d( 6.000000024821047079E-02, -9.699999999505402037E-01),
           Vector2d( 5.000000001436286784E-02, -9.699999999987333199E-01),
        ]
        targetTri = Polygon()
        initializePolygon(targetTri, targetTriVertices, [[2,1], [0,2], [1,0], ])

        # Now turn the clipper triangle into list of clip planes
        clipperTriArea, clipperTriCentroid = moments(clipperTri)

        v0 = clipperTriVertices[0]
        v1 = clipperTriVertices[1]
        v2 = clipperTriVertices[2]

        vc = (v0 + v1 + v2)*(1.0/3.0)

        def getClipPlane(v0, v1, vc):
            edgeVector  = v1 - v0
            halfEdgePos = 0.5*(v1 + v0)
            unitVector = Vector2d( edgeVector.y, -edgeVector.x ).unitVector()

            # Make sure the clipperTri face normal vectors point inward
            centroidDir = vc - halfEdgePos
            if unitVector.dot(centroidDir) < 0.0:
               unitVector *= -1.0

            return Plane2d(v0, unitVector)

        # note the funny order here. The error mode depends on the order.
        clipPlanes = []
        clipPlanes.append(getClipPlane(v0, v1, vc))
        clipPlanes.append(getClipPlane(v2, v0, vc))
        clipPlanes.append(getClipPlane(v1, v2, vc))

        # clip the polygon
        clipPolygon(targetTri, clipPlanes)

        # test passes as long as collapseDegenerates doesn't throw
        collapseDegenerates(targetTri, 1.0e-12)

    #---------------------------------------------------------------------------
    # Split a (convex) polygon into triangles.
    #---------------------------------------------------------------------------
    def testSplitIntoTriangles(self):
        for points in self.convexPointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            tris = splitIntoTriangles(poly)
            vol0, centroid0 = moments(poly)
            volTris = 0.0
            centroidTris = Vector2d()
            for inds in tris:
                assert len(inds) == 3
                a = (poly[inds[1]].position - poly[inds[0]].position).crossmag(poly[inds[2]].position - poly[inds[0]].position)
                assert a >= 0.0
                volTris += a
                centroidTris += a*(poly[inds[0]].position + poly[inds[1]].position + poly[inds[2]].position)
            volTris *= 0.5
            centroidTris /= 6.0*volTris
            assert abs(volTris - vol0) < 1.0e-20
            assert (centroidTris - centroid0).magnitude() < 1.0e-20

    #---------------------------------------------------------------------------
    # extractFaces
    #---------------------------------------------------------------------------
    def testExtractFaces(self):
        for points in self.pointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            n = len(points)
            answer = [[i, (i+1)%n] for i in xrange(n)]
            faces = extractFaces(poly)
            assert len(faces) == len(answer)
            for face in faces:
                assert len(face) == 2
                assert face in answer

    #---------------------------------------------------------------------------
    # commonFaceClips (convex)
    #---------------------------------------------------------------------------
    def testCommonFaceClipsConvex(self):
        for points in self.convexPointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            for i in xrange(self.ntests):
                p0 = Vector2d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                norm1 = Vector2d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector2d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = [Plane2d(p0,  norm1, 10),
                           Plane2d(p0,  norm2, 20)]
                chunk1 = Polygon(poly)
                clipPolygon(chunk1, planes1)
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
        for points in self.nonconvexPointSets:
            poly = Polygon()
            initializePolygon(poly, points, vertexNeighbors(points))
            for i in xrange(self.ntests):
                p0 = Vector2d(rangen.uniform(0.0, 1.0),
                              rangen.uniform(0.0, 1.0))
                norm1 = Vector2d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                norm2 = Vector2d(rangen.uniform(-1.0, 1.0), 
                                 rangen.uniform(-1.0, 1.0)).unitVector()
                planes1 = [Plane2d(p0,  norm1, 10),
                           Plane2d(p0,  norm2, 20)]
                chunk1 = Polygon(poly)
                clipPolygon(chunk1, planes1)
                faces1 = extractFaces(chunk1)
                clips1 = commonFaceClips(chunk1, faces1)
                assert len(clips1) == len(faces1)
                for clip in clips1:
                    assert len(clip) in (0, 1, 2)
                    for iclip in clip:
                        assert iclip in (10, 20)

    #---------------------------------------------------------------------------
    # Pathological degeneracy reported by Branson
    #---------------------------------------------------------------------------
    def testDegenerateClip1(self):
        targetTriVertices = [
            Vector2d( 6.104166666666669405E-01,  1.000000000000000000E+00),
            Vector2d( 6.100000000000003197E-01,  1.000000000000000000E+00),
            Vector2d( 6.104166666666669405E-01,  0.000000000000000000E+00),
        ]
        targetTri = Polygon()
        initializePolygon(targetTri, targetTriVertices, [[2,1], [0,2], [1,0], ])
        vol0, centroid0 = moments(targetTri)
        clipPlanes = [
            Plane2d(6.103886381791467919E-01, Vector2d(-9.999540830238576872E-01, -9.582893295645468490E-03)),
        ]
        clipPolygon(targetTri, clipPlanes)
        vol1, centroid1 = moments(targetTri)
        self.failUnless(vol1 == 0 and centroid1 == Vector2d(0,0),
                        "Degenerate test 1 failed: %s %s %s %s" % (vol0, centroid0, vol1, centroid1))

if __name__ == "__main__":
    unittest.main()
