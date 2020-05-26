import sys
sys.path.append("../test")

from math import *
from PolyClipper import *
from PolyClipperTestUtilities import *
from testPolyClipper2d import *

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
n = len(notched_points)
notched_neighbors = [[(i - 1) % n, (i + 1) % n] for i in xrange(n)]

#-------------------------------------------------------------------------------
# Unclipped object
#-------------------------------------------------------------------------------
poly = Polygon()
initializePolygon(poly, notched_points, notched_neighbors)
writePolyOBJ(poly, "notched_polygon.obj")
print "Starting poly: ", poly

#-------------------------------------------------------------------------------
# Clip 1
#-------------------------------------------------------------------------------
initializePolygon(poly, notched_points, notched_neighbors)
planes = [Plane2d(2.5, Vector2d(-1, 0.5).unitVector(), 10)]
clipPolygon(poly, planes)
writePolyOBJ(poly, "notched_polygon_clip1.obj")
print "Single clip: ", poly

#-------------------------------------------------------------------------------
# Clip 2
#-------------------------------------------------------------------------------
initializePolygon(poly, notched_points, notched_neighbors)
planes = [Plane2d(-2.5, Vector2d(1, -0.5).unitVector(), 20)]
clipPolygon(poly, planes)
writePolyOBJ(poly, "notched_polygon_clip2.obj")
print "Reverse clip: ", poly

#-------------------------------------------------------------------------------
# Clip 3
#-------------------------------------------------------------------------------
initializePolygon(poly, notched_points, notched_neighbors)
planes = [Plane2d(2.5, Vector2d(-1, 0.5).unitVector(), 10),
          Plane2d(-1.5, Vector2d(1, 5).unitVector(), 30)]
clipPolygon(poly, planes)
writePolyOBJ(poly, "notched_polygon_clip3.obj")
print "Double clip: ", poly
