import sys
sys.path.append("../test")

from math import *
from PolyClipper import *
from PolyClipperTestUtilities import *
from testPolyClipper3d import *

#-------------------------------------------------------------------------------
# Unclipped object
#-------------------------------------------------------------------------------
poly = Polyhedron()
initializePolyhedron(poly, notched_points, notched_neighbors)
writePolyOBJ(poly, "notched_polyhedron.obj")
print "Starting poly: ", list(poly)

#-------------------------------------------------------------------------------
# Clip 1
#-------------------------------------------------------------------------------
poly = Polyhedron()
initializePolyhedron(poly, notched_points, notched_neighbors)
planes = [Plane3d(Vector3d(3, 1, 0), Vector3d(-1, 0.5, -1.5).unitVector(), 10)]
clipPolyhedron(poly, planes)
writePolyOBJ(poly, "notched_polyhedron_clip1.obj")
print "Single clip: ", list(poly)

#-------------------------------------------------------------------------------
# Clip 2
#-------------------------------------------------------------------------------
poly = Polyhedron()
initializePolyhedron(poly, notched_points, notched_neighbors)
planes = [Plane3d(Vector3d(3, 1, 0), Vector3d(1, -0.5, 1.5).unitVector(), 10)]
clipPolyhedron(poly, planes)
writePolyOBJ(poly, "notched_polyhedron_clip2.obj")
print "Reverse clip: ", list(poly)

#-------------------------------------------------------------------------------
# Clip 3
#-------------------------------------------------------------------------------
poly = Polyhedron()
initializePolyhedron(poly, notched_points, notched_neighbors)
planes = [Plane3d(Vector3d(3, 1, 0), Vector3d(-1, 0.5, -1.5).unitVector(), 10),
          Plane3d(Vector3d(1, 1, 0), Vector3d(1, 0, -1).unitVector(), 30)]
clipPolyhedron(poly, planes)
writePolyOBJ(poly, "notched_polyhedron_clip3.obj")
print "Double clip: ", list(poly)
