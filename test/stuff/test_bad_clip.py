import os, os.path
from PolyClipper import *

for filename in os.listdir("."):
    if filename.endswith(".bin"):
        print "--------------------------------------------------------------------------------"
        print "Reading ", filename
        with open(filename, 'rb') as f:
            binstuff = f.read()
        buf  = vector_of_char(binstuff)

        poly, itr = deserialize_Polyhedron(0, buf)
        planes, itr = deserialize_vector_of_Plane3d(itr, buf)

        print "Read poly:\n", polyhedron2string(poly)
        print "  moments: ", moments(poly)

        clipPolyhedron(poly, planes)

        print "After clipping:\n", polyhedron2string(poly)
        print "  moments: ", moments(poly)
