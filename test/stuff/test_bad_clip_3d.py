import os, os.path
from PolyClipper import *

for filename in os.listdir("3d"):
    if filename.endswith(".bin"):
        print("--------------------------------------------------------------------------------")
        print("Reading ", filename)
        with open(os.path.join("3d", filename), 'rb') as f:
            buf = f.read()

        poly, itr = deserialize_Polyhedron(0, buf)
        planes, itr = deserialize_vector_of_Plane3d(itr, buf)

        print("Read poly:\n", polyhedron2string(poly))
        print("  moments: ", moments(poly))
        print("Read planes:\n", planes)

        clipPolyhedron(poly, planes)

        print("After clipping:\n", polyhedron2string(poly))
        print("  moments: ", moments(poly))
