PolyClipper
==============

PolyClipper is a C++ reimplementation of the geometric clipping operations in the [R3D](https://github.com/devonmpowell/r3d) library originally written by Devon Powell, as documented in the paper
[Powell & Abell (2015)](http://www.sciencedirect.com/science/article/pii/S0021999115003563).

The main focus here is on clipping polygons (in 2D (x,y) coordinates) and polyhedra (in 3D (x,y,z) coordinates) with planes, returning new polygons/polyhedra as the result of this clipping.  The input polygons/polyhedra may be non-convex and arbitrarily complex, but the only clipping operation supported is with planes.  This is equivalent to intersecting one arbitrary (not necessarily convex) polygon/polyhedron with a convex polygon/polyhedron.

PolyClipper reimplements these clipping operations from R3D for two reasons:
  * PolyClipper removes the hard-coded size limitations of R3D on the number of vertices/complexity of the polygons and polyhedra.
  * PolyClipper also removes the assumption that each vertex in 3D has exactly three neighbors (as well as the related limitation of two neighbors in 2D) -- the number of neighbors per vertex is now arbitrary.  This also removes the complexity of requiring degenerate/redundant vertices.

Note PolyClipper currently does not provide the generalized voxelization or arbitrary integrals over polygons/polyhedra as provided in R3D.   These would be straightforward to add, but were not necessary for the authors needs from the library, which is to generalize the clipping algorithms.  The only method of this sort provided by PolyClipper is the ability to do the zeroth and first moment integrals over the polygons/polyhedra.

PolyClipper currently provides both C++ and Python interfaces.

Documentation
-------------

PolyClipper is documented at [readthedocs](https://polyclipper.readthedocs.io/en/latest/).

Note the source for this documentation is embedded in the PolyClipper repository under docs/.

Contributions
-------------

Contributions are welcome, and should be provided as pull requests to the main repository.  Note all contributions must be provided under the same license for distribution (in this case the BSD license).

License
-------

PolyClipper is released under the [BSD license](https://github.com/LLNL/PolyClipper/blob/master/LICENSE).

LLNL-CODE-811676

SPDX-License-Identifier: BSD-3
