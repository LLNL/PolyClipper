############
Introduction
############

PolyClipper is a C++ reimplementation of the geometric clipping operations in the `R3D <https://github.com/devonmpowell/r3d>`_ library originally written by Devon Powell, as documented in the paper
`Powell & Abell (2015) <http://www.sciencedirect.com/science/article/pii/S0021999115003563>`_.

PolyClipper reimplements these clipping operations for two reasons:
  * PolyClipper removes the hard-coded size limitations of R3D on the number of vertices/complexity of the polygons and polyhedra.
  * PolyClipper also removes the assumption that each vertex in 3D has exactly three neighbors (as well as the related limitation of two neighbors in 2D) -- the number of neighbors per vertex is now arbitrary.  This also removes the complexity of requiring degenerate/redundant vertices.

Note PolyClipper currently does not provide the generalized voxelization or arbitrary integrals over polygons/polyhedra as provided in R3D.   These would be straightforward to add, but were not necessary for my (i.e., Mike Owen's) needs from the library, which is to generalize the clipping algorithms.  The only method of this sort provided by PolyClipper is the ability to do the zeroth and first moment integrals over the polygons/polyhedra.
