import PolyClipper
from collections.abc import Iterable

#-------------------------------------------------------------------------------
# Fuzzy comparisons.
#-------------------------------------------------------------------------------
def fuzzyEqual(lhs, rhs,
               fuzz = 1.0e-5):
    if isinstance(lhs, Iterable):
        assert isinstance(rhs, Iterable) and len(lhs) == len(rhs)
        return min([fuzzyEqual(x, y, fuzz) for (x, y) in zip(lhs, rhs)])
    else:
        return abs(lhs - rhs)/max(1.0, abs(lhs) + abs(rhs)) < fuzz;

#-------------------------------------------------------------------------------
# Write an OBJ (vertex-facet labeled) shape file from a polygon/polyhedron.
# NOTE!!  These things seem to count from 1 rather than 0 for indices!
#-------------------------------------------------------------------------------
def writePolyOBJ(poly, filename, forceTriangles=False):
    with open(filename, "w") as f:
        facets = PolyClipper.extractFaces(poly)

        if isinstance(poly, PolyClipper.Polygon):
            for v in poly:
                f.write("v %g %g 0.0\n" % (v.position.x, v.position.y))

            # There could be multiple loops in a polygon, so break
            # those up into individual facets for obj output.
            # First, sort topologically so the loops are contiguous
            for i in range(len(facets) - 1):
                v1 = facets[i][1]
                j = i + 1
                while j < len(facets) and v1 != facets[j][0]:
                    j += 1
                if j < len(facets):
                    facets[i+1], facets[j] = facets[j], facets[i+1]

            # Now we can write the facets.
            f.write("f")
            i = 0
            for i in range(len(facets)):
                f.write(" %i" % (facets[i][0] + 1))
                if ((i + 1) < len(facets) and
                    facets[i][1] != facets[i+1][0]):
                    f.write("\nf")
            f.write("\n")

        else:
            for v in poly:
                f.write("v %g %g %g\n" % (v.position.x, v.position.y, v.position.z))

            for facet in facets:
                f.write("f")
                if forceTriangles and len(facet > 3):
                    for j in range(1, len(facet)-1):
                        f.write(" %i %i %i\n" % (0, j, j+1))
                else:
                    for i in facet:
                        f.write(" %i" % (i + 1))
                f.write("\n")

    return

