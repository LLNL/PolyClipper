from collections import Iterable

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
