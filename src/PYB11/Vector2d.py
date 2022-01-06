from PYB11Generator import *

class Vector2d:
    "The PolyClipper (x,y) geometric vector type"

    #---------------------------------------------------------------------------
    # Constructors
    #---------------------------------------------------------------------------
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                X = "double",
                Y = "double"):
        "Construct with (X,Y) values"

    #---------------------------------------------------------------------------
    # Operators
    #---------------------------------------------------------------------------
    def __eq__(self):
        return

    def __imul__(self, rhs="double()"):
        return

    def __idiv__(self, rhs="double()"):
        return

    def __iadd__(self):
        return

    def __isub__(self):
        return

    def __mul__(self, rhs="double()"):
        return

    def __rmul__(self, rhs="double()"):
        return

    def __div__(self, rhs="double()"):
        return

    def __add__(self):
        return

    def __sub__(self):
        return

    def __neg__(self):
        return

    @PYB11implementation("[](Vector2d &s) { return py::make_iterator(&s.x, &s.y+1u); }, py::keep_alive<0,1>()")
    def __iter__(self):
        "Python iteration through a Vector."

    #---------------------------------------------------------------------------
    # Methods
    #---------------------------------------------------------------------------
    @PYB11const
    def dot(self, rhs="const Vector2d&"):
        "Return the dot product"
        return "double"

    @PYB11const
    def crossmag(self, rhs="const Vector2d&"):
        "Return the z-value of the cross product as a scalar"
        return "double"

    @PYB11const
    def magnitude(self):
        "Return the magnitude of the vector"
        return "double"

    @PYB11const
    def magnitude2(self):
        "Return the square of the magnitude of the vector"
        return "double"

    @PYB11const
    def unitVector(self):
        "Return a unit vector in the same direction as this one"
        return "Vector2d"

    @PYB11implementation('''[](const Vector2d& self) { return "(" + std::to_string(self.x) + ", " + std::to_string(self.y) + ")"; }''')
    def __repr__(self):
        return

    #---------------------------------------------------------------------------
    # Attributes
    #---------------------------------------------------------------------------
    x = PYB11readwrite()
    y = PYB11readwrite()
