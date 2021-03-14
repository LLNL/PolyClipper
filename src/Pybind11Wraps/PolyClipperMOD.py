"""
PolyClipper module.

Binds the PolyClipper geometry operations for clipping polygons & polyhedra.
"""

from PYB11Generator import *
import types

# Include files.
PYB11includes = ['"polyclipper2d.hh"',
                 '"polyclipper3d.hh"',
                 '"polyclipper_vector2d.hh"',
                 '"polyclipper_vector3d.hh"',
                 '"polyclipper_plane.hh"']

PYB11namespaces = ["PolyClipper"]

PYB11preamble = """
using Polygon = std::vector<PolyClipper::Vertex2d<>>;
using Polyhedron = std::vector<PolyClipper::Vertex3d<>>;
using Plane2d = PolyClipper::Plane<PolyClipper::internal::VectorAdapter<PolyClipper::Vector2d>>;
using Plane3d = PolyClipper::Plane<PolyClipper::internal::VectorAdapter<PolyClipper::Vector3d>>;
"""

#-------------------------------------------------------------------------------
# Import the types defined in this module
#-------------------------------------------------------------------------------
from Vector2d import *
from Vector3d import *
from Vertex2d import *
from Vertex3d import *
from Plane import *

#-------------------------------------------------------------------------------
# Polygon & Polyhedron
#-------------------------------------------------------------------------------
Polygon    = PYB11_bind_vector("PolyClipper::Vertex2d<>", opaque=True, local=True)
Polyhedron = PYB11_bind_vector("PolyClipper::Vertex3d<>", opaque=True, local=True)

# We also use vector<char>
vector_of_char = PYB11_bind_vector("char", opaque=True, local=True)

#-------------------------------------------------------------------------------
# Plane
#-------------------------------------------------------------------------------
Plane2d = PYB11TemplateClass(Plane, template_parameters="internal::VectorAdapter<Vector2d>")
Plane3d = PYB11TemplateClass(Plane, template_parameters="internal::VectorAdapter<Vector3d>")

#-------------------------------------------------------------------------------
# Polygon methods.
#-------------------------------------------------------------------------------
def initializePolygon(poly = "Polygon&",
                      positions = "const std::vector<Vector2d>&",
                      neighbors = "const std::vector<std::vector<int>>&"):
    "Initialize a PolyClipper::Polygon from vertex positions and vertex neighbors."
    return "void"

@PYB11cppname("polygon2string<>")
def polygon2string(poly = "Polygon&"):
    "Return a formatted string representation for a PolyClipper::Polygon."
    return "std::string"

@PYB11implementation("""[](const Polygon& self) {
                                                  double zerothMoment;
                                                  Vector2d firstMoment;
                                                  moments(zerothMoment, firstMoment, self);
                                                  return py::make_tuple(zerothMoment, firstMoment);
                                                }""")
@PYB11pycppname("moments")
def momentsPolygon(poly = "const Polygon&"):
    "Compute the zeroth and first moment of a PolyClipper::Polygon."
    return "py::tuple"

def clipPolygon(poly = "Polygon&",
                planes = "const std::vector<Plane2d>&"):
    "Clip a PolyClipper::Polygon with a collection of planes."
    return "void"

@PYB11pycppname("collapseDegenerates")
def collapseDegeneratesPolygon(poly = "Polygon&",
                               tol = "const double"):
    "Collapse edges in a PolyClipper::Polygon below the given tolerance."
    return "void"

def extractFaces(poly = "const Polygon&"):
    "Compute the faces (as pairs of vertex indices) for the Polygon"
    return "std::vector<std::vector<int>>"

def commonFaceClips(poly = "const Polygon&",
                    faces = "const std::vector<std::vector<int>>&"):
    "Find the common clipping planes for each face"
    return "std::vector<std::set<int>>"

def splitIntoTriangles(poly = "const Polygon&",
                       tol = ("const double", "0.0")):
    """Split a PolyClipper::Polygon into triangles.
The result is returned as a vector<vector<int>>, where each inner vector is a triple of
ints representing vertex indices in the input Polygon."""
    return "std::vector<std::vector<int>>"

#-------------------------------------------------------------------------------
# Polyhedron methods.
#-------------------------------------------------------------------------------
def initializePolyhedron(poly = "Polyhedron&",
                         positions = "const std::vector<Vector3d>&",
                         neighbors = "const std::vector<std::vector<int>>&"):
    "Initialize a PolyClipper::Polyhedron from vertex positions and vertex neighbors."
    return "void"

@PYB11cppname("polyhedron2string<>")
def polyhedron2string(poly = "Polyhedron&"):
    "Return a formatted string representation for a PolyClipper::Polyhedron."
    return "std::string"

@PYB11implementation("""[](const Polyhedron& self) {
                                                     double zerothMoment;
                                                     Vector3d firstMoment;
                                                     moments(zerothMoment, firstMoment, self);
                                                     return py::make_tuple(zerothMoment, firstMoment);
                                                   }""")
@PYB11pycppname("moments")
def momentsPolyhedron(poly = "const Polyhedron&"):
    "Compute the zeroth and first moment of a PolyClipper::Polyhedron."
    return "py::tuple"

def clipPolyhedron(poly = "Polyhedron&",
                   planes = "const std::vector<Plane3d>&"):
    "Clip a PolyClipper::Polyhedron with a collection of planes."
    return "void"

@PYB11pycppname("collapseDegenerates")
def collapseDegeneratesPolyhedron(poly = "Polyhedron&",
                                  tol = "const double"):
    "Collapse edges in a PolyClipper::Polyhedron below the given tolerance."
    return "void"

@PYB11pycppname("extractFaces")
def extractFacesPolyhedron(poly = "const Polyhedron&"):
    "Compute the faces (as pairs of vertex indices) for the Polyhedron"
    return "std::vector<std::vector<int>>"

@PYB11pycppname("commonFaceClips")
def commonFaceClipsPolyhedron(poly = "const Polyhedron&",
                              faces = "const std::vector<std::vector<int>>&"):
    "Find the common clipping planes for each face"
    return "std::vector<std::set<int>>"

def splitIntoTetrahedra(poly = "const Polyhedron&",
                        tol = ("const double", "0.0")):
    """Split a PolyClipper::Polyhedron into tetrahedra.
The result is returned as a vector<vector<int>>, where each inner vector is a set of four
ints representing vertex indices in the input Polyhedron."""
    return "std::vector<std::vector<int>>"

#-------------------------------------------------------------------------------
# Serialization
#-------------------------------------------------------------------------------
@PYB11cppname("internal::serialize")
def serialize_double(val = "const double",
                     buffer = "std::vector<char>&"):
    "Serialize a double"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_int(val = "const int",
                     buffer = "std::vector<char>&"):
    "Serialize an int"
    return "void"

@PYB11implementation("[](const int val, std::vector<char>& buffer) { size_t uval = size_t(val); internal::serialize(uval, buffer); }")
def serialize_size_t(val = "const int",
                     buffer = "std::vector<char>&"):
    "Serialize a size_t (via Python int)"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_string(val = "const std::string&",
                     buffer = "std::vector<char>&"):
    "Serialize a std::string"
    return "void"

@PYB11cppname("internal::serialize<internal::VectorAdapter<Vector2d>>")
def serialize_Vector2d(val = "const Vector2d&",
                       buffer = "std::vector<char>&"):
    "Serialize a Vector2d"
    return "void"

@PYB11cppname("internal::serialize<internal::VectorAdapter<Vector3d>>")
def serialize_Vector3d(val = "const Vector3d&",
                       buffer = "std::vector<char>&"):
    "Serialize a Vector3d"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_Vertex2d(val = "const Vertex2d<>&",
                       buffer = "std::vector<char>&"):
    "Serialize a Vertex2d"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_Vertex3d(val = "const Vertex3d<>&",
                       buffer = "std::vector<char>&"):
    "Serialize a Vertex3d"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_Polygon(val = "const std::vector<Vertex2d<>>&",
                       buffer = "std::vector<char>&"):
    "Serialize a Polygon"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_Polyhedron(val = "const std::vector<Vertex3d<>>&",
                         buffer = "std::vector<char>&"):
    "Serialize a Polyhedron"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_Plane2d(val = "const Plane2d&",
                      buffer = "std::vector<char>&"):
    "Serialize a Plane2d"
    return "void"

@PYB11cppname("internal::serialize")
def serialize_Plane3d(val = "const Plane3d&",
                      buffer = "std::vector<char>&"):
    "Serialize a Plane3d"
    return "void"

#-------------------------------------------------------------------------------
# Deserialization
# I'll get cute here and use python code generation to make the various versions
#-------------------------------------------------------------------------------
for (cppname, pyname, template_arg) in (("double", "double", ""),
                                        ("int", "int", ""),
                                        ("size_t", "size_t", ""),
                                        ("std::string", "string", ""),
                                        ("Vector2d", "Vector2d", "<internal::VectorAdapter<Vector2d>>"),
                                        ("Vector3d", "Vector3d", "<internal::VectorAdapter<Vector3d>>"),
                                        ("Vertex2d<>", "Vertex2d", ""),
                                        ("Vertex3d<>", "Vertex3d", ""),
                                        ("Polygon", "Polygon", ""),
                                        ("Polyhedron", "Polyhedron", ""),
                                        ("Plane2d", "Plane2d", "<internal::VectorAdapter<Vector2d>>"),
                                        ("Plane3d", "Plane3d", "<internal::VectorAdapter<Vector3d>>")):
    exec("""
@PYB11implementation("[](int itr_pos, const std::vector<char>& buffer) {{ auto itr = buffer.begin() + itr_pos; {cppname} val; internal::deserialize{template_arg}(val, itr, buffer.end()); return py::make_tuple(val, int(std::distance(buffer.begin(), itr))); }}")
def deserialize_{pyname}(itr_pos = "int",
                    buffer = "const std::vector<char>&"):
    "Deserialize {cppname} starting at itr_pos in buffer.  Returns ({cppname}, new_itr_pos)"
    return "py::tuple"

""".format(cppname = cppname,
           pyname = pyname,
           template_arg = template_arg))

