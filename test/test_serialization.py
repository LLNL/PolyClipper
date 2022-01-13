#ATS:test(SELF, label="PolyClipper serialization tests")

import unittest
import string

from PolyClipper import *

# Create a global random number generator.
import random
rangen = random.Random()

#-------------------------------------------------------------------------------
class TestPolyClipperSerialization(unittest.TestCase):

    #---------------------------------------------------------------------------
    # double
    #---------------------------------------------------------------------------
    def test_double(self):
        buf = bytes()
        x = rangen.uniform(-1e100, 1e100)
        buf = serialize_double(x, buf)
        y, bufend = deserialize_double(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # int
    #---------------------------------------------------------------------------
    def test_int(self):
        buf = bytes()
        x = rangen.randint(-2**30, 2**30)
        buf = serialize_int(x, buf)
        y, bufend = deserialize_int(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # string
    #---------------------------------------------------------------------------
    def test_string(self):
        buf = bytes()
        x = "".join(random.choices(string.ascii_uppercase +
                                   string.ascii_lowercase +
                                   string.digits, k=1024))
        buf = serialize_string(x, buf)
        y, bufend = deserialize_string(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Vector2d
    #---------------------------------------------------------------------------
    def test_Vector2d(self):
        buf = bytes()
        x = Vector2d(rangen.uniform(-1e100, 1e100),
                     rangen.uniform(-1e100, 1e100))
        buf = serialize_Vector2d(x, buf)
        y, bufend = deserialize_Vector2d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Vector3d
    #---------------------------------------------------------------------------
    def test_Vector3d(self):
        buf = bytes()
        x = Vector3d(rangen.uniform(-1e100, 1e100),
                     rangen.uniform(-1e100, 1e100),
                     rangen.uniform(-1e100, 1e100))
        buf = serialize_Vector3d(x, buf)
        y, bufend = deserialize_Vector3d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Vertex2d
    #---------------------------------------------------------------------------
    def test_Vertex2d(self):
        buf = bytes()
        x = Vertex2d()
        x.position = Vector2d(rangen.uniform(-1e100, 1e100),
                              rangen.uniform(-1e100, 1e100))
        x.neighbors = (1, 2)
        x.comp = rangen.choice((-2, -1, 0, 1, 2))
        x.ID = rangen.randint(-2**30, 2**30)
        buf = serialize_Vertex2d(x, buf)
        y, bufend = deserialize_Vertex2d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Vertex3d
    #---------------------------------------------------------------------------
    def test_Vertex3d(self):
        buf = bytes()
        x = Vertex3d()
        x.position = Vector3d(rangen.uniform(-1e100, 1e100),
                              rangen.uniform(-1e100, 1e100),
                              rangen.uniform(-1e100, 1e100))
        x.neighbors = [1, 2, 10, 100, 14]
        x.comp = rangen.choice((-2, -1, 0, 1, 2))
        x.ID = rangen.randint(-2**30, 2**30)
        buf = serialize_Vertex3d(x, buf)
        y, bufend = deserialize_Vertex3d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Polygon
    #---------------------------------------------------------------------------
    def test_Polygon(self):
        buf = bytes()
        x = Polygon()
        N = 20
        for i in range(N):
            v = Vertex2d()
            v.position = Vector2d(rangen.uniform(-1e100, 1e100),
                                  rangen.uniform(-1e100, 1e100))
            v.neighbors = ((i - 1) % N,
                           (i + 1) % N)
            v.comp = rangen.choice((-2, -1, 0, 1, 2))
            v.ID = i
            x.append(v)
        assert len(x) == N
        buf = serialize_Polygon(x, buf)
        y, bufend = deserialize_Polygon(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Polyhedron
    #---------------------------------------------------------------------------
    def test_Polyhedron(self):
        buf = bytes()
        x = Polyhedron()
        N = 20
        for i in range(N):
            v = Vertex3d()
            v.position = Vector3d(rangen.uniform(-1e100, 1e100),
                                  rangen.uniform(-1e100, 1e100),
                                  rangen.uniform(-1e100, 1e100))
            v.neighbors = rangen.choices(range(N), k=3)
            v.comp = rangen.choice((-2, -1, 0, 1, 2))
            v.ID = i
            x.append(v)
        assert len(x) == N
        buf = serialize_Polyhedron(x, buf)
        y, bufend = deserialize_Polyhedron(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Plane2d
    #---------------------------------------------------------------------------
    def test_Plane2d(self):
        buf = bytes()
        x = Plane2d(p = Vector2d(rangen.uniform(-1e100, 1e100),
                                 rangen.uniform(-1e100, 1e100)),
                    nhat = Vector2d(rangen.uniform(-1e100, 1e100),
                                    rangen.uniform(-1e100, 1e100)).unitVector(),
                    id = rangen.randint(-2**30, 2**30))
        buf = serialize_Plane2d(x, buf)
        y, bufend = deserialize_Plane2d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Plane3d
    #---------------------------------------------------------------------------
    def test_Plane3d(self):
        buf = bytes()
        x = Plane3d(p = Vector3d(rangen.uniform(-1e100, 1e100),
                                 rangen.uniform(-1e100, 1e100),
                                 rangen.uniform(-1e100, 1e100)),
                    nhat = Vector3d(rangen.uniform(-1e100, 1e100),
                                    rangen.uniform(-1e100, 1e100),
                                    rangen.uniform(-1e100, 1e100)).unitVector(),
                    id = rangen.randint(-2**30, 2**30))
        buf = serialize_Plane3d(x, buf)
        y, bufend = deserialize_Plane3d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # vector<Plane2d>
    #---------------------------------------------------------------------------
    def test_vector_of_Plane2d(self):
        buf = bytes()
        N = 20
        x = [Plane2d(p = Vector2d(rangen.uniform(-1e100, 1e100),
                                  rangen.uniform(-1e100, 1e100)),
                     nhat = Vector2d(rangen.uniform(-1e100, 1e100),
                                     rangen.uniform(-1e100, 1e100)).unitVector(),
                     id = rangen.randint(-2**30, 2**30)) for i in range(N)]
        assert len(x) == N
        buf = serialize_vector_of_Plane2d(x, buf)
        y, bufend = deserialize_vector_of_Plane2d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

    #---------------------------------------------------------------------------
    # Plane3d
    #---------------------------------------------------------------------------
    def test_Plane3d(self):
        buf = bytes()
        N = 20
        x = [Plane3d(p = Vector3d(rangen.uniform(-1e100, 1e100),
                                  rangen.uniform(-1e100, 1e100),
                                  rangen.uniform(-1e100, 1e100)),
                     nhat = Vector3d(rangen.uniform(-1e100, 1e100),
                                     rangen.uniform(-1e100, 1e100),
                                     rangen.uniform(-1e100, 1e100)).unitVector(),
                     id = rangen.randint(-2**30, 2**30)) for i in range(N)]
        assert len(x) == N
        buf = serialize_vector_of_Plane3d(x, buf)
        y, bufend = deserialize_vector_of_Plane3d(0, buf)
        self.assertTrue(x == y,
                        "%s != %s" % (x, y))

#-------------------------------------------------------------------------------
# Run the tests
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
