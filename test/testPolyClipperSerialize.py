import unittest, random, string
from PolyClipper import *

# Create a global random number generator.
rangen = random.Random()

ntests = 100

#-------------------------------------------------------------------------------
# Serialize/deserialize unit tests
#-------------------------------------------------------------------------------
class TestPolyClipperSerialize(unittest.TestCase):

    # double
    def test_double(self):
        vals0 = [rangen.uniform(-100.0, 100.0) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_double(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_double(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # int
    def test_int(self):
        vals0 = [rangen.randint(-100, 100) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_int(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_int(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # size_t
    def test_size_t(self):
        vals0 = [rangen.randint(0, 1000) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_size_t(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_size_t(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # string
    def test_string(self):
        stuff = string.ascii_letters + string.digits + string.punctuation
        vals0 = [''.join(rangen.choice(stuff) for j in xrange(20)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_string(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_string(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Vector2d
    def test_Vector2d(self):
        vals0 = [Vector2d(rangen.uniform(-100.0, 100.0),
                          rangen.uniform(-100.0, 100.0)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Vector2d(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Vector2d(itr, buf)
            vals1.append(x)
        assert vals1 == vals0


    # Vector3d
    def test_Vector3d(self):
        vals0 = [Vector3d(rangen.uniform(-100.0, 100.0),
                          rangen.uniform(-100.0, 100.0),
                          rangen.uniform(-100.0, 100.0)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Vector3d(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Vector3d(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Vertex2d
    def test_Vertex2d(self):
        vals0 = [Vertex2d(pos = Vector2d(rangen.uniform(-100.0, 100.0),
                                         rangen.uniform(-100.0, 100.0)),
                          c = rangen.randint(-1, -1)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Vertex2d(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Vertex2d(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Vertex3d
    def test_Vertex3d(self):
        vals0 = [Vertex3d(pos = Vector3d(rangen.uniform(-100.0, 100.0),
                                         rangen.uniform(-100.0, 100.0),
                                         rangen.uniform(-100.0, 100.0)),
                          c = rangen.randint(-1, -1)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Vertex3d(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Vertex3d(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Polygon (these will be toplogically invalid, but OK for our test)
    def test_Polygon(self):
        vals0 = [Polygon([Vertex2d(pos = Vector2d(rangen.uniform(-100.0, 100.0),
                                                  rangen.uniform(-100.0, 100.0)),
                                   c = rangen.randint(-1, -1)) for i in xrange(ntests)]) for j in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Polygon(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Polygon(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Polyhedron (these will be toplogically invalid, but OK for our test)
    def test_Polyhedron(self):
        vals0 = [Polyhedron([Vertex3d(pos = Vector3d(rangen.uniform(-100.0, 100.0),
                                                     rangen.uniform(-100.0, 100.0),
                                                     rangen.uniform(-100.0, 100.0)),
                                      c = rangen.randint(-1, -1)) for i in xrange(ntests)]) for j in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Polyhedron(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Polyhedron(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Plane2d
    def test_Plane2d(self):
        vals0 = [Plane2d(p = Vector2d(rangen.uniform(-100.0, 100.0),
                                      rangen.uniform(-100.0, 100.0)),
                         nhat = Vector2d(rangen.uniform(-1.0, 1.0),
                                         rangen.uniform(-1.0, 1.0)).unitVector(),
                         id = rangen.randint(0, 1000)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Plane2d(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Plane2d(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # Plane3d
    def test_Plane3d(self):
        vals0 = [Plane3d(p = Vector3d(rangen.uniform(-100.0, 100.0),
                                      rangen.uniform(-100.0, 100.0),
                                      rangen.uniform(-100.0, 100.0)),
                         nhat = Vector3d(rangen.uniform(-1.0, 1.0),
                                         rangen.uniform(-1.0, 1.0),
                                         rangen.uniform(-1.0, 1.0)).unitVector(),
                         id = rangen.randint(0, 1000)) for i in xrange(ntests)]
        buf = vector_of_char()
        for x in vals0:
            serialize_Plane3d(x, buf)
        itr = 0
        vals1 = []
        for i in xrange(ntests):
            x, itr = deserialize_Plane3d(itr, buf)
            vals1.append(x)
        assert vals1 == vals0

    # # std::vector<Plane2d>
    # def test_vector_of_Plane2d(self):
    #     vals0 = [vector_of_Plane2d([Plane2d(p = Vector2d(rangen.uniform(-100.0, 100.0),
    #                                                      rangen.uniform(-100.0, 100.0)),
    #                                         nhat = Vector2d(rangen.uniform(-1.0, 1.0),
    #                                                         rangen.uniform(-1.0, 1.0)).unitVector(),
    #                                         id = rangen.randint(0, 1000)) for i in xrange(ntests)]) for j in xrange(ntests)]
    #     buf = vector_of_char()
    #     for x in vals0:
    #         serialize_vector_of_Plane2d(x, buf)
    #     itr = 0
    #     vals1 = []
    #     for i in xrange(ntests):
    #         x, itr = deserialize_vector_of_Plane2d(itr, buf)
    #         vals1.append(x)
    #     assert vals1 == vals0

    # # std::vector<Plane3d>
    # def test_vector_of_Plane3d(self):
    #     vals0 = [vector_of_Plane3d([Plane3d(p = Vector3d(rangen.uniform(-100.0, 100.0),
    #                                                      rangen.uniform(-100.0, 100.0),
    #                                                      rangen.uniform(-100.0, 100.0)),
    #                                         nhat = Vector3d(rangen.uniform(-1.0, 1.0),
    #                                                         rangen.uniform(-1.0, 1.0),
    #                                                         rangen.uniform(-1.0, 1.0)).unitVector(),
    #                                         id = rangen.randint(0, 1000)) for i in xrange(ntests)]) for j in xrange(ntests)]
    #     buf = vector_of_char()
    #     for x in vals0:
    #         serialize_vector_of_Plane3d(x, buf)
    #     itr = 0
    #     vals1 = []
    #     for i in xrange(ntests):
    #         x, itr = deserialize_vector_of_Plane3d(itr, buf)
    #         vals1.append(x)
    #     assert vals1 == vals0

#-------------------------------------------------------------------------------
# Execute the tests
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    unittest.main()
