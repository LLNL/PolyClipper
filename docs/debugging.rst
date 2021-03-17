############################################
A few notes about debugging with PolyClipper
############################################

PolyClipper can be built for debugging using the standard Cmake configuration ``-DCMAKE_BUILD_TYPE=Debug``.  When configured this way, this activates a custom assertion method in PolyClipper (``PCASSERT``)  defined in ``polyclipper_utilities.hh``.  Since PolyClipper in C++ is entirely inlined header functions, these assertions are active for any compilation using including PolyClipper headers that do not define the ``-DNDEBUG`` option on the compile line.

In general compiling with these assertions active will have a serious runtime cost.  However, when trying to track down bugs or problems using PolyCliperr code these internal checks can be invaluable.  There are two aspects in particular developers should be aware of.

* PolyClipper defines a custom exception type, ``PolyClipperError``, which can be screened for with standard C++ try/catch clauses.  This is the exception that will be raised by any ``PCASSERT`` violation.

* PolyClipper includes some handy serialization methods for packaging up polygons/polyhedra and planes when exceptions occur.  Many of the ``PCASSERT`` statements will result in this serialized state being written out to a local file (in binary format), which can then be picked up by standalone PolyClipper to reproduce the error condition.  The source file ``test/stuff/test_bad_clip.py`` is one example of how to pick up one of these serialized binary blobs using standalone PolyClipper's Python interface to reproduce 3D polyhedral clipping problems.  When reporting problems with PolyClipper these blobs capturing the exact input that reproduces the error can be extremely useful in getting such problems resolved.

