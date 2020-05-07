TOP = /usr/gapps/Spheral/toss_3_x86_64_ib/exp
PYTHON = $(TOP)/bin/python
CXX = g++
AR = ar
RANLIB = ranlib
CXXFLAGS = -fpic -fexceptions -std=c++11 -march=native -I src -I$(PYTHONTOP)/include -I$(PYTHONTOP)/include/python2.7
OPTFLAGS = -O3 -DNDEBUG

all:	debug

SRC = polyclipper2d.cc polyclipper3d.cc
OBJS = $(subst .cc,.o,$(SRC))

%.o:	src/%.cc
	$(CXX) $(CXXFLAGS) -c $< -o $(*F).o

debug:	$(OBJS)
	$(AR) ru libPolyClipper.a $(OBJS)
	$(RANLIB) libPolyClipper.a

clean:
	rm -f *.o *.a
