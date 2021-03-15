//------------------------------------------------------------------------------
// Methods for serializing/dserializing PolyClipper types
//------------------------------------------------------------------------------
#ifndef __PolyClipper_serialize__
#define __PolyClipper_serialize__

#include "polyclipper_utilities.hh"

namespace PolyClipper {

// Forward declarations
template<typename VA> class Vertex2d;
template<typename VA> class Vertex3d;

namespace internal {

//------------------------------------------------------------------------------
// Serialize a double
//------------------------------------------------------------------------------
inline
void
serialize(const double val,
          std::vector<char>& buffer) {
  const auto n = sizeof(double);
  const char* data = reinterpret_cast<const char*>(&val);
  std::copy(data, data + n, std::back_inserter(buffer));
}

//------------------------------------------------------------------------------
// Serialize an int
//------------------------------------------------------------------------------
inline
void
serialize(const int val,
          std::vector<char>& buffer) {
  const auto n = sizeof(int);
  const char* data = reinterpret_cast<const char*>(&val);
  std::copy(data, data + n, std::back_inserter(buffer));
}

//------------------------------------------------------------------------------
// Serialize a size_t
//------------------------------------------------------------------------------
inline
void
serialize(const size_t val,
          std::vector<char>& buffer) {
  const auto n = sizeof(size_t);
  const char* data = reinterpret_cast<const char*>(&val);
  std::copy(data, data + n, std::back_inserter(buffer));
}

//------------------------------------------------------------------------------
// Serialize a std::string
//------------------------------------------------------------------------------
inline
void
serialize(const std::string& val,
          std::vector<char>& buffer) {
  const auto n = val.size();
  serialize(n, buffer);
  std::copy(val.begin(), val.end(), std::back_inserter(buffer));
}

//------------------------------------------------------------------------------
// Serialize a Vector
//------------------------------------------------------------------------------
template<typename VA>
inline
void
serialize(const typename VA::VECTOR& val,
          std::vector<char>& buffer) {
  const auto triple_vals = VA::get_triple(val);
  serialize(triple_vals[0], buffer);
  serialize(triple_vals[1], buffer);
  serialize(triple_vals[2], buffer);
}

//------------------------------------------------------------------------------
// Serialize a Vertex2d
//------------------------------------------------------------------------------
template<typename VA>
inline
void
serialize(const Vertex2d<VA>& val,
          std::vector<char>& buffer) {
  serialize<VA>(val.position, buffer);
  serialize(val.neighbors.first, buffer);
  serialize(val.neighbors.second, buffer);
  serialize(val.comp, buffer);
  serialize(val.ID, buffer);
  serialize(val.clips.size(), buffer);
  for (const auto x: val.clips) serialize(x, buffer);
}

//------------------------------------------------------------------------------
// Serialize a Vertex3d
//------------------------------------------------------------------------------
template<typename VA>
inline
void
serialize(const Vertex3d<VA>& val,
          std::vector<char>& buffer) {
  serialize<VA>(val.position, buffer);
  serialize(val.neighbors.size(), buffer);
  for (const auto x: val.neighbors) serialize(x, buffer);
  serialize(val.comp, buffer);
  serialize(val.ID, buffer);
  serialize(val.clips.size(), buffer);
  for (const auto x: val.clips) serialize(x, buffer);
}

//------------------------------------------------------------------------------
// Serialize a polygon
//------------------------------------------------------------------------------
template<typename VA>
inline
void
serialize(const std::vector<Vertex2d<VA>>& val,
          std::vector<char>& buffer) {
  serialize(val.size(), buffer);
  for (const auto& x: val) serialize(x, buffer);
}

//------------------------------------------------------------------------------
// Serialize a polyhedron
//------------------------------------------------------------------------------
template<typename VA>
inline
void
serialize(const std::vector<Vertex3d<VA>>& val,
          std::vector<char>& buffer) {
  serialize(val.size(), buffer);
  for (const auto& x: val) serialize(x, buffer);
}

//------------------------------------------------------------------------------
// Serialize a plane
//------------------------------------------------------------------------------
template<typename VA>
inline
void
serialize(const Plane<VA>& val,
          std::vector<char>& buffer) {
  serialize(val.dist, buffer);
  serialize<VA>(val.normal, buffer);
  serialize(val.ID, buffer);
}

//------------------------------------------------------------------------------
// Deserialize a double
//------------------------------------------------------------------------------
inline
void
deserialize(double& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  const auto n = sizeof(double);
  char* data = reinterpret_cast<char*>(&val);
  std::copy(itr, itr + n, data);
  itr += n;
  PCASSERT(itr <= endBuffer);
}

//------------------------------------------------------------------------------
// Deserialize an int
//------------------------------------------------------------------------------
inline
void
deserialize(int& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  const auto n = sizeof(int);
  char* data = reinterpret_cast<char*>(&val);
  std::copy(itr, itr + n, data);
  itr += n;
  PCASSERT(itr <= endBuffer);
}

//------------------------------------------------------------------------------
// Deserialize a size_t
//------------------------------------------------------------------------------
inline
void
deserialize(size_t& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  const auto n = sizeof(size_t);
  char* data = reinterpret_cast<char*>(&val);
  std::copy(itr, itr + n, data);
  itr += n;
  PCASSERT(itr <= endBuffer);
}

//------------------------------------------------------------------------------
// Deserialize a std::string
//------------------------------------------------------------------------------
inline
void
deserialize(std::string& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  size_t n;
  deserialize(n, itr, endBuffer);
  val.resize(n);
  std::copy(itr, itr+n, val.begin());
  itr += n;
  PCASSERT(itr <= endBuffer);
}

//------------------------------------------------------------------------------
// Deserialize a Vector
//------------------------------------------------------------------------------
template<typename VA>
inline
void
deserialize(typename VA::VECTOR& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  std::array<double, 3> triple_vals;
  deserialize(triple_vals[0],itr, endBuffer);
  deserialize(triple_vals[1],itr, endBuffer);
  deserialize(triple_vals[2],itr, endBuffer);
  VA::set_triple(val, triple_vals);
}

//------------------------------------------------------------------------------
// Deserialize a Vertex2d
//------------------------------------------------------------------------------
template<typename VA>
inline
void
deserialize(Vertex2d<VA>& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  size_t n;
  int x;
  deserialize<VA>(val.position, itr, endBuffer);
  deserialize(val.neighbors.first, itr, endBuffer);
  deserialize(val.neighbors.second, itr, endBuffer);
  deserialize(val.comp, itr, endBuffer);
  deserialize(val.ID, itr, endBuffer);
  deserialize(n, itr, endBuffer);
  val.clips.clear();
  for (auto i = 0; i < n; ++i) {
    deserialize(x, itr, endBuffer);
    val.clips.insert(x);
  }
}

//------------------------------------------------------------------------------
// Deserialize a Vertex3d
//------------------------------------------------------------------------------
template<typename VA>
inline
void
deserialize(Vertex3d<VA>& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  size_t n;
  int x;
  deserialize<VA>(val.position, itr, endBuffer);
  deserialize(n, itr, endBuffer);
  val.neighbors.resize(n);
  for (auto i = 0; i < n; ++i) deserialize(val.neighbors[i], itr, endBuffer);
  deserialize(val.comp, itr, endBuffer);
  deserialize(val.ID, itr, endBuffer);
  deserialize(n, itr, endBuffer);
  val.clips.clear();
  for (auto i = 0; i < n; ++i) {
    deserialize(x, itr, endBuffer);
    val.clips.insert(x);
  }
}

//------------------------------------------------------------------------------
// Deserialize a polygon
//------------------------------------------------------------------------------
template<typename VA>
inline
void
deserialize(std::vector<Vertex2d<VA>>& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  val.clear();
  size_t n;
  Vertex2d<VA> x;
  for (auto i = 0; i < n; ++i) {
    deserialize(x, itr, endBuffer);
    val.push_back(x);
  }
}

//------------------------------------------------------------------------------
// Deserialize a polyhedron
//------------------------------------------------------------------------------
template<typename VA>
inline
void
deserialize(std::vector<Vertex3d<VA>>& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  val.clear();
  size_t n;
  Vertex3d<VA> x;
  for (auto i = 0; i < n; ++i) {
    deserialize(x, itr, endBuffer);
    val.push_back(x);
  }
}

//------------------------------------------------------------------------------
// Deserialize a plane
//------------------------------------------------------------------------------
template<typename VA>
inline
void
deserialize(Plane<VA>& val,
            std::vector<char>::const_iterator& itr,
            const std::vector<char>::const_iterator& endBuffer) {
  deserialize(val.dist, itr, endBuffer);
  deserialize<VA>(val.normal, itr, endBuffer);
  deserialize(val.ID, itr, endBuffer);
}

}
}

#endif
