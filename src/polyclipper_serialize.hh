//---------------------------------PolyClipper--------------------------------//
// Methods for serializing/dserializing PolyClipper types
//----------------------------------------------------------------------------//
#ifndef __PolyClipper_serialize__
#define __PolyClipper_serialize__

#include "polyclipper_utilities.hh"

namespace PolyClipper {

// Forward declarations
template<typename VA> class Vertex2d;
template<typename VA> class Vertex3d;

namespace internal {

// Serialize
void serialize(const double val, std::string& buffer);
void serialize(const int val, std::string& buffer);
void serialize(const size_t val, std::string& buffer);
void serialize(const std::string& val, std::string& buffer);
template<typename VA> void serialize(const typename VA::VECTOR& val, std::string& buffer);
template<typename VA> void serialize(const Vertex2d<VA>& val, std::string& buffer);
template<typename VA> void serialize(const Vertex3d<VA>& val, std::string& buffer);
template<typename VA> void serialize(const std::vector<Vertex2d<VA>>& val, std::string& buffer);
template<typename VA> void serialize(const std::vector<Vertex3d<VA>>& val, std::string& buffer);
template<typename VA> void serialize(const Plane<VA>& val, std::string& buffer);
template<typename VA> void serialize(const std::vector<Plane<VA>>& val, std::string& buffer);

// Deserialize
void deserialize(double& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
void deserialize(int& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
void deserialize(size_t& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
void deserialize(std::string& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(typename VA::VECTOR& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(Vertex2d<VA>& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(Vertex3d<VA>& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(std::vector<Vertex2d<VA>>& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(std::vector<Vertex3d<VA>>& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(Plane<VA>& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);
template<typename VA> void deserialize(std::vector<Plane<VA>>& val, std::string::const_iterator& itr, const std::string::const_iterator& endBuffer);

}
}

#include "polyclipper_serializeImpl.hh"

#endif
