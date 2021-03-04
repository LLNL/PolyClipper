//------------------------------------------------------------------------------
// A collection of internal utilities common in PolyClipper.
//
// These largely repoduce things in Spheral to make PolyClipper standalone.
//------------------------------------------------------------------------------
#ifndef __polyclipper_utilities__
#define __polyclipper_utilities__

#include "polyclipper_vector2d.hh"
#include "polyclipper_vector3d.hh"
#include "polyclipper_plane.hh"

#include <iostream>
#include <ostream>
#include <vector>
#include <stdexcept>

//------------------------------------------------------------------------------
// Define an assert command with optional message
//------------------------------------------------------------------------------
#ifndef NDEBUG
#   define PCASSERT2(condition, message) \
    do { \
        if (! (condition)) { \
            std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                      << " line " << __LINE__ << ": " << message << std::endl; \
            throw std::runtime_error("PolyClipper ERROR");                                              \
        } \
    } while (false)
#else
   define PCASSERT2(condition, message) do { } while (false)
#endif

#define PCASSERT(condition) PCASSERT2(condition, #condition)

namespace PolyClipper {
namespace internal {

//------------------------------------------------------------------------------
// Return the sign of the argument determined as follows:
//   
//    x > 0 -> sgn0(x) =  1
//    x = 0 -> sgn0(x) =  0
//    x < 0 -> sgn0(x) = -1
//------------------------------------------------------------------------------
inline
double
sgn0(const double x) {
  return (x > 0.0 ?  1.0 :
          x < 0.0 ? -1.0 :
          0.0);
}

//------------------------------------------------------------------------------
// Return the sign of the argument determined as follows:
//   
//    x >= 0 -> sgn(x) =  1
//    x <  0 -> sgn(x) = -1
//------------------------------------------------------------------------------
inline
double
sgn(const double x) {
  return (x >= 0.0 ?  1.0 : -1.0);
}

//------------------------------------------------------------------------------
// Fuzzy comparisons.
//------------------------------------------------------------------------------
template<typename DataType>
inline
bool fuzzyEqual(const DataType& lhs, const DataType& rhs,
                const double fuzz = 1.0e-15) {
  return std::abs(lhs - rhs)/std::max(1.0, std::abs(lhs) + std::abs(rhs)) < fuzz;
}

template<typename DataType>
inline
bool fuzzyLessThanOrEqual(const DataType& lhs, const DataType& rhs,
                          const double fuzz = 1.0e-15) {
  return lhs < rhs || fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
inline
bool fuzzyGreaterThanOrEqual(const DataType& lhs, const DataType& rhs,
                             const double fuzz = 1.0e-15) {
  return lhs > rhs || fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
inline
bool distinctlyLessThan(const DataType& lhs, const DataType& rhs,
                        const double fuzz = 1.0e-15) {
  return lhs < rhs && !fuzzyEqual(lhs, rhs, fuzz);
}

template<typename DataType>
inline
bool distinctlyGreaterThan(const DataType& lhs, const DataType& rhs,
                           const double fuzz = 1.0e-15) {
  return lhs > rhs && !fuzzyEqual(lhs, rhs, fuzz);
}

//------------------------------------------------------------------------------
// safeInv
//
// Take the inverse of a number, turning it to zero if the argument is 0.0.
//------------------------------------------------------------------------------
template<typename Value>
inline
Value
safeInv(const Value& x,
        const double fuzz = 1.0e-30) {
  return sgn(x)/std::max(fuzz, std::abs(x));
}

//------------------------------------------------------------------------------
// removeElements
//
// Removes the specified elements by index from a std::vector.
// This is needed because there doesn't seem to be a good solution for this in
// the STL.  The std::vector::erase method involves N^2 operations when you're
// removing many elements, while the std::remove and std::remove_if do not work
// for removing elements by index/iterator.
//------------------------------------------------------------------------------
template<typename Value, typename index_t>
inline
void
removeElements(std::vector<Value>& vec,
	       const std::vector<index_t>& elements) {

  // Is there anything to do?
  if (elements.size() > 0) {

    const index_t originalSize = vec.size();
    const index_t newSize = originalSize - elements.size();

    // Remove the elements.
    // We prefer not to use the vector::erase here 'cause if we're removing
    // many elements the copy and move behaviour of erase can make this
    // an N^2 thing.  Yuck!
    typename std::vector<index_t>::const_iterator delItr = elements.begin();
    typename std::vector<index_t>::const_iterator endItr = elements.end();
    index_t i = *delItr;
    index_t j = i + 1;
    ++delItr;
    while (j != originalSize and delItr != endItr) {
      if (j == *delItr) {
        ++delItr;
        ++j;
      } else {
        vec[i] = vec[j];
        ++i;
        ++j;
      }
    }
    if (j != originalSize) std::copy(vec.begin() + j, vec.end(), vec.begin() + i);

    // Resize vec to it's new size.
    vec.erase(vec.begin() + newSize, vec.end());
  }
}

//------------------------------------------------------------------------------
// Compare a plane and point.
//------------------------------------------------------------------------------
template<typename VA>
inline
int compare(const Plane<VA>& plane,
            const typename VA::VECTOR& point) {
  const auto sgndist = plane.dist + VA::dot(plane.normal, point);
  if (std::abs(sgndist) < 1.0e-10) return 0;
  return sgn0(sgndist);
}

//------------------------------------------------------------------------------
// Intersect a line-segment with a plane.
//------------------------------------------------------------------------------
template<typename VA>
inline
typename VA::VECTOR
segmentPlaneIntersection(const typename VA::VECTOR& a,         // line-segment begin
                         const typename VA::VECTOR& b,         // line-segment end
                         const Plane<VA>& plane) { // plane
  const auto asgndist = plane.dist + VA::dot(plane.normal, a);
  const auto bsgndist = plane.dist + VA::dot(plane.normal, b);
  PCASSERT(asgndist != bsgndist);
  return VA::div(VA::sub(VA::mul(a, bsgndist), VA::mul(b, asgndist)), bsgndist - asgndist);
}

}
}

//------------------------------------------------------------------------------
// Vector2d ostream
//------------------------------------------------------------------------------
inline
std::ostream&
operator<<(std::ostream& os, const PolyClipper::Vector2d& vec) {
  os << "( " << vec.x << " " << vec.y << ")";
  return os;
}

//------------------------------------------------------------------------------
// double * Vector2d
//------------------------------------------------------------------------------
inline
PolyClipper::Vector2d
operator*(const double lhs, const PolyClipper::Vector2d& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Vector3d ostream
//------------------------------------------------------------------------------
inline
std::ostream&
operator<<(std::ostream& os, const PolyClipper::Vector3d& vec) {
  os << "( " << vec.x << " " << vec.y << " " << vec.z << ")";
  return os;
}

//------------------------------------------------------------------------------
// double * Vector3d
//------------------------------------------------------------------------------
inline
PolyClipper::Vector3d
operator*(const double lhs, const PolyClipper::Vector3d& rhs) {
  return rhs*lhs;
}

//------------------------------------------------------------------------------
// Plane ostream
//------------------------------------------------------------------------------
template<typename VA>
inline
std::ostream&
operator<<(std::ostream& os, const PolyClipper::Plane<VA>& plane) {
  os << "Plane[ " << plane.dist << " " << VA::str(plane.normal) << " " << plane.ID << "]";
  return os;
}

#endif
