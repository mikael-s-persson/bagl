// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_NUMERIC_VALUES_H_
#define BAGL_BAGL_NUMERIC_VALUES_H_

#include <concepts>
#include <limits>

namespace bagl {

// This generic type reports various numeric values for some type. In the
// general case, numeric values simply treat their maximum value as infinity
// and the default-constructed value as 0.
//
// Specializations of this template can redefine the notions of zero and
// infinity for various types. For example, the class is specialized for
// floating point types to use the built in notion of infinity.
template <typename T>
struct numeric_values {
  using value_type = T;
  static T zero() { return T{}; }
  static T infinity() { return (std::numeric_limits<T>::max)(); }
};

// Specializations for floating point types refer to 0.0 and their infinity
// value defined by numeric_limits.
template <typename T>
struct numeric_values_floating_point {
  using value_type = T;
  static T zero() { return T{0.0}; }
  static T infinity() { return std::numeric_limits<T>::infinity(); }
};
template <>
struct numeric_values<float> : numeric_values_floating_point<float> {};
template <>
struct numeric_values<double> : numeric_values_floating_point<double> {};
template <>
struct numeric_values<long double> : numeric_values_floating_point<long double> {};

}  // namespace bagl

#endif  // BAGL_BAGL_NUMERIC_VALUES_H_
