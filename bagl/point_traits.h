// Copyright 2004, 2005 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_POINT_TRAITS_H_
#define BAGL_BAGL_POINT_TRAITS_H_

#include <cstddef>

namespace bagl {

template <typename Point>
struct point_traits {
  // The type of each component of the point
  using component_type = typename Point::component_type;

  // The number of dimensions in the point
  static constexpr std::size_t dimensions = Point::dimensions;
};

}  // namespace bagl

#endif  // BAGL_BAGL_POINT_TRAITS_H_
