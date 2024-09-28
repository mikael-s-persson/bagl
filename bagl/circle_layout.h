// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CIRCLE_LAYOUT_H_
#define BAGL_BAGL_CIRCLE_LAYOUT_H_

#include <cmath>
#include <numbers>
#include <utility>

#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/topology.h"

namespace bagl {

// Layout the graph with the vertices at the points of a regular n-polygon.
//
// The distance from the center of the polygon to each point is
// determined by the `radius` parameter. The `position` parameter
// must be an Lvalue Property Map whose value type is a class type
/// containing `x` and `y` members that will be set to the `x` and
// `y` coordinates.
template <typename VertexListGraph, typename PositionMap, typename Radius>
void circle_graph_layout(const VertexListGraph& g, PositionMap position, Radius radius) {
  static_assert(property_traits_value_t<PositionMap>::dimensions >= 2);

  const std::size_t n = num_vertices(g);
  const double two_pi_over_n = 2.0 * std::numbers::pi / n;
  std::size_t i = 0;
  for (auto v : vertices(g)) {
    position[v][0] = radius * std::cos(static_cast<double>(i) * two_pi_over_n);
    position[v][1] = radius * std::sin(static_cast<double>(i) * two_pi_over_n);
    ++i;
  }
}
}  // namespace bagl

#endif  // BAGL_BAGL_CIRCLE_LAYOUT_H_
