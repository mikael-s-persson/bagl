// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ADJACENCY_RANGE_H_
#define BAGL_BAGL_ADJACENCY_RANGE_H_

#include <ranges>

#include "bagl/graph_traits.h"

namespace bagl {

template <typename OutEdgeRange, typename Graph>
auto adjacency_range(OutEdgeRange&& r, const Graph& g) {
  return std::views::transform(r, [&g](const auto& e) { return target(e, g); });
}

template <typename InEdgeRange, typename Graph>
auto inv_adjacency_range(InEdgeRange&& r, const Graph& g) {
  return std::views::transform(r, [&g](const auto& e) { return source(e, g); });
}

}  // namespace bagl

#endif  // BAGL_BAGL_ADJACENCY_RANGE_H_
