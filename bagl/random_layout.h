// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_RANDOM_LAYOUT_H_
#define BAGL_BAGL_RANDOM_LAYOUT_H_

#include <random>
#include <type_traits>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"

namespace bagl {

template <typename Topology, concepts::VertexListGraph Graph, concepts::WritableVertexPropertyMap<Graph> PositionMap>
void random_graph_layout(const Graph& g, PositionMap position_map, const Topology& topology) {
  for (auto v : vertices(g)) {
    put(position_map, v, topology.random_point());
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_RANDOM_LAYOUT_H_
