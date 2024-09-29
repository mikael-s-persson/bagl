// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_FIND_FLOW_COST_H_
#define BAGL_BAGL_FIND_FLOW_COST_H_

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

template <concepts::EdgeListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> Capacity,
          concepts::ReadableEdgePropertyMap<Graph> ResidualCapacity, concepts::ReadableEdgePropertyMap<Graph> Weight>
auto find_flow_cost(const Graph& g, Capacity capacity, ResidualCapacity residual_capacity, Weight weight) {
  using Cost = property_traits_value_t<Weight>;

  Cost cost = 0;
  for (auto e : edges(g)) {
    if (get(capacity, e) > Cost(0)) {
      cost += (get(capacity, e) - get(residual_capacity, e)) * get(weight, e);
    }
  }
  return cost;
}

template <class Graph>
auto find_flow_cost(const Graph& g) {
  return find_flow_cost(g, get(edge_capacity, g), get(edge_residual_capacity, g), get(edge_weight, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_FIND_FLOW_COST_H_
