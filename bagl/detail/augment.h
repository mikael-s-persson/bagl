// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_AUGMENT_H_
#define BAGL_BAGL_DETAIL_AUGMENT_H_

#include <limits>
#include <utility>

#include "bagl/filtered_graph.h"
#include "bagl/graph_traits.h"

namespace bagl::augment_detail {

template <concepts::Graph Graph, concepts::ReadableEdgePropertyMap<Graph> ResCapMap>
auto residual_graph(const Graph& g, ResCapMap residual_capacity) {
  return filtered_graph<const Graph, is_residual_edge<ResCapMap> >(g, is_residual_edge<ResCapMap>(residual_capacity));
}

template <concepts::Graph Graph, concepts::ReadableVertexPropertyMap<Graph> PredEdgeMap,
          concepts::ReadWriteEdgePropertyMap<Graph> ResCapMap, concepts::ReadableEdgePropertyMap<Graph> RevEdgeMap>
void augment(const Graph& g, graph_vertex_descriptor_t<Graph> src, graph_vertex_descriptor_t<Graph> sink, PredEdgeMap p,
             ResCapMap residual_capacity, RevEdgeMap reverse_edge) {
  using FlowValue = property_traits_value_t<ResCapMap>;

  // find minimum residual capacity along the augmenting path
  auto u = graph_traits<Graph>::null_vertex();
  FlowValue delta = std::numeric_limits<FlowValue>::max();
  graph_edge_descriptor_t<Graph> e = get(p, sink);
  do {
    delta = std::min(delta, get(residual_capacity, e));
    u = source(e, g);
    e = get(p, u);
  } while (u != src);

  // push delta units of flow along the augmenting path
  e = get(p, sink);
  do {
    put(residual_capacity, e, get(residual_capacity, e) - delta);
    put(residual_capacity, get(reverse_edge, e), get(residual_capacity, get(reverse_edge, e)) + delta);
    u = source(e, g);
    e = get(p, u);
  } while (u != src);
}

}  // namespace bagl::augment_detail

#endif  // BAGL_BAGL_DETAIL_AUGMENT_H_
