// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

// This algorithm is described in "Network Flows: Theory, Algorithms, and
// Applications"
// by Ahuja, Magnanti, Orlin.

#ifndef BAGL_BAGL_CYCLE_CANCELING_H_
#define BAGL_BAGL_CYCLE_CANCELING_H_

#include <numeric>

#include "bagl/bellman_ford_shortest_paths.h"
#include "bagl/detail/augment.h"
#include "bagl/find_flow_cost.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace cycle_detail {

template <typename PredEdgeMap, typename Vertex>
class record_edge_map_and_cycle_vertex : public predecessor_recorder_on_tree_edge<PredEdgeMap> {
  using Base = predecessor_recorder_on_tree_edge<PredEdgeMap>;

 public:
  record_edge_map_and_cycle_vertex(PredEdgeMap pred, Vertex& v) : Base(pred), v_(&v) {}

  template <typename Graph, typename Edge>
  void edge_not_minimized(Edge e, const Graph& g) const {
    std::size_t n = num_vertices(g) + 1;

    // edge e is not minimized but does not have to be on the negative
    // weight cycle to find vertex on negative wieight cycle we move n+1
    // times backword in the PredEdgeMap graph.
    while (n > 0) {
      e = get(Base::predecessor, source(e, g));
      --n;
    }
    *v_ = source(e, g);
  }

 private:
  Vertex* v_;
};

}  // namespace cycle_detail

template <concepts::VertexListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> Weight,
          concepts::ReadableEdgePropertyMap<Graph> Reversed, concepts::ReadWriteEdgePropertyMap<Graph> ResidualCapacity,
          concepts::ReadWriteVertexPropertyMap<Graph> Pred, concepts::ReadWriteVertexPropertyMap<Graph> Distance>
requires std::same_as<graph_edge_descriptor_t<Graph>, property_traits_value_t<Pred>> &&
    std::convertible_to<property_traits_value_t<Weight>, property_traits_value_t<Distance>>
void cycle_canceling(const Graph& g,
                     Weight weight,                       // get(edge_weight, g)
                     Reversed rev,                        // get(edge_reverse, g)
                     ResidualCapacity residual_capacity,  // get(edge_residual_capacity, g)
                     Pred pred, Distance distance) {
  auto gres = augment_detail::residual_graph(g, residual_capacity);

  using ResGraph = decltype(gres);
  using Edge = graph_edge_descriptor_t<ResGraph>;
  using Vertex = graph_vertex_descriptor_t<ResGraph>;

  for (auto v : vertices(g)) {
    put(pred, v, Edge());
    put(distance, v, 0);
  }

  std::size_t n = num_vertices(g);
  Vertex cycle_start;

  while (!bellman_ford_shortest_paths(
      gres, make_bellman_visitor(cycle_detail::record_edge_map_and_cycle_vertex<Pred, Vertex>(pred, cycle_start)), n,
      weight, get(vertex_predecessor, gres), distance)) {
    augment_detail::augment(g, cycle_start, cycle_start, pred, residual_capacity, rev);

    for (auto v : vertices(g)) {
      put(pred, v, Edge());
      put(distance, v, 0);
    }
  }
}

template <concepts::VertexListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> Weight,
          concepts::ReadableEdgePropertyMap<Graph> Reversed, concepts::ReadWriteEdgePropertyMap<Graph> ResidualCapacity,
          concepts::ReadableVertexPropertyMap<Graph> Pred>
requires std::same_as<graph_edge_descriptor_t<Graph>, property_traits_value_t<Pred>>
void cycle_canceling(const Graph& g, Weight weight, Reversed rev, ResidualCapacity residual_capacity, Pred pred) {
  cycle_canceling(g, weight, rev, residual_capacity, pred,
                  vector_property_map(num_vertices(g), get(vertex_index, g), property_traits_value_t<Weight>{}).ref());
}

template <concepts::VertexListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> Weight,
          concepts::ReadableEdgePropertyMap<Graph> Reversed, concepts::ReadWriteEdgePropertyMap<Graph> ResidualCapacity,
          concepts::ReadableVertexPropertyMap<Graph> Distance>
requires std::convertible_to<property_traits_value_t<Weight>, property_traits_value_t<Distance>>
void cycle_canceling(const Graph& g, Weight weight, Reversed rev, ResidualCapacity residual_capacity,
                     Distance distance) {
  cycle_canceling(g, weight, rev, residual_capacity,
                  vector_property_map(num_vertices(g), get(vertex_index, g), graph_edge_descriptor_t<Graph>{}).ref(),
                  distance);
}

template <concepts::VertexListGraph Graph, concepts::ReadableEdgePropertyMap<Graph> Weight,
          concepts::ReadableEdgePropertyMap<Graph> Reversed, concepts::ReadWriteEdgePropertyMap<Graph> ResidualCapacity>
void cycle_canceling(const Graph& g, Weight weight, Reversed rev, ResidualCapacity residual_capacity) {
  cycle_canceling(g, weight, rev, residual_capacity,
                  vector_property_map(num_vertices(g), get(vertex_index, g), graph_edge_descriptor_t<Graph>{}).ref(),
                  vector_property_map(num_vertices(g), get(vertex_index, g), property_traits_value_t<Weight>{}).ref());
}

}  // namespace bagl

#endif  // BAGL_BAGL_CYCLE_CANCELING_H_
