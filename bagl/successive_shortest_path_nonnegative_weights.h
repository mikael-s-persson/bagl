// Copyright 2013 University of Warsaw.
// Authors: Piotr Wygocki
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SUCCESSIVE_SHORTEST_PATH_NONNEGATIVE_WEIGHTS_H_
#define BAGL_BAGL_SUCCESSIVE_SHORTEST_PATH_NONNEGATIVE_WEIGHTS_H_

#include <numeric>

#include "bagl/detail/augment.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

// This algorithm is described in "Network Flows: Theory, Algorithms, and
// Applications"
// by Ahuja, Magnanti, Orlin.

namespace bagl {

namespace successive_shortest_path_detail {

template <concepts::Graph G, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::ReadableVertexPropertyMap<G> Distance>
class map_reduced_weight : public put_get_helper<map_reduced_weight<G, Weight, Distance> > {
 public:
  using value_type = property_traits_value_t<Weight>;
  using reference = value_type;
  using key_type = graph_edge_descriptor_t<G>;

  map_reduced_weight(const G& g, Weight w, Distance d) : g_(g), weight_(w), distance_(d) {}

  reference operator[](key_type e) const {
    return get(distance_, source(e, g_)) - get(distance_, target(e, g_)) + get(weight_, e);
  }

 private:
  const G& g_;
  Weight weight_;
  Distance distance_;
};

}  // namespace successive_shortest_path_detail

// See overloads below for simpler calls and sensible defaults.

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Capacity,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacity, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::ReadableEdgePropertyMap<G> Reversed, concepts::ReadableVertexIndexMap<G> VertexIndex,
          concepts::ReadWriteVertexPropertyMap<G> Pred, concepts::ReadWriteVertexPropertyMap<G> Distance,
          concepts::ReadWriteVertexPropertyMap<G> Distance2>
void successive_shortest_path_nonnegative_weights(const G& g, graph_vertex_descriptor_t<G> s,
                                                  graph_vertex_descriptor_t<G> t,
                                                  Capacity capacity,                   // get(edge_capacity, g)
                                                  ResidualCapacity residual_capacity,  // get(edge_residual_capacity, g)
                                                  Weight weight,                       // get(edge_weight, g)
                                                  Reversed rev,                        // get(edge_reverse, g)
                                                  VertexIndex index,                   // get(vertex_index, g)
                                                  Pred pred, Distance distance, Distance2 distance_prev) {
  auto gres = augment_detail::residual_graph(g, residual_capacity);
  using edge_descriptor = graph_edge_descriptor_t<G>;

  for (auto e : edges(g)) {
    put(residual_capacity, e, get(capacity, e));
  }

  for (auto v : vertices(g)) {
    put(distance_prev, v, 0);
  }

  while (true) {
    for (auto v : vertices(g)) {
      put(pred, v, edge_descriptor());
    }

    dijkstra_shortest_paths(gres, s, vector_property_map(num_vertices(g), index, graph_traits<G>::null_vertex()),
                            distance, successive_shortest_path_detail::map_reduced_weight(gres, weight, distance_prev),
                            index, make_dijkstra_visitor(edge_predecessor_recorder_on_edge_relaxed(pred)));

    if (get(pred, t) == edge_descriptor()) {
      break;
    }

    for (auto v : vertices(g)) {
      put(distance_prev, v, get(distance_prev, v) + get(distance, v));
    }

    augment_detail::augment(g, s, t, pred, residual_capacity, rev);
  }
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Capacity,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacity, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::ReadableEdgePropertyMap<G> Reversed, concepts::ReadableVertexIndexMap<G> VertexIndex,
          concepts::ReadWriteVertexPropertyMap<G> Pred, concepts::ReadWriteVertexPropertyMap<G> Distance>
void successive_shortest_path_nonnegative_weights(const G& g, graph_vertex_descriptor_t<G> s,
                                                  graph_vertex_descriptor_t<G> t, Capacity capacity,
                                                  ResidualCapacity residual_capacity, Weight weight, Reversed rev,
                                                  VertexIndex index, Pred pred, Distance distance) {
  using D = property_traits_value_t<Weight>;
  successive_shortest_path_nonnegative_weights(g, s, t, capacity, residual_capacity, weight, rev, index, pred, distance,
                                               vector_property_map(num_vertices(g), index, D{}));
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Capacity,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacity, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::ReadableEdgePropertyMap<G> Reversed, concepts::ReadableVertexIndexMap<G> VertexIndex,
          concepts::ReadWriteVertexPropertyMap<G> Pred>
void successive_shortest_path_nonnegative_weights(const G& g, graph_vertex_descriptor_t<G> s,
                                                  graph_vertex_descriptor_t<G> t, Capacity capacity,
                                                  ResidualCapacity residual_capacity, Weight weight, Reversed rev,
                                                  VertexIndex index, Pred pred) {
  using D = property_traits_value_t<Weight>;
  successive_shortest_path_nonnegative_weights(g, s, t, capacity, residual_capacity, weight, rev, index, pred,
                                               vector_property_map(num_vertices(g), index, D{}));
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Capacity,
          concepts::ReadWriteEdgePropertyMap<G> ResidualCapacity, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::ReadableEdgePropertyMap<G> Reversed, concepts::ReadableVertexIndexMap<G> VertexIndex>
void successive_shortest_path_nonnegative_weights(const G& g, graph_vertex_descriptor_t<G> s,
                                                  graph_vertex_descriptor_t<G> t, Capacity capacity,
                                                  ResidualCapacity residual_capacity, Weight weight, Reversed rev,
                                                  VertexIndex index) {
  successive_shortest_path_nonnegative_weights(
      g, s, t, capacity, residual_capacity, weight, rev, index,
      vector_property_map(num_vertices(g), index, graph_edge_descriptor_t<G>{}));
}

template <concepts::VertexAndEdgeListGraph G>
void successive_shortest_path_nonnegative_weights(G& g, graph_vertex_descriptor_t<G> s,
                                                  graph_vertex_descriptor_t<G> t) {
  successive_shortest_path_nonnegative_weights(g, s, t, get(edge_capacity, g), get(edge_residual_capacity, g),
                                               get(edge_weight, g), get(edge_reverse, g), get(vertex_index, g));
}

}  // namespace bagl
#endif  // BAGL_BAGL_SUCCESSIVE_SHORTEST_PATH_NONNEGATIVE_WEIGHTS_H_
