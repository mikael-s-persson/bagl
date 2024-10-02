// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_JOHNSON_ALL_PAIRS_SHORTEST_H_
#define BAGL_BAGL_JOHNSON_ALL_PAIRS_SHORTEST_H_

#include <type_traits>

#include "bagl/adjacency_list.h"
#include "bagl/bellman_ford_shortest_paths.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_traits.h"
#include "bagl/numeric_values.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::BasicMatrix<property_traits_value_t<Weight>> DistanceMatrix,
          concepts::ReadableVertexIndexMap<G> VertexID, concepts::PropertyComparator<Weight> Compare,
          concepts::PropertyCombinator<Weight> Combine>
bool johnson_all_pairs_shortest_paths(const G& g1, DistanceMatrix& D, VertexID id1, Weight w1, Compare compare,
                                      Combine combine, property_traits_value_t<Weight> inf,
                                      property_traits_value_t<Weight> zero) {
  using DT = property_traits_value_t<Weight>;

  constexpr bool is_undirected = is_undirected_graph_v<G>;

  using Graph2 = adjacency_list<vec_s, vec_s, directed_s, property<vertex_distance_t, DT>,
                                property<edge_weight_t, DT, property<edge_weight2_t, DT>>>;

  Graph2 g2(num_vertices(g1) + 1);
  auto w = get(edge_weight, g2);
  auto w_hat = get(edge_weight2, g2);
  auto d = get(vertex_distance, g2);
  auto id2 = get(vertex_index, g2);

  // Construct g2 where V[g2] = V[g1] U {s}
  //   and  E[g2] = E[g1] U {(s,v)| v in V[g1]}
  std::vector<graph_vertex_descriptor_t<G>> verts1(num_vertices(g1) + 1);
  auto s = *vertices(g2).begin();
  {
    int i = 1;
    for (auto v : vertices(g1)) {
      auto [e, z] = add_edge(s, get(id1, v) + 1, g2);
      put(w, e, zero);
      verts1[i] = v;
      ++i;
    }
    for (auto e : edges(g1)) {
      auto [e2, z] = add_edge(get(id1, source(e, g1)) + 1, get(id1, target(e, g1)) + 1, g2);
      put(w, e2, get(w1, e));
      if constexpr (is_undirected) {
        std::tie(e2, z) = add_edge(get(id1, target(e, g1)) + 1, get(id1, source(e, g1)) + 1, g2);
        put(w, e2, get(w1, e));
      }
    }
  }

  auto h = vector_property_map(num_vertices(g2), id2, DT{});

  for (auto v : vertices(g2)) {
    put(d, v, inf);
  }

  put(d, s, zero);
  // Using the non-named parameter versions of bellman_ford and
  // dijkstra for portability reasons.
  dummy_property_map pred;
  bellman_visitor<> bvis;
  if (bellman_ford_shortest_paths(g2, num_vertices(g2), w, pred, d, combine, compare, bvis)) {
    for (auto v : vertices(g2)) {
      put(h, v, get(d, v));
    }
    // Reweight the edges to remove negatives
    for (auto e : edges(g2)) {
      auto a = source(e, g2);
      auto b = target(e, g2);
      put(w_hat, e, combine((get(h, a) - get(h, b)), get(w, e)));
    }
    for (auto u : vertices(g2)) {
      dijkstra_visitor<> dvis;
      dijkstra_shortest_paths(g2, u, pred, d, w_hat, id2, compare, combine, inf, zero, dvis);
      for (auto v : vertices(g2)) {
        if (u != s && v != s) {
          D[get(id2, u) - 1][get(id2, v) - 1] = combine((get(h, v) - get(h, u)), get(d, v));
        }
      }
    }
    return true;
  }
  return false;
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::BasicMatrix<property_traits_value_t<Weight>> DistanceMatrix,
          concepts::ReadableVertexIndexMap<G> VertexID>
bool johnson_all_pairs_shortest_paths(const G& g1, DistanceMatrix& D, VertexID id1, Weight w1) {
  using DT = property_traits_value_t<Weight>;
  return johnson_all_pairs_shortest_paths(g1, D, id1, w1, std::less<>(), closed_plus<>(),
                                          numeric_values<DT>::infinity(), numeric_values<DT>::zero());
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableEdgePropertyMap<G> Weight,
          concepts::BasicMatrix<property_traits_value_t<Weight>> DistanceMatrix>
bool johnson_all_pairs_shortest_paths(const G& g1, DistanceMatrix& D, Weight w1) {
  using DT = property_traits_value_t<Weight>;
  return johnson_all_pairs_shortest_paths(g1, D, get(vertex_index, g1), w1, std::less<>(), closed_plus<>(),
                                          numeric_values<DT>::infinity(), numeric_values<DT>::zero());
}

template <concepts::VertexAndEdgeListGraph G, typename DistanceMatrix>
bool johnson_all_pairs_shortest_paths(const G& g1, DistanceMatrix& D) {
  return johnson_all_pairs_shortest_paths(g1, D, get(edge_weight, g1));
}

}  // namespace bagl

#endif  // BAGL_BAGL_JOHNSON_ALL_PAIRS_SHORTEST_H_
