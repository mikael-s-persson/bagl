// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BELLMAN_FORD_SHORTEST_PATHS_H_
#define BAGL_BAGL_BELLMAN_FORD_SHORTEST_PATHS_H_

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/relax.h"
#include "bagl/visitors.h"

namespace bagl {

//  This file implements the function
//
//  template <class EdgeListGraph, class Size, class P, class T, class R>
//  bool bellman_ford_shortest_paths(EdgeListGraph& g, Size N,
//     const bgl_named_params<P, T, R>& params)

namespace concepts {

template <typename V, typename G>
concept BellmanFordVisitor = std::copy_constructible<V> &&
    requires(const V& vis, const G& g, graph_edge_descriptor_t<G> e) {
  vis.examine_edge(e, g);
  vis.edge_relaxed(e, g);
  vis.edge_not_relaxed(e, g);
  vis.edge_minimized(e, g);
  vis.edge_not_minimized(e, g);
};

}  // namespace concepts

template <typename Visitors = null_visitors>
class bellman_visitor {
 public:
  bellman_visitor() = default;
  explicit bellman_visitor(Visitors vis) : vis_(vis) {}

  template <typename Edge, typename Graph>
  void examine_edge(Edge e, Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void edge_relaxed(Edge e, Graph& g) {
    visitors_detail::invoke_edge_relaxed_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void edge_not_relaxed(Edge e, Graph& g) {
    visitors_detail::invoke_edge_not_relaxed_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void edge_minimized(Edge e, Graph& g) {
    visitors_detail::invoke_edge_minimized_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void edge_not_minimized(Edge e, Graph& g) {
    visitors_detail::invoke_edge_not_minimized_on_all(vis_, e, g);
  }

 protected:
  Visitors vis_{};
};

template <typename... Visitors>
auto make_bellman_visitor(Visitors&&... vis) {
  if constexpr (sizeof...(Visitors) == 0) {
    return bellman_visitor<>();
  } else {
    return bellman_visitor<std::tuple<std::decay_t<Visitors>...>>(std::forward<Visitors>(vis)...);
  }
}
using default_bellman_visitor = bellman_visitor<>;

template <concepts::EdgeListGraph G, concepts::BellmanFordVisitor<G> V,
          concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, class BinaryFunction,
          class BinaryPredicate>
bool bellman_ford_shortest_paths(G& g, V vis, std::size_t n, WeightMap weight, PredecessorMap pred, DistanceMap distance,
                                 BinaryFunction combine, BinaryPredicate compare) {
  for (std::size_t k = 0; k < n; ++k) {
    bool at_least_one_edge_relaxed = false;
    for (auto e : edges(g)) {
      vis.examine_edge(e, g);
      if (relax(e, g, weight, pred, distance, combine, compare)) {
        at_least_one_edge_relaxed = true;
        vis.edge_relaxed(e, g);
      } else {
        vis.edge_not_relaxed(e, g);
      }
    }
    if (!at_least_one_edge_relaxed) {
      break;
    }
  }

  for (auto e : edges(g)) {
    if (compare(combine(get(distance, source(e, g)), get(weight, e)), get(distance, target(e, g)))) {
      vis.edge_not_minimized(e, g);
      return false;
    }
    vis.edge_minimized(e, g);
  }

  return true;
}

template <concepts::EdgeListGraph G, concepts::BellmanFordVisitor<G> V,
          concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap>
bool bellman_ford_shortest_paths(G& g, V vis, std::size_t n, WeightMap weight, PredecessorMap pred, DistanceMap distance) {
  return bellman_ford_shortest_paths(g, vis, n, weight, pred, distance, std::plus<>(), std::less<>());
}

template <concepts::VertexAndEdgeListGraph G, concepts::BellmanFordVisitor<G> V,
          concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap, class BinaryFunction,
          class BinaryPredicate>
bool bellman_ford_shortest_paths(G& g, V vis, graph_vertex_descriptor_t<G> start, std::size_t n, WeightMap weight,
                                 PredecessorMap pred, DistanceMap distance, BinaryFunction combine,
                                 BinaryPredicate compare) {
  using weight_type = property_traits_value_t<WeightMap>;
  for (auto v : vertices(g)) {
    put(distance, v, (std::numeric_limits<weight_type>::max)());
    put(pred, v, v);
  }
  put(distance, start, weight_type(0));
  return bellman_ford_shortest_paths(g, vis, n, weight, pred, distance, combine, compare);
}

template <concepts::VertexAndEdgeListGraph G, concepts::BellmanFordVisitor<G> V,
          concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadWriteVertexPropertyMap<G> PredecessorMap,
          concepts::ReadWriteVertexPropertyMap<G> DistanceMap>
bool bellman_ford_shortest_paths(G& g, V vis, graph_vertex_descriptor_t<G> start, std::size_t n, WeightMap weight,
                                 PredecessorMap pred, DistanceMap distance) {
  return bellman_ford_shortest_paths(g, vis, start, n, weight, pred, distance, std::plus<>(), std::less<>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_BELLMAN_FORD_SHORTEST_PATHS_H_
