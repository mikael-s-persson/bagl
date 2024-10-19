// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BETWEENNESS_CENTRALITY_H_
#define BAGL_BAGL_BETWEENNESS_CENTRALITY_H_

#include <algorithm>
#include <functional>
#include <limits>
#include <queue>
#include <stack>
#include <tuple>
#include <type_traits>
#include <vector>

#include "bagl/breadth_first_search.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/relax.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace bc_detail {

// Customized visitor passed to Dijkstra's algorithm by Brandes'
// betweenness centrality algorithm. This visitor is responsible for
// keeping track of the order in which vertices are discovered, the
// predecessors on the shortest path(s) to a vertex, and the number
// of shortest paths.
template <concepts::Graph G, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::MutableLvalueVertexPropertyMap<G> IncomingMap, concepts::ReadableVertexPropertyMap<G> DistanceMap,
          concepts::ReadWriteVertexPropertyMap<G> PathCountMap>
struct brandes_dijkstra_visitor : public bfs_visitor<> {
  using vertex_descriptor = graph_vertex_descriptor_t<G>;
  using edge_descriptor = graph_edge_descriptor_t<G>;

  brandes_dijkstra_visitor(std::stack<vertex_descriptor>& ordered_vertices, WeightMap weight, IncomingMap incoming,
                           DistanceMap distance, PathCountMap path_count)
      : ordered_vertices_(ordered_vertices),
        weight_(weight),
        incoming_(incoming),
        distance_(distance),
        path_count_(path_count) {}

  // Whenever an edge e = (v, w) is relaxed, the incoming edge list
  // for w is set to {(v, w)} and the shortest path count of w is set
  // to the number of paths that reach{v}.
  void edge_relaxed(edge_descriptor e, const G& g) {
    auto v = source(e, g);
    auto w = target(e, g);
    incoming_[w].clear();
    incoming_[w].push_back(e);
    put(path_count_, w, get(path_count_, v));
  }

  // If an edge e = (v, w) was not relaxed, it may still be the case
  // that we've found more equally-short paths, so include {(v, w)} in
  // the incoming edges of w and add all of the shortest paths to v to
  // the shortest path count of w.
  void edge_not_relaxed(edge_descriptor e, const G& g) {
    auto v = source(e, g);
    auto w = target(e, g);
    auto d_v = get(distance_, v);
    auto d_w = get(distance_, w);
    auto w_e = get(weight_, e);

    closed_plus<decltype(d_v)> combine;
    if (d_w == combine(d_v, w_e)) {
      put(path_count_, w, get(path_count_, w) + get(path_count_, v));
      incoming_[w].push_back(e);
    }
  }

  // Keep track of vertices as they are reached
  void examine_vertex(vertex_descriptor w, const G&) { ordered_vertices_.push(w); }

 private:
  std::stack<vertex_descriptor>& ordered_vertices_;
  WeightMap weight_;
  IncomingMap incoming_;
  DistanceMap distance_;
  PathCountMap path_count_;
};

// Function object that calls Dijkstra's shortest paths algorithm
// using the Dijkstra visitor for the Brandes betweenness centrality
// algorithm.
template <typename WeightMap>
struct brandes_dijkstra_shortest_paths {
  explicit brandes_dijkstra_shortest_paths(WeightMap weight_map) : weight_map_(weight_map) {}

  template <concepts::VertexListGraph G, concepts::MutableLvalueVertexPropertyMap<G> IncomingMap,
            concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadWriteVertexPropertyMap<G> PathCountMap,
            concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
  requires concepts::ReadableEdgePropertyMap<WeightMap, G>
  void operator()(G& g, graph_vertex_descriptor_t<G> s, std::stack<graph_vertex_descriptor_t<G>>& ov,
                  IncomingMap incoming, DistanceMap distance, PathCountMap path_count, VertexIndexMap vertex_index) {
    using Visitor = brandes_dijkstra_visitor<G, WeightMap, IncomingMap, DistanceMap, PathCountMap>;
    Visitor visitor(ov, weight_map_, incoming, distance, path_count);

    using D = property_traits_value_t<DistanceMap>;
    dijkstra_shortest_paths(g, s, make_vector_property_map<graph_vertex_descriptor_t<G>>(vertex_index), distance,
                            weight_map_, vertex_index, std::less<>{}, std::plus<>{}, std::numeric_limits<D>::max(), D{},
                            visitor);
  }

 private:
  WeightMap weight_map_;
};

// Function object that invokes breadth-first search for the
// unweighted form of the Brandes betweenness centrality algorithm.
struct brandes_unweighted_shortest_paths {
  // Customized visitor passed to breadth-first search, which
  // records predecessor and the number of shortest paths to each
  // vertex.
  template <concepts::IncidenceGraph G, concepts::MutableLvalueVertexPropertyMap<G> IncomingMap,
            concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadWriteVertexPropertyMap<G> PathCountMap>
  struct visitor_type : public bfs_visitor<> {
    using edge_descriptor = typename graph_traits<G>::edge_descriptor;
    using vertex_descriptor = typename graph_traits<G>::vertex_descriptor;

    visitor_type(IncomingMap incoming, DistanceMap distance, PathCountMap path_count,
                 std::stack<vertex_descriptor>& ordered_vertices)
        : incoming_(incoming), distance_(distance), path_count_(path_count), ordered_vertices_(ordered_vertices) {}

    // Keep track of vertices as they are reached
    void examine_vertex(vertex_descriptor v, G&) { ordered_vertices_.push(v); }

    // Whenever an edge e = (v, w) is labelled a tree edge, the
    // incoming edge list for w is set to {(v, w)} and the shortest
    // path count of w is set to the number of paths that reach {v}.
    void tree_edge(edge_descriptor e, G& g) {
      vertex_descriptor v = source(e, g);
      vertex_descriptor w = target(e, g);
      put(distance_, w, get(distance_, v) + 1);

      put(path_count_, w, get(path_count_, v));
      incoming_[w].push_back(e);
    }

    // If an edge e = (v, w) is not a tree edge, it may still be the
    // case that we've found more equally-short paths, so include
    // (v, w) in the incoming edge list of w and add all of the
    // shortest paths to v to the shortest path count of w.
    void non_tree_edge(edge_descriptor e, G& g) {
      vertex_descriptor v = source(e, g);
      vertex_descriptor w = target(e, g);
      if (get(distance_, w) == get(distance_, v) + 1) {
        put(path_count_, w, get(path_count_, w) + get(path_count_, v));
        incoming_[w].push_back(e);
      }
    }

   private:
    IncomingMap incoming_;
    DistanceMap distance_;
    PathCountMap path_count_;
    std::stack<vertex_descriptor>& ordered_vertices_;
  };

  template <concepts::VertexListGraph G, concepts::MutableLvalueVertexPropertyMap<G> IncomingMap,
            concepts::ReadWriteVertexPropertyMap<G> DistanceMap, concepts::ReadWriteVertexPropertyMap<G> PathCountMap,
            concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
  void operator()(G& g, typename graph_traits<G>::vertex_descriptor s,
                  std::stack<typename graph_traits<G>::vertex_descriptor>& ov, IncomingMap incoming,
                  DistanceMap distance, PathCountMap path_count, VertexIndexMap vertex_index) {
    using Vertex = typename graph_traits<G>::vertex_descriptor;

    visitor_type<G, IncomingMap, DistanceMap, PathCountMap> visitor(incoming, distance, path_count, ov);

    buffer_queue<Vertex> q;
    breadth_first_visit(g, s, q, visitor,
                        vector_property_map(num_vertices(g), vertex_index, color_traits<default_color_type>::white()));
  }
};

// When the edge centrality map is a dummy property map, no
// initialization is needed.
template <std::ranges::input_range R>
void init_centrality_map(R&& /*unused*/, dummy_property_map /*unused*/) {}

// When we have a real edge centrality map, initialize all of the
// centralities to zero.
template <std::ranges::input_range R, typename Centrality>
void init_centrality_map(R&& keys, Centrality centrality_map) {
  using CentralityValue = property_traits_value_t<Centrality>;
  for (auto k : keys) {
    put(centrality_map, k, CentralityValue(0));
  }
}

// When the edge centrality map is a dummy property map, no update
// is performed.
template <typename Key, typename T>
void update_centrality(dummy_property_map /*unused*/, const Key& /*unused*/, const T& /*unused*/) {}

// When we have a real edge centrality map, add the value to the map
template <typename CentralityMap, typename Key, typename T>
void update_centrality(CentralityMap centrality_map, Key k, const T& x) {
  put(centrality_map, k, get(centrality_map, k) + x);
}

template <std::ranges::input_range R>
void divide_centrality_by_two(R&& /*unused*/, dummy_property_map /*unused*/) {}

template <std::ranges::input_range R, typename Centrality>
void divide_centrality_by_two(R&& keys, Centrality centrality_map) {
  using CentralityValue = property_traits_value_t<Centrality>;
  CentralityValue two(2);
  for (auto k : keys) {
    put(centrality_map, k, get(centrality_map, k) / two);
  }
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadWriteEdgePropertyMap<G> EdgeCentralityMap,
          concepts::MutableLvalueVertexPropertyMap<G> IncomingMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadWriteVertexPropertyMap<G> DependencyMap, concepts::ReadWriteVertexPropertyMap<G> PathCountMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap, typename ShortestPaths>
void brandes_betweenness_centrality_impl(const G& g,
                                         CentralityMap centrality,  // C_B
                                         EdgeCentralityMap edge_centrality_map,
                                         IncomingMap incoming,      // P
                                         DistanceMap distance,      // d
                                         DependencyMap dependency,  // delta
                                         PathCountMap path_count,   // sigma
                                         VertexIndexMap vertex_index, ShortestPaths shortest_paths) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using DependencyType = property_traits_value_t<DependencyMap>;

  // Initialize centrality
  init_centrality_map(vertices(g), centrality);
  init_centrality_map(edges(g), edge_centrality_map);

  std::stack<Vertex> ordered_vertices;
  for (auto s : vertices(g)) {
    // Initialize for this iteration
    for (auto w : vertices(g)) {
      incoming[w].clear();
      put(path_count, w, 0);
      put(dependency, w, 0);
    }
    put(path_count, s, 1);

    // Execute the shortest paths algorithm. This will be either
    // Dijkstra's algorithm or a customized breadth-first search,
    // depending on whether the graph is weighted or unweighted.
    shortest_paths(g, s, ordered_vertices, incoming, distance, path_count, vertex_index);

    while (!ordered_vertices.empty()) {
      Vertex w = ordered_vertices.top();
      ordered_vertices.pop();

      auto w_pcount = static_cast<DependencyType>(get(path_count, w));

      for (auto vw : incoming[w]) {
        Vertex v = source(vw, g);
        auto v_pcount = static_cast<DependencyType>(get(path_count, v));
        auto factor = v_pcount / w_pcount;
        factor *= (DependencyType{1} + get(dependency, w));
        put(dependency, v, get(dependency, v) + factor);
        update_centrality(edge_centrality_map, vw, factor);
      }

      if (w != s) {
        update_centrality(centrality, w, get(dependency, w));
      }
    }
  }

  if constexpr (is_undirected_graph_v<G>) {
    divide_centrality_by_two(vertices(g), centrality);
    divide_centrality_by_two(edges(g), edge_centrality_map);
  }
}

}  // namespace bc_detail

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadWriteEdgePropertyMap<G> EdgeCentralityMap,
          concepts::MutableLvalueVertexPropertyMap<G> IncomingMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadWriteVertexPropertyMap<G> DependencyMap, concepts::ReadWriteVertexPropertyMap<G> PathCountMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
void brandes_betweenness_centrality(const G& g,
                                    CentralityMap centrality,  // C_B
                                    EdgeCentralityMap edge_centrality_map,
                                    IncomingMap incoming,      // P
                                    DistanceMap distance,      // d
                                    DependencyMap dependency,  // delta
                                    PathCountMap path_count,   // sigma
                                    VertexIndexMap vindex) {
  bc_detail::brandes_unweighted_shortest_paths shortest_paths;

  bc_detail::brandes_betweenness_centrality_impl(g, centrality, edge_centrality_map, incoming, distance, dependency,
                                                 path_count, vindex, shortest_paths);
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadWriteEdgePropertyMap<G> EdgeCentralityMap,
          concepts::MutableLvalueVertexPropertyMap<G> IncomingMap, concepts::ReadWriteVertexPropertyMap<G> DistanceMap,
          concepts::ReadWriteVertexPropertyMap<G> DependencyMap, concepts::ReadWriteVertexPropertyMap<G> PathCountMap,
          concepts::ReadableVertexPropertyMap<G> VertexIndexMap, concepts::ReadableEdgePropertyMap<G> WeightMap>
void brandes_betweenness_centrality(const G& g,
                                    CentralityMap centrality,  // C_B
                                    EdgeCentralityMap edge_centrality_map,
                                    IncomingMap incoming,      // P
                                    DistanceMap distance,      // d
                                    DependencyMap dependency,  // delta
                                    PathCountMap path_count,   // sigma
                                    VertexIndexMap vindex, WeightMap weight_map) {
  bc_detail::brandes_dijkstra_shortest_paths<WeightMap> shortest_paths(weight_map);

  bc_detail::brandes_betweenness_centrality_impl(g, centrality, edge_centrality_map, incoming, distance, dependency,
                                                 path_count, vindex, shortest_paths);
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadWriteEdgePropertyMap<G> EdgeCentralityMap, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void brandes_betweenness_centrality(const G& g, CentralityMap centrality, EdgeCentralityMap edge_centrality_map,
                                    VertexIndexMap vindex) {
  using Edge = graph_edge_descriptor_t<G>;
  using OtherCentralityMap =
      std::conditional_t<std::is_same_v<CentralityMap, dummy_property_map>, EdgeCentralityMap, CentralityMap>;
  using OtherCentralityType = property_traits_value_t<OtherCentralityMap>;

  const std::size_t num_v_in_g = num_vertices(g);

  auto incoming = make_vector_property_map<std::vector<Edge>>(num_v_in_g, vindex);
  auto distance = make_vector_property_map<OtherCentralityType>(num_v_in_g, vindex);
  auto dependency = make_vector_property_map<OtherCentralityType>(num_v_in_g, vindex);
  auto path_count = make_vector_property_map<std::size_t>(num_v_in_g, vindex);

  brandes_betweenness_centrality(g, centrality, edge_centrality_map, incoming, distance, dependency, path_count,
                                 vindex);
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap>
void brandes_betweenness_centrality(const G& g, CentralityMap centrality) {
  brandes_betweenness_centrality(g, centrality, dummy_property_map(), get(vertex_index, g));
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadWriteEdgePropertyMap<G> EdgeCentralityMap>
void brandes_betweenness_centrality(const G& g, CentralityMap centrality, EdgeCentralityMap edge_centrality_map) {
  brandes_betweenness_centrality(g, centrality, edge_centrality_map, get(vertex_index, g));
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadWriteEdgePropertyMap<G> EdgeCentralityMap, concepts::ReadableEdgePropertyMap<G> WeightMap,
          concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void brandes_betweenness_centrality_weighted(const G& g, CentralityMap centrality,
                                             EdgeCentralityMap edge_centrality_map, WeightMap weight_map,
                                             VertexIndexMap vindex) {
  using Edge = graph_edge_descriptor_t<G>;
  using OtherCentralityMap =
      std::conditional_t<std::is_same_v<CentralityMap, dummy_property_map>, EdgeCentralityMap, CentralityMap>;
  using OtherCentralityType = property_traits_value_t<OtherCentralityMap>;

  const std::size_t num_v_in_g = num_vertices(g);

  auto incoming = make_vector_property_map<std::vector<Edge>>(num_v_in_g, vindex);
  auto distance = make_vector_property_map<OtherCentralityType>(num_v_in_g, vindex);
  auto dependency = make_vector_property_map<OtherCentralityType>(num_v_in_g, vindex);
  auto path_count = make_vector_property_map<std::size_t>(num_v_in_g, vindex);

  brandes_betweenness_centrality(g, centrality, edge_centrality_map, incoming, distance, dependency, path_count, vindex,
                                 weight_map);
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap, concepts::ReadableVertexIndexMap<G> VertexIndexMap>
void brandes_betweenness_centrality_weighted(const G& g, CentralityMap centrality, WeightMap weight_map,
                                             VertexIndexMap vindex) {
  brandes_betweenness_centrality_weighted(g, centrality, dummy_property_map(), weight_map, vindex);
}

template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap,
          concepts::ReadableEdgePropertyMap<G> WeightMap>
void brandes_betweenness_centrality_weighted(const G& g, CentralityMap centrality, WeightMap weight_map) {
  brandes_betweenness_centrality_weighted(g, centrality, dummy_property_map(), weight_map, get(vertex_index, g));
}

// Converts "absolute" betweenness centrality (as computed by the
// brandes_betweenness_centrality algorithm) in the centrality map
// into "relative" centrality. The result is placed back into the
// given centrality map.
template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap>
void relative_betweenness_centrality(const G& g, CentralityMap centrality) {
  using CentralityValue = property_traits_value_t<CentralityMap>;

  std::size_t n = num_vertices(g);
  CentralityValue factor = CentralityValue{2} / static_cast<CentralityValue>(n * n - 3 * n + 2);
  for (auto v : vertices(g)) {
    put(centrality, v, factor * get(centrality, v));
  }
}

// Compute the central point dominance of a graph.
template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> CentralityMap>
auto central_point_dominance(const G& g, CentralityMap centrality) {
  using CentralityValue = property_traits_value_t<CentralityMap>;

  std::size_t n = num_vertices(g);

  // Find max centrality
  CentralityValue max_centrality{};
  for (auto v : vertices(g)) {
    max_centrality = std::max(max_centrality, get(centrality, v));
  }

  // Compute central point dominance
  CentralityValue sum{};
  for (auto v : vertices(g)) {
    sum += (max_centrality - get(centrality, v));
  }
  return sum / (n - 1);
}

}  // namespace bagl

#endif  // BAGL_BAGL_BETWEENNESS_CENTRALITY_H_
