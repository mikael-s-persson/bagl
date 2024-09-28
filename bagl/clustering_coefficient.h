// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CLUSTERING_COEFFICIENT_H_
#define BAGL_BAGL_CLUSTERING_COEFFICIENT_H_

#include <iterator>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/lookup_edge.h"

namespace bagl {
namespace clustering_detail {
template <concepts::AdjacencyGraph G>
std::size_t count_edges(const G& g, graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v) {
  if constexpr (is_undirected_graph_v<G>) {
    return lookup_edge(u, v, g).second ? 1 : 0;
  } else {
    return (lookup_edge(u, v, g).second ? 1 : 0) + (lookup_edge(v, u, g).second ? 1 : 0);
  }
}
}  // namespace clustering_detail

template <concepts::AdjacencyGraph G, typename Vertex>
std::size_t num_paths_through_vertex(const G& g, Vertex v) {
  const std::size_t k = num_adjacent_vertices(v, g);
  if constexpr (is_undirected_graph_v<G>) {
    return (k * (k - 1)) / 2;
  } else {
    return k * (k - 1);
  }
}

template <concepts::AdjacencyGraph G>
std::size_t num_triangles_on_vertex(const G& g, graph_vertex_descriptor_t<G> v) {
  std::size_t count = 0;
  auto adj_rg = adjacent_vertices(v, g);
  for (auto i = adj_rg.begin(), i_end = adj_rg.end(); i != i_end; ++i) {
    for (auto j = std::next(i); j != i_end; ++j) {
      count += clustering_detail::count_edges(g, *i, *j);
    }
  }
  return count;
}

template <typename T, concepts::AdjacencyGraph G>
T clustering_coefficient(const G& g, graph_vertex_descriptor_t<G> v) {
  const T zero = T{0};
  const T routes = T{num_paths_through_vertex(g, v)};
  return (routes > zero) ? T{num_triangles_on_vertex(g, v)} / routes : zero;
}

template <concepts::AdjacencyGraph G>
double clustering_coefficient(const G& g, graph_vertex_descriptor_t<G> v) {
  return clustering_coefficient<double>(g, v);
}

template <concepts::VertexListGraph G, concepts::WritableVertexPropertyMap<G> ClusteringMap>
property_traits_value_t<ClusteringMap> all_clustering_coefficients(const G& g, ClusteringMap cm) {
  using Coefficient = property_traits_value_t<ClusteringMap>;

  Coefficient sum = Coefficient{0};
  for (auto v : vertices(g)) {
    const auto cc = clustering_coefficient<Coefficient>(g, v);
    put(cm, v, cc);
    sum += cc;
  }
  return sum / Coefficient{num_vertices(g)};
}

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> ClusteringMap>
property_traits_value_t<ClusteringMap> mean_clustering_coefficient(const G& g, ClusteringMap cm) {
  using Coefficient = property_traits_value_t<ClusteringMap>;

  Coefficient cc = Coefficient{0};
  for (auto v : vertices(g)) {
    cc += get(cm, v);
  }
  return cc / Coefficient{num_vertices(g)};
}

}  // namespace bagl

#endif  // BAGL_BAGL_CLUSTERING_COEFFICIENT_H_
