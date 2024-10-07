// Copyright 2005 The Trustees of Indiana University.
//  Authors: Alex Breuer
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_STATS_H_
#define BAGL_BAGL_GRAPH_STATS_H_

#include <algorithm>
#include <cassert>
#include <map>
#include <ranges>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/subgraph.h"

namespace bagl {

template <concepts::EdgeListGraph G>
std::size_t num_dup_edges(G& g) {
  std::vector<graph_edge_descriptor_t<G>> all_edges;
  all_edges.reserve(num_edges(g));
  for (auto e : edges(g)) {
    all_edges.push_back(e);
  }

  std::ranges::sort(all_edges, [&g](const auto& lhs_e, const auto& rhs_e) {
    return std::pair{source(lhs_e, g), target(lhs_e, g)} < std::pair{source(rhs_e, g), target(rhs_e, g)};
  });
  auto uniq_all_edges = std::ranges::unique(all_edges);

  return num_edges(g) - uniq_all_edges.size();
}

template <concepts::VertexListGraph G>
requires concepts::AdjacencyGraph<G> std::map<std::size_t, std::size_t> dup_edge_dist(G& g) {
  std::map<std::size_t, std::size_t> dist;
  std::vector<graph_vertex_descriptor_t<G>> front_neighbors;
  for (auto v : vertices(g)) {
    front_neighbors.clear();
    for (auto u : adjacent_vertices(v, g)) {
      front_neighbors.push_back(u);
    }

    std::ranges::sort(front_neighbors);
    auto uniq_front_neighbors = std::ranges::unique(front_neighbors);
    dist[out_degree(v, g) - uniq_front_neighbors.size()] += 1;
  }
  return dist;
}

template <concepts::VertexListGraph G>
requires concepts::IncidenceGraph<G> std::map<std::size_t, std::size_t> degree_dist(G& g) {
  std::map<std::size_t, std::size_t> dist;
  for (auto v : vertices(g)) {
    dist[out_degree(v, g)] += 1;
  }

  return dist;
}

template <concepts::VertexListGraph G>
requires concepts::IncidenceGraph<G> std::map<std::size_t, double> weight_degree_dist(G& g) {
  std::map<std::size_t, double> dist;
  std::map<std::size_t, double> n;
  using edge_map_type = property_map_const_t<G, edge_weight_t>;
  using edge_weight_type = property_traits_value_t<edge_map_type>;

  auto em = get(edge_weight, g);
  for (auto v : vertices(g)) {
    edge_weight_type tmp = 0;
    for (auto e : out_edges(v, g)) {
      tmp += em[e];
    }
    n[out_degree(v, g)] += 1.;
    dist[out_degree(v, g)] += tmp;
  }

  for (auto& [deg, val] : dist) {
    assert(n[deg] != 0);
    val /= n[deg];
  }

  return dist;
}

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPH_STATS_H_
