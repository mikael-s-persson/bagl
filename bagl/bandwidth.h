// Copyright 2001 Jeremy Siek
// Copyright 2002 Marc Wintermantel
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BANDWIDTH_H_
#define BAGL_BAGL_BANDWIDTH_H_

#include <algorithm>  // for std::min and std::max
#include <cmath>      // for std::abs

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

template <concepts::IncidenceGraph G, concepts::ReadablePropertyMap<graph_vertex_descriptor_t<G>> VertexIndexMap>
std::size_t ith_bandwidth(graph_vertex_descriptor_t<G> i, const G& g, VertexIndexMap index) {
  using std::abs;
  std::size_t b = 0;
  for (auto e : out_edges(i, g)) {
    int f_i = get(index, i);
    int f_j = get(index, target(e, g));
    b = std::max(b, static_cast<std::size_t>(abs(f_i - f_j)));
  }
  return b;
}

template <concepts::IncidenceGraph G>
requires concepts::VertexIndexGraph<G> std::size_t ith_bandwidth(graph_vertex_descriptor_t<G> i, const G& g) {
  return ith_bandwidth(i, g, get(vertex_index, g));
}

template <concepts::EdgeListGraph G, concepts::ReadablePropertyMap<graph_vertex_descriptor_t<G>> VertexIndexMap>
std::size_t bandwidth(const G& g, VertexIndexMap index) {
  using std::abs;
  std::size_t b = 0;
  for (auto e : edges(g)) {
    int f_i = get(index, source(e, g));
    int f_j = get(index, target(e, g));
    b = std::max(b, static_cast<std::size_t>(abs(f_i - f_j)));
  }
  return b;
}

template <concepts::EdgeListGraph G>
std::size_t bandwidth(const G& g) {
  return bandwidth(g, get(vertex_index, g));
}

template <typename Graph, typename VertexIndexMap>
std::size_t edgesum(const Graph& g, VertexIndexMap index_map) {
  using std::abs;
  std::size_t sum = 0;
  for (auto e : edges(g)) {
    auto f_u = get(index_map, source(e, g));
    auto f_v = get(index_map, target(e, g));
    sum += abs(f_u - f_v);
  }
  return sum;
}

}  // namespace bagl

#endif  // BAGL_BAGL_BANDWIDTH_H_
