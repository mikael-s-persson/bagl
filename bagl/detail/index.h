// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_INDEX_H_
#define BAGL_BAGL_DETAIL_INDEX_H_

#include <type_traits>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"

// The structures in this module are responsible for selecting and defining
// types for accessing a builting index map. Note that the selection of these
// types requires the Graph parameter to model either VertexIndexGraph or
// EdgeIndexGraph.

namespace bagl::index_detail {
template <typename Graph>
struct vertex_indexer {
  using index_type = vertex_index_t;
  using map_type = property_map_t<Graph, vertex_index_t>;
  using const_map_type = property_map_const_t<Graph, vertex_index_t>;
  using value_type = property_traits_value_t<map_type>;
  using key_type = graph_vertex_descriptor_t<Graph>;

  static const_map_type index_map(const Graph& g) { return get(vertex_index, g); }
  static map_type index_map(Graph& g) { return get(vertex_index, g); }

  static value_type index(key_type k, const Graph& g) { return get(vertex_index, g, k); }
};

template <typename Graph>
struct edge_indexer {
  using index_type = edge_index_t;
  using map_type = property_map_t<Graph, edge_index_t>;
  using const_map_type = property_map_const_t<Graph, edge_index_t>;
  using value_type = property_traits_value_t<map_type>;
  using key_type = graph_edge_descriptor_t<Graph>;

  static const_map_type index_map(const Graph& g) { return get(edge_index, g); }
  static map_type index_map(Graph& g) { return get(edge_index, g); }

  static value_type index(key_type k, const Graph& g) { return get(edge_index, g, k); }
};

// NOTE: The Graph parameter MUST be a model of VertexIndexGraph or
// VertexEdgeGraph - whichever type Key is selecting.
template <typename Graph, typename Key>
struct choose_indexer {
  using indexer_type = std::conditional_t<std::is_same_v<Key, graph_vertex_descriptor_t<Graph>>, vertex_indexer<Graph>,
                                          edge_indexer<Graph>>;
  using index_type = typename indexer_type::index_type;
};
}  // namespace bagl::index_detail

#endif  // BAGL_BAGL_DETAIL_INDEX_H_
