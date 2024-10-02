// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2006 The Trustees of Indiana University.
// Copyright (C) 2001 Vladimir Prus <ghost@cs.msu.su>
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek, Douglas Gregor
// The mutating functions (add_edge, etc.) were added by Vladimir Prus.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_VECTOR_AS_GRAPH_H_
#define BAGL_BAGL_VECTOR_AS_GRAPH_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <ranges>
#include <utility>
#include <vector>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

// This module implements the VertexListGraph concept using a
// std::vector as the "back-bone" of the graph (the vector *is* the
// graph object). The edge-lists type of the graph is templated, so the
// user can choose any STL container, so long as the value_type of the
// container is convertible to the size_type of the vector. For now any
// graph properties must be stored seperately.

namespace bagl {

struct vector_as_graph_traversal_tag : public virtual vertex_list_graph_tag,
                                       public virtual adjacency_graph_tag,
                                       public virtual incidence_graph_tag {};

struct vector_as_graph_edge {
  std::size_t m_source = std::numeric_limits<std::size_t>::max();
  std::size_t m_target = std::numeric_limits<std::size_t>::max();
};

template <typename EdgeList, typename Alloc>
struct graph_traits<std::vector<EdgeList, Alloc> > {
  using V = typename EdgeList::value_type;
  using vertex_descriptor = std::size_t;
  using edge_descriptor = vector_as_graph_edge;
  using directed_category = directed_tag;
  using edge_parallel_category = allow_parallel_edge_tag;
  using traversal_category = vector_as_graph_traversal_tag;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  using degree_size_type = std::size_t;

  static std::size_t null_vertex() { return std::numeric_limits<std::size_t>::max(); }
};

template <typename EdgeList, typename Alloc>
auto out_edges(std::size_t u, const std::vector<EdgeList, Alloc>& g) {
  return std::ranges::ref_view(g[u]) | std::views::transform([u](std::size_t v) { return vector_as_graph_edge{u, v}; });
}

template <typename EdgeList, typename Alloc>
std::size_t out_degree(std::size_t u, const std::vector<EdgeList, Alloc>& g) {
  return g[u].size();
}

template <typename EdgeList, typename Alloc>
auto adjacent_vertices(std::size_t u, const std::vector<EdgeList, Alloc>& g) {
  return std::ranges::ref_view(g[u]);
}

template <typename EdgeList, typename Alloc>
std::size_t source(vector_as_graph_edge e, const std::vector<EdgeList, Alloc>& /*g*/) {
  return e.m_source;
}

template <typename EdgeList, typename Alloc>
std::size_t target(vector_as_graph_edge e, const std::vector<EdgeList, Alloc>& /*g*/) {
  return e.m_target;
}

template <typename EdgeList, typename Alloc>
auto vertices(const std::vector<EdgeList, Alloc>& v) {
  return std::ranges::iota_view<std::size_t, std::size_t>(0, v.size());
}

template <typename EdgeList, typename Alloc>
std::size_t num_vertices(const std::vector<EdgeList, Alloc>& v) {
  return v.size();
}

template <typename EdgeList, typename Alloc>
auto add_edge(std::size_t u, std::size_t v, std::vector<EdgeList, Alloc>& g) {
  g[u].insert(g[u].end(), v);
  return std::pair{vector_as_graph_edge{u, v}, true};
}

template <typename EdgeList, typename Alloc>
auto edge(std::size_t u, std::size_t v, std::vector<EdgeList, Alloc>& g) {
  for (auto i : g[u]) {
    if (*i == v) {
      return std::pair{vector_as_graph_edge{u, v}, true};
    }
  }
  return std::pair{vector_as_graph_edge{}, false};
}

template <typename EdgeList, typename Alloc>
void remove_edge(std::size_t u, std::size_t v, std::vector<EdgeList, Alloc>& g) {
  auto i = std::remove(g[u].begin(), g[u].end(), v);
  if (i != g[u].end()) {
    g[u].erase(i, g[u].end());
  }
}

template <typename EdgeList, typename Alloc>
void remove_edge(vector_as_graph_edge e, std::vector<EdgeList, Alloc>& g) {
  remove_edge(e.m_source, e.m_target, g);
}

template <typename EdgeList, typename Alloc, typename Predicate>
void remove_edge_if(Predicate p, std::vector<EdgeList, Alloc>& g) {
  for (std::size_t u = 0; u < g.size(); ++u) {
    auto i = std::remove_if(g[u].begin(), g[u].end(), [u, &p](std::size_t v) { return p(vector_as_graph_edge{u, v}); });
    if (i != g[u].end()) {
      g[u].erase(i, g[u].end());
    }
  }
}

template <typename EdgeList, typename Alloc>
std::size_t add_vertex(std::vector<EdgeList, Alloc>& g) {
  g.resize(g.size() + 1);
  return g.size() - 1;
}

template <typename EdgeList, typename Alloc>
void clear_vertex(std::size_t u, std::vector<EdgeList, Alloc>& g) {
  g[u].clear();
  for (std::size_t i = 0; i < g.size(); ++i) {
    remove_edge(i, u, g);
  }
}

template <typename EdgeList, typename Alloc>
void remove_vertex(std::size_t u, std::vector<EdgeList, Alloc>& g) {
  clear_vertex(u, g);
  g.erase(g.begin() + u);
  for (std::size_t i = 0; i < g.size(); ++i) {
    for (auto& v : g[i]) {
      // after clear_vertex v is never equal to u
      if (v > u) {
        --v;
      }
    }
  }
}

template <typename EdgeList, typename Alloc>
struct property_map<std::vector<EdgeList, Alloc>, vertex_index_t> {
  using type = identity_property_map;
  using const_type = type;
};

template <typename EdgeList, typename Alloc>
identity_property_map get(vertex_index_t /*unused*/, const std::vector<EdgeList, Alloc>& /*unused*/) {
  return {};
}

template <typename EdgeList, typename Alloc>
identity_property_map get(vertex_index_t /*unused*/, std::vector<EdgeList, Alloc>& /*unused*/) {
  return {};
}
}  // namespace bagl

#endif  // BAGL_BAGL_VECTOR_AS_GRAPH_H_
