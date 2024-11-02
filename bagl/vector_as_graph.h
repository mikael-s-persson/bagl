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
  std::size_t u = std::numeric_limits<std::size_t>::max();
  std::size_t v = std::numeric_limits<std::size_t>::max();

  bool operator==(const vector_as_graph_edge& rhs) const { return u == rhs.u && v == rhs.v; }
  bool operator!=(const vector_as_graph_edge& rhs) const { return !(*this == rhs); }
};

template <std::ranges::random_access_range Vector>
struct vector_as_graph {
  Vector& data;
  explicit vector_as_graph(Vector& a_data) : data(a_data) {}

  using vertex_descriptor = std::size_t;
  using edge_descriptor = vector_as_graph_edge;
  using directed_category = directed_tag;
  using edge_parallel_category = allow_parallel_edge_tag;
  using traversal_category = vector_as_graph_traversal_tag;

  static std::size_t null_vertex() { return std::numeric_limits<std::size_t>::max(); }

  // Indexing operator. Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
  const auto& operator[](std::size_t v) const { return data[v]; }

  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  const auto& operator[](const edge_descriptor& e) const { return data[e.u][e.v]; }
};

template <typename Vector>
auto out_edges(std::size_t u, const vector_as_graph<Vector>& g) {
  return std::ranges::ref_view(g.data[u]) | std::views::transform([u](std::size_t v) {
           return vector_as_graph_edge{u, v};
         });
}

template <typename Vector>
std::size_t out_degree(std::size_t u, const vector_as_graph<Vector>& g) {
  return g.data[u].size();
}

template <typename Vector>
auto adjacent_vertices(std::size_t u, const vector_as_graph<Vector>& g) {
  return std::ranges::ref_view(g.data[u]);
}

template <typename Vector>
std::size_t source(vector_as_graph_edge e, const vector_as_graph<Vector>& /*g*/) {
  return e.u;
}

template <typename Vector>
std::size_t target(vector_as_graph_edge e, const vector_as_graph<Vector>& /*g*/) {
  return e.v;
}

template <typename Vector>
auto vertices(const vector_as_graph<Vector>& g) {
  return std::ranges::iota_view<std::size_t, std::size_t>(0, g.data.size());
}

template <typename Vector>
std::size_t num_vertices(const vector_as_graph<Vector>& g) {
  return g.data.size();
}

template <typename Vector>
auto add_edge(std::size_t u, std::size_t v, vector_as_graph<Vector>& g) {
  g.data[u].insert(g.data[u].end(), v);
  return std::pair{vector_as_graph_edge{u, v}, true};
}

template <typename Vector>
auto edge(std::size_t u, std::size_t v, const vector_as_graph<Vector>& g) {
  for (auto i : g.data[u]) {
    if (i == v) {
      return std::pair{vector_as_graph_edge{u, v}, true};
    }
  }
  return std::pair{vector_as_graph_edge{}, false};
}

template <typename Vector>
void remove_edge(std::size_t u, std::size_t v, vector_as_graph<Vector>& g) {
  auto i = std::remove(g.data[u].begin(), g.data[u].end(), v);
  if (i != g.data[u].end()) {
    g.data[u].erase(i, g.data[u].end());
  }
}

template <typename Vector>
void remove_edge(vector_as_graph_edge e, vector_as_graph<Vector>& g) {
  remove_edge(e.u, e.v, g);
}

template <typename Vector, typename Predicate>
void remove_edge_if(Predicate p, vector_as_graph<Vector>& g) {
  for (std::size_t u = 0; u < g.data.size(); ++u) {
    auto& elist = g.data[u];
    auto i = std::remove_if(elist.begin(), elist.end(), [u, &p](std::size_t v) {
      return p(vector_as_graph_edge{u, v});
    });
    if (i != elist.end()) {
      elist.erase(i, elist.end());
    }
  }
}

template <typename Vector>
std::size_t add_vertex(vector_as_graph<Vector>& g) {
  g.data.resize(g.data.size() + 1);
  return g.data.size() - 1;
}

template <typename Vector>
void clear_vertex(std::size_t u, vector_as_graph<Vector>& g) {
  g.data[u].clear();
  for (std::size_t i = 0; i < g.data.size(); ++i) {
    remove_edge(i, u, g);
  }
}

template <typename Vector>
void remove_vertex(std::size_t u, vector_as_graph<Vector>& g) {
  clear_vertex(u, g);
  g.data.erase(g.data.begin() + u);
  for (auto& elist : g.data) {
    for (auto& v : elist) {
      // after clear_vertex v is never equal to u
      if (v > u) {
        --v;
      }
    }
  }
}

template <typename Vector>
struct property_map<vector_as_graph<Vector>, vertex_index_t> {
  using type = identity_property_map;
  using const_type = type;
};

template <typename Vector>
identity_property_map get(vertex_index_t /*unused*/, const vector_as_graph<Vector>& /*unused*/) {
  return {};
}

template <typename Vector>
std::size_t get(vertex_index_t /*unused*/, const vector_as_graph<Vector>& /*unused*/, std::size_t v) {
  return v;
}

}  // namespace bagl

#endif  // BAGL_BAGL_VECTOR_AS_GRAPH_H_
