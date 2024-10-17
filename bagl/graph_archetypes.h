// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_ARCHETYPES_H_
#define BAGL_BAGL_GRAPH_ARCHETYPES_H_

#include <ranges>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

namespace graph_archetypes_detail {
struct null_graph_archetype {
  struct traversal_category {};
};
}  // namespace graph_archetypes_detail

// Integer parameter is only used to create multiple distinct types in tests.
template <int I = 0>
struct semiregular_archetype {
  bool operator==(semiregular_archetype<I>) const { return true; }
  bool operator!=(semiregular_archetype<I>) const { return false; }
};

//===========================================================================
template <typename Vertex, typename Directed, typename ParallelCategory,
          typename Base = graph_archetypes_detail::null_graph_archetype>
struct incidence_graph_archetype : public Base {
  using base_trav_cat = typename Base::traversal_category;
  struct traversal_category : public incidence_graph_tag, public base_trav_cat {};
  using vertex_descriptor = Vertex;
  using degree_size_type = std::size_t;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  struct edge_descriptor {
    edge_descriptor() = default;
    explicit edge_descriptor(const graph_archetypes_detail::null_graph_archetype& /*unused*/) {}
    bool operator==(const edge_descriptor& /*unused*/) const { return false; }
    bool operator!=(const edge_descriptor& /*unused*/) const { return false; }
  };

  using directed_category = Directed;
  using edge_parallel_category = ParallelCategory;

  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
};
template <typename V, typename D, typename P, typename B>
V source(const typename incidence_graph_archetype<V, D, P, B>::edge_descriptor& /*unused*/,
         const incidence_graph_archetype<V, D, P, B>& g) {
  return V{};
}
template <typename V, typename D, typename P, typename B>
V target(const typename incidence_graph_archetype<V, D, P, B>::edge_descriptor& /*unused*/,
         const incidence_graph_archetype<V, D, P, B>& g) {
  return V{};
}

template <typename V, typename D, typename P, typename B>
auto out_edges(const V& /*unused*/, const incidence_graph_archetype<V, D, P, B>& /*unused*/) {
  return std::ranges::empty_view<typename incidence_graph_archetype<V, D, P, B>::edge_descriptor>{};
}

template <typename V, typename D, typename P, typename B>
int out_degree(const V& /*unused*/, const incidence_graph_archetype<V, D, P, B>& /*unused*/) {
  return 0;
}

//===========================================================================
template <typename Vertex, typename Directed, typename ParallelCategory,
          typename Base = graph_archetypes_detail::null_graph_archetype>
struct adjacency_graph_archetype : public Base {
  using base_trav_cat = typename Base::traversal_category;
  struct traversal_category : public adjacency_graph_tag, public base_trav_cat {};
  using vertex_descriptor = Vertex;
  using degree_size_type = std::size_t;
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  using edge_descriptor = void;

  using directed_category = Directed;
  using edge_parallel_category = ParallelCategory;

  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
};

template <typename V, typename D, typename P, typename B>
auto adjacent_vertices(const V& /*unused*/, const adjacency_graph_archetype<V, D, P, B>& /*unused*/) {
  return std::ranges::empty_view<typename adjacency_graph_archetype<V, D, P, B>::vertex_descriptor>{};
}

template <typename V, typename D, typename P, typename B>
int out_degree(const V& /*unused*/, const adjacency_graph_archetype<V, D, P, B>& /*unused*/) {
  return 0;
}

//===========================================================================
template <typename Vertex, typename Directed, typename ParallelCategory,
          typename Base = graph_archetypes_detail::null_graph_archetype>
struct vertex_list_graph_archetype : public Base {
  using Incidence = incidence_graph_archetype<Vertex, Directed, ParallelCategory>;
  using Adjacency = adjacency_graph_archetype<Vertex, Directed, ParallelCategory>;

  using base_trav_cat = typename Base::traversal_category;
  struct traversal_category : public vertex_list_graph_tag, public base_trav_cat {};
  using vertex_descriptor = Vertex;
  using degree_size_type = std::size_t;
  using edge_descriptor = typename Incidence::edge_descriptor;

  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;

  using directed_category = Directed;
  using edge_parallel_category = ParallelCategory;

  static vertex_descriptor null_vertex() { return vertex_descriptor(); }
};

template <typename V, typename D, typename P, typename B>
auto vertices(const vertex_list_graph_archetype<V, D, P, B>& /*unused*/) {
  return std::ranges::empty_view<typename vertex_list_graph_archetype<V, D, P, B>::vertex_descriptor>{};
}

template <typename V, typename D, typename P, typename B>
int num_vertices(const vertex_list_graph_archetype<V, D, P, B>& /*unused*/) {
  return 0;
}

// ambiguously inherited from incidence graph and adjacency graph
template <typename V, typename D, typename P, typename B>
int out_degree(const V& /*unused*/, const vertex_list_graph_archetype<V, D, P, B>& /*unused*/) {
  return 0;
}

//===========================================================================

struct property_graph_archetype_tag {};

template <typename GraphArchetype, typename Property, typename ValueArch>
struct property_graph_archetype : public GraphArchetype {
  using graph_tag = property_graph_archetype_tag;
  using vertex_property_type = ValueArch;
  using edge_property_type = ValueArch;
};

struct choose_edge_property_map_archetype {
  template <typename Graph, typename Property, typename Tag>
  struct bind_ {
    using type = mutable_lvalue_property_map_archetype<typename Graph::edge_descriptor, Property>;
    using const_type = lvalue_property_map_archetype<typename Graph::edge_descriptor, Property>;
  };
};
template <>
struct edge_property_selector<property_graph_archetype_tag> {
  using type = choose_edge_property_map_archetype;
};

struct choose_vertex_property_map_archetype {
  template <typename Graph, typename Property, typename Tag>
  struct bind_ {
    using type = mutable_lvalue_property_map_archetype<typename Graph::vertex_descriptor, Property>;
    using const_type = lvalue_property_map_archetype<typename Graph::vertex_descriptor, Property>;
  };
};

template <>
struct vertex_property_selector<property_graph_archetype_tag> {
  using type = choose_vertex_property_map_archetype;
};

template <typename G, typename P, typename V>
property_map_t<property_graph_archetype<G, P, V>, P> get(P /*unused*/, property_graph_archetype<G, P, V>& /*unused*/) {
  return {};
}

template <typename G, typename P, typename V>
property_map_const_t<property_graph_archetype<G, P, V>, P> get(P /*unused*/,
                                                               const property_graph_archetype<G, P, V>& /*unused*/) {
  return {};
}

template <typename G, typename P, typename K, typename V>
auto get(P p, const property_graph_archetype<G, P, V>& g, K k) {
  return get(get(p, g), k);
}

template <typename G, typename P, typename V, typename Key>
void put(P p, property_graph_archetype<G, P, V>& g, const Key& key, const V& value) {
  auto pmap = get(p, g);
  put(pmap, key, value);
}

struct color_value_archetype {
  color_value_archetype() = default;
  bool operator==(const color_value_archetype& /*unused*/) const { return true; }
  bool operator!=(const color_value_archetype& /*unused*/) const { return true; }
};
template <>
struct color_traits<color_value_archetype> {
  static color_value_archetype white() { return color_value_archetype{}; }
  static color_value_archetype gray() { return color_value_archetype{}; }
  static color_value_archetype black() { return color_value_archetype{}; }
};

template <typename T>
class buffer_archetype {
 public:
  using value_type = T;
  void push(const T& /*unused*/) {}
  void pop() {}
  [[nodiscard]] T& top() { return value; }
  [[nodiscard]] const T& top() const { return value; }
  [[nodiscard]] bool empty() const { return true; }
  T value{};
};

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPH_ARCHETYPES_H_
