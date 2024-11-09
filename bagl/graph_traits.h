// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_TRAITS_H_
#define BAGL_BAGL_GRAPH_TRAITS_H_

#include <tuple>
#include <type_traits>
#include <utility>

#include "bagl/detail/container_generators.h"
#include "bagl/has_trait_member.h"
#include "bagl/property.h"

namespace bagl {

namespace graph_traits_detail {

BAGL_GRAPH_HAS_TRAIT_MEMBER(vertices_size_type, std::size_t)
BAGL_GRAPH_HAS_TRAIT_MEMBER(edges_size_type, std::size_t)
BAGL_GRAPH_HAS_TRAIT_MEMBER(degree_size_type, std::size_t)

BAGL_GRAPH_HAS_TRAIT_MEMBER(graph_property_type, no_property)
BAGL_GRAPH_HAS_TRAIT_MEMBER(edge_property_type, no_property)
BAGL_GRAPH_HAS_TRAIT_MEMBER(vertex_property_type, no_property)

}  // namespace graph_traits_detail

template <typename G>
struct graph_traits {
  using vertex_descriptor = typename G::vertex_descriptor;
  using edge_descriptor = typename G::edge_descriptor;

  using directed_category = typename G::directed_category;
  using edge_parallel_category = typename G::edge_parallel_category;
  using traversal_category = typename G::traversal_category;

  using vertices_size_type = graph_traits_detail::get_vertices_size_type_or_not<G>;
  using edges_size_type = graph_traits_detail::get_edges_size_type_or_not<G>;
  using degree_size_type = graph_traits_detail::get_degree_size_type_or_not<G>;

  static vertex_descriptor null_vertex() { return G::null_vertex(); }
};

template <typename G>
struct graph_traits<const G> : graph_traits<G> {};

template <typename G>
struct graph_traits<G&> : graph_traits<G> {};

// Shorthands to avoid the "typename graph_traits<...>::foo_bar" syntax everywhere.
template <typename G>
using graph_vertex_descriptor_t = typename graph_traits<G>::vertex_descriptor;
template <typename G>
using graph_edge_descriptor_t = typename graph_traits<G>::edge_descriptor;
template <typename G>
using graph_directed_category_t = typename graph_traits<G>::directed_category;
template <typename G>
using graph_edge_parallel_category_t = typename graph_traits<G>::edge_parallel_category;
template <typename G>
using graph_traversal_category_t = typename graph_traits<G>::traversal_category;
template <typename G>
using graph_vertices_size_type_t = typename graph_traits<G>::vertices_size_type;
template <typename G>
using graph_edges_size_type_t = typename graph_traits<G>::edges_size_type;
template <typename G>
using graph_degree_size_type_t = typename graph_traits<G>::degree_size_type;

template <typename G>
using graph_vertex_range_t = decltype(vertices(std::declval<G>()));
template <typename G>
using graph_out_edge_range_t = decltype(out_edges(std::declval<graph_vertex_descriptor_t<G>>(), std::declval<G>()));
template <typename G>
using graph_in_edge_range_t = decltype(in_edges(std::declval<graph_vertex_descriptor_t<G>>(), std::declval<G>()));
template <typename G>
using graph_edge_range_t = decltype(edges(std::declval<G>()));
template <typename G>
using graph_adjacency_range_t = decltype(adjacent_vertices(std::declval<graph_vertex_descriptor_t<G>>(), std::declval<G>()));

// Hash functors for vertex or edge descriptors.
// Specialize the graph_descriptor_hash class for descriptors not supported in container_generators.h.
template <typename Desc>
struct graph_descriptor_hash {
  using type = container_detail::desc_hasher;
};
template <typename Desc>
using graph_descriptor_hash_t = typename graph_descriptor_hash<Desc>::type;

// directed_category tags
struct directed_tag {};
struct undirected_tag {};
struct bidirectional_tag : public directed_tag {};

/** Return true if the given graph is directed. */
template <typename G>
constexpr bool is_directed_graph_v = std::is_convertible_v<typename graph_traits<G>::directed_category, directed_tag>;

template <typename G>
bool is_directed(const G&) {
  return is_directed_graph_v<G>;
}

/** Return true if the given graph is undirected. */
template <typename G>
constexpr bool is_undirected_graph_v = !is_directed_graph_v<G>;

template <typename G>
bool is_undirected(const G&) {
  return is_undirected_graph_v<G>;
}

// edge_parallel_category tags
struct allow_parallel_edge_tag {};
struct disallow_parallel_edge_tag {};

template <typename G>
constexpr bool allows_parallel_edges_v =
    std::is_convertible_v<typename graph_traits<G>::edge_parallel_category, allow_parallel_edge_tag>;

template <typename G>
bool allows_parallel_edges(const G&) {
  return allows_parallel_edges_v<G>;
}

template <typename G>
constexpr bool is_multigraph_v = allows_parallel_edges_v<G>;

// traversal_category tags
struct incidence_graph_tag {};
struct adjacency_graph_tag {};
struct bidirectional_graph_tag : virtual incidence_graph_tag {};
struct vertex_list_graph_tag {};
struct edge_list_graph_tag {};
struct adjacency_matrix_tag {};

// Parallel traversal_category tags
struct distributed_graph_tag {};
struct distributed_vertex_list_graph_tag {};
struct distributed_edge_list_graph_tag {};

/** @name Traversal Category Traits
 * These traits classify graph types by their supported methods of
 * vertex and edge traversal.
 */
//@{
template <typename G>
constexpr bool is_incidence_graph_v =
    std::is_convertible_v<typename graph_traits<G>::traversal_category, incidence_graph_tag>;

template <typename G>
constexpr bool is_bidirectional_graph_v =
    std::is_convertible_v<typename graph_traits<G>::traversal_category, bidirectional_graph_tag>;

template <typename G>
constexpr bool is_vertex_list_graph_v =
    std::is_convertible_v<typename graph_traits<G>::traversal_category, vertex_list_graph_tag>;

template <typename G>
constexpr bool is_edge_list_graph_v =
    std::is_convertible_v<typename graph_traits<G>::traversal_category, edge_list_graph_tag>;

template <typename G>
constexpr bool is_adjacency_matrix_v =
    std::is_convertible_v<typename graph_traits<G>::traversal_category, adjacency_matrix_tag>;
//@}

/** @name Directed Graph Traits
 * These metafunctions are used to fully classify directed vs. undirected
 * graphs. Recall that an undirected graph is also bidirectional, but it
 * cannot be both undirected and directed at the same time.
 */
//@{
template <typename G>
constexpr bool is_directed_unidirectional_graph_v = is_directed_graph_v<G> && !is_bidirectional_graph_v<G>;

template <typename G>
constexpr bool is_directed_bidirectional_graph = is_directed_graph_v<G>&& is_bidirectional_graph_v<G>;
//@}

template <typename G>
using graph_property_type = graph_traits_detail::get_graph_property_type_or_not<G>;
template <typename G>
using edge_property_type = graph_traits_detail::get_edge_property_type_or_not<G>;
template <typename G>
using vertex_property_type = graph_traits_detail::get_vertex_property_type_or_not<G>;

template <typename G>
using graph_bundle_type = typename G::graph_bundled;

template <typename G>
using vertex_bundle_type = typename G::vertex_bundled;

template <typename G>
using edge_bundle_type = typename G::edge_bundled;

namespace graph_traits_detail {

template <typename G, typename Descriptor>
struct bundled_result {
  using Vertex = typename graph_traits<G>::vertex_descriptor;

 public:
  using type = std::conditional_t<std::is_same_v<Descriptor, Vertex>, vertex_bundle_type<G>, edge_bundle_type<G>>;
};

template <typename G>
struct bundled_result<G, graph_bundle_t> {
  using type = graph_bundle_type<G>;
};

}  // namespace graph_traits_detail

namespace graph_traits_detail {
// A helper metafunction for determining whether or not a type is
// bundled.
template <typename T>
constexpr bool is_no_bundle_v = std::is_same_v<T, no_property>;
}  // namespace graph_traits_detail

/** @name Graph Property Traits
 * These metafunctions (along with those above), can be used to access the
 * vertex and edge properties (bundled or otherwise) of vertices and
 * edges.
 */
//@{
template <typename G>
constexpr bool has_graph_property_v = !std::is_same_v<graph_property_type<G>, no_property>;

template <typename G>
constexpr bool has_bundled_graph_property_v = !std::is_same_v<graph_bundle_type<G>, no_property>;

template <typename G>
constexpr bool has_vertex_property_v = !std::is_same_v<vertex_property_type<G>, no_property>;

template <typename G>
constexpr bool has_bundled_vertex_property_v = !std::is_same_v<vertex_bundle_type<G>, no_property>;

template <typename G>
constexpr bool has_edge_property_v = !std::is_same_v<edge_property_type<G>, no_property>;

template <typename G>
constexpr bool has_bundled_edge_property_v = !std::is_same_v<edge_bundle_type<G>, no_property>;
//@}

/* Some helper functions for dealing with pairs as edges */
template <class T, class G>
T source(std::pair<T, T> p, const G&) {
  return p.first;
}

template <class T, class G>
T target(std::pair<T, T> p, const G&) {
  return p.second;
}

template <typename G>
auto get_default_starting_vertex(const G& g) {
  auto v_rg = vertices(g);
  return (v_rg.empty()) ? graph_traits<G>::null_vertex() : *v_rg.begin();
}

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPH_TRAITS_H_
