// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2009, Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_CONCEPTS_H_
#define BAGL_BAGL_GRAPH_CONCEPTS_H_

#include <concepts>
#include <ranges>
#include <type_traits>

#include "bagl/buffer_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/numeric_values.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl::concepts {

template <typename G>
concept Graph = std::regular<graph_vertex_descriptor_t<G>> && std::semiregular<graph_edge_descriptor_t<G>> &&
    std::semiregular<graph_directed_category_t<G>> && std::semiregular<graph_edge_parallel_category_t<G>> &&
    std::semiregular<graph_traversal_category_t<G>>;

template <typename R, typename G>
concept EdgeRange = std::ranges::input_range<R> && std::convertible_to<std::ranges::range_value_t<R>, graph_edge_descriptor_t<G>>;

template <typename R, typename G>
concept VertexRange = std::ranges::input_range<R> && std::convertible_to<std::ranges::range_value_t<R>, graph_vertex_descriptor_t<G>>;

template <typename G>
concept IncidenceGraph =
    Graph<G> && std::regular<graph_edge_descriptor_t<G>> &&
    !std::is_same_v<graph_degree_size_type_t<G>, void> &&
    std::convertible_to<graph_traversal_category_t<G>, incidence_graph_tag> &&
    requires(const G& g, graph_vertex_descriptor_t<G> u) {
  { out_edges(u, g) } -> EdgeRange<G>;
  { out_degree(u, g) } -> std::integral;
}
&&requires(const G& g, graph_edge_descriptor_t<G> e) {
  { source(e, g) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
  { target(e, g) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
};

template <typename G>
concept BidirectionalGraph = IncidenceGraph<G> &&
    std::convertible_to<graph_traversal_category_t<G>, bidirectional_graph_tag> &&
    requires(const G& g, graph_vertex_descriptor_t<G> u) {
  { in_edges(u, g) } -> EdgeRange<G>;
  { in_degree(u, g) } -> std::integral;
  { degree(u, g) } -> std::integral;
};

template <typename G>
concept AdjacencyGraph = Graph<G> &&
    std::convertible_to<graph_traversal_category_t<G>, adjacency_graph_tag> &&
    requires(const G& g, graph_vertex_descriptor_t<G> u) {
  { adjacent_vertices(u, g) } -> VertexRange<G>;
};

template <typename G>
concept VertexListGraph = Graph<G> &&
    std::integral<graph_vertices_size_type_t<G>> &&
    std::convertible_to<graph_traversal_category_t<G>, vertex_list_graph_tag> &&
    requires(const G& g) {
  { vertices(g) } -> VertexRange<G>;
};

template <typename G>
concept EdgeListGraph = Graph<G> && std::regular<graph_edge_descriptor_t<G>> &&
    std::integral<graph_edges_size_type_t<G>> &&
    std::convertible_to<graph_traversal_category_t<G>, edge_list_graph_tag> &&
    requires(const G& g) {
  { edges(g) } -> EdgeRange<G>;
  { num_edges(g) } -> std::integral;
} && requires(const G& g, graph_edge_descriptor_t<G> e) {
  { source(e, g) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
  { target(e, g) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
};

template <typename G>
concept VertexAndEdgeListGraph = VertexListGraph<G> && EdgeListGraph<G>;

template <typename G>
concept EdgeMutableGraph = Graph<G> && std::regular<graph_edge_descriptor_t<G>> &&
    requires(G& g, graph_vertex_descriptor_t<G> u, graph_edge_descriptor_t<G> e) {
  std::tie(e, std::declval<bool&>()) = add_edge(u, u, g);
  remove_edge(u, u, g);
  remove_edge(e, g);
  clear_vertex(u, g);
};

template <typename G>
concept VertexMutableGraph = Graph<G> && requires(G& g, graph_vertex_descriptor_t<G> u) {
  { add_vertex(g) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
  remove_vertex(u, g);
};

template <typename G>
concept MutableGraph = EdgeMutableGraph<G> && VertexMutableGraph<G>;

template <class edge_descriptor>
struct dummy_edge_predicate {
  bool operator()(const edge_descriptor&) const { return false; }
};

template <typename G>
concept MutableIncidenceGraph = MutableGraph<G> &&
    requires(G& g, graph_vertex_descriptor_t<G> u) {
  remove_out_edge_if(
      u, [](graph_edge_descriptor_t<G>) { return false; }, g);
};

template <typename G>
concept MutableBidirectionalGraph = MutableIncidenceGraph<G> && requires(G& g, graph_vertex_descriptor_t<G> u) {
  remove_in_edge_if(
      u, [](graph_edge_descriptor_t<G>) { return false; }, g);
};

template <typename G>
concept MutableEdgeListGraph = EdgeMutableGraph<G> && requires(G& g) {
  remove_edge_if([](graph_edge_descriptor_t<G>) { return false; }, g);
};

template <typename G>
concept VertexMutablePropertyGraph = VertexMutableGraph<G> && requires(G& g, graph_vertex_descriptor_t<G> v, const vertex_property_type<G>& vp) {
  { add_vertex(g, vp) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
  { add_vertex(g, std::move(vp)) } -> std::convertible_to<graph_vertex_descriptor_t<G>>;
  remove_vertex(v, g, &vp);
};

template <typename G>
concept EdgeMutablePropertyGraph = EdgeMutableGraph<G> &&
    requires(G& g, graph_vertex_descriptor_t<G> u, graph_edge_descriptor_t<G> e, const edge_property_type<G>& ep) {
  { add_edge(u, u, g, ep) } -> std::convertible_to<std::pair<graph_edge_descriptor_t<G>, bool>>;
  { add_edge(u, u, g, std::move(ep)) } -> std::convertible_to<std::pair<graph_edge_descriptor_t<G>, bool>>;
  remove_edge(u, u, g, &ep);
  remove_edge(e, g, &ep);
};

template <typename G>
concept MutablePropertyGraph = EdgeMutablePropertyGraph<G> && VertexMutablePropertyGraph<G>;

template <typename G>
concept AdjacencyMatrix = Graph<G> && std::regular<graph_edge_descriptor_t<G>> &&
    requires(const G& g, graph_vertex_descriptor_t<G> u) {
  { edge(u, u, g) } -> std::convertible_to<std::pair<graph_edge_descriptor_t<G>, bool>>;
};

template <typename G, typename Key, typename Property>
concept ReadablePropertyGraph = Graph<G> && ReadablePropertyMap<property_map_const_t<G, Property>, Key> &&
    requires(const G& g, const Key& k) {
  { get(Property{}, g) } -> std::convertible_to<property_map_const_t<G, Property>>;
  { get(Property{}, g, k) } -> std::convertible_to<property_traits_value_t<property_map_const_t<G, Property>>>;
};

template <typename G, typename Key, typename Property>
concept PropertyGraph = ReadablePropertyGraph<G, Key, Property> &&
    ReadWritePropertyMap<property_map_t<G, Property>, Key> &&
    requires(G& g, const Key& k, property_traits_value_t<property_map_t<G, Property>>&& pval) {
  { get(Property{}, g) } -> std::convertible_to<property_map_t<G, Property>>;
  { get(Property{}, g, k) } -> std::convertible_to<property_traits_value_t<property_map_t<G, Property>>>;
  put(Property{}, g, k, std::move(pval));
};

template <typename G, typename Key, typename Property>
concept LvaluePropertyGraph = ReadablePropertyGraph<G, Key, Property> &&
    LvaluePropertyMap<property_map_const_t<G, Property>, Key> &&
    requires(G& g, const Key& k, property_traits_value_t<property_map_t<G, Property>>&& pval) {
  { get(Property{}, g, k) } -> std::convertible_to<property_traits_value_t<property_map_t<G, Property>>>;
  put(Property{}, g, k, std::move(pval));
};

// The *IndexGraph concepts are "semantic" graph concpepts. These can be
// applied to describe any graph that has an index map that can be accessed
// using the get(*_index, g) method. For example, adjacency lists with
// VertexSet == vecS are implicitly models of this concept.

template <typename G>
concept VertexIndexGraph = Graph<G> &&
    ReadablePropertyMap<property_map_const_t<G, vertex_index_t>, graph_vertex_descriptor_t<G>> && requires(G& g) {
  renumber_vertex_indices(g);
};

template <typename G>
concept EdgeIndexGraph = Graph<G> && std::regular<graph_edge_descriptor_t<G>> &&
    ReadablePropertyMap<property_map_const_t<G, edge_index_t>, graph_edge_descriptor_t<G>> && requires(G& g) {
  renumber_edge_indices(g);
};

template <typename T>
concept ColorValue = std::regular<T> && requires {
  { color_traits<T>::white() } -> std::convertible_to<T>;
  { color_traits<T>::gray() } -> std::convertible_to<T>;
  { color_traits<T>::black() } -> std::convertible_to<T>;
};

template <typename T, typename V>
concept BasicMatrix = requires(const T& m, int i) {
  { m[i][i] } -> std::convertible_to<V>;
}
&&requires(T& m, int i) {
  { m[i][i] } -> std::convertible_to<V&>;
};

// The following concepts describe aspects of numeric values and measure
// functions. We're extending the notion of numeric values to include
// emulation for zero and infinity.
template <typename T>
concept NumericValue = std::regular<T> && requires {
  { numeric_values<T>::zero() } -> std::convertible_to<T>;
  { numeric_values<T>::infinity() } -> std::convertible_to<T>;
};

namespace graph_concepts_detail {

BAGL_GRAPH_HAS_TRAIT_MEMBER(degree_type, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(vertex_type, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(distance_type, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(result_type, void)

}  // namespace graph_concepts_detail

template <typename M, typename G>
concept DegreeMeasure = requires(const M& m, const G& g, graph_concepts_detail::get_vertex_type_or_not<M> v) {
  { m(v, g) } -> std::convertible_to<graph_concepts_detail::get_degree_type_or_not<M>>;
};

template <typename M, typename G>
concept DistanceMeasure = requires(const M& m, const G& g, graph_concepts_detail::get_distance_type_or_not<M> d) {
  { m(d, g) } -> std::convertible_to<graph_concepts_detail::get_result_type_or_not<M>>;
};

}  // namespace bagl::concepts

#endif  // BAGL_BAGL_GRAPH_CONCEPTS_H_
