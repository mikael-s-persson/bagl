// Copyright 2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_MUTABILITY_TRAITS_H_
#define BAGL_BAGL_GRAPH_MUTABILITY_TRAITS_H_

#include <type_traits>

namespace bagl {

// The mutabiltiy categories classify graphs by their mutating operations
// on the edge and vertex sets. This is a substantially more refined
// categorization than the MutableGraph and MutablePropertyGraph denote.
// Currently, this framework is only used in the graph tests to help
// dispatch test to the correct places. However, there are probably some
// constructive or destructive algorithms (i.e., graph generators) that
// may use these to describe requirements on graph inputs.

struct add_vertex_tag {};
struct add_vertex_property_tag : virtual add_vertex_tag {};
struct add_edge_tag {};
struct add_edge_property_tag : virtual add_edge_tag {};
struct remove_vertex_tag {};
struct remove_edge_tag {};

struct mutable_vertex_graph_tag : virtual add_vertex_tag, virtual remove_vertex_tag {};
struct mutable_vertex_property_graph_tag : virtual add_vertex_property_tag, virtual remove_vertex_tag {};

struct mutable_edge_graph_tag : virtual add_edge_tag, virtual remove_edge_tag {};
struct mutable_edge_property_graph_tag : virtual add_edge_property_tag, virtual remove_edge_tag {};

struct mutable_graph_tag : virtual mutable_vertex_graph_tag, virtual mutable_edge_graph_tag {};
struct mutable_property_graph_tag : virtual mutable_vertex_property_graph_tag,
                                    virtual mutable_edge_property_graph_tag {};

// Some graphs just don't like to be torn down. Note this only restricts
// teardown to the set of vertices, not the vertex set.
// TODO: Find a better name for this tag.
struct add_only_property_graph_tag : virtual add_vertex_property_tag, virtual mutable_edge_property_graph_tag {};

/**
 * The graph_mutability_traits provide methods for determining the
 * interfaces supported by graph classes for adding and removing vertices
 * and edges.
 */
template <typename Graph>
struct graph_mutability_traits {
  using category = typename Graph::mutability_category;
};
template <typename Graph>
using graph_mutability_category_t = typename graph_mutability_traits<Graph>::category;

template <typename Graph>
constexpr bool graph_has_add_vertex_v = std::is_convertible_v<graph_mutability_category_t<Graph>, add_vertex_tag>;

template <typename Graph>
constexpr bool graph_has_add_vertex_with_property_v =
    std::is_convertible_v<graph_mutability_category_t<Graph>, add_vertex_property_tag>;

template <typename Graph>
constexpr bool graph_has_remove_vertex_v = std::is_convertible_v<graph_mutability_category_t<Graph>, remove_vertex_tag>;

template <typename Graph>
constexpr bool graph_has_add_edge_v = std::is_convertible_v<graph_mutability_category_t<Graph>, add_edge_tag>;

template <typename Graph>
constexpr bool graph_has_add_edge_with_property_v =
    std::is_convertible_v<graph_mutability_category_t<Graph>, add_edge_property_tag>;

template <typename Graph>
constexpr bool graph_has_remove_edge_v = std::is_convertible_v<graph_mutability_category_t<Graph>, remove_edge_tag>;

template <typename Graph>
constexpr bool is_mutable_vertex_graph_v = graph_has_add_vertex_v<Graph>&& graph_has_remove_vertex_v<Graph>;

template <typename Graph>
constexpr bool is_mutable_vertex_property_graph_v =
    graph_has_add_vertex_with_property_v<Graph>&& graph_has_remove_vertex_v<Graph>;

template <typename Graph>
constexpr bool is_mutable_edge_graph_v = graph_has_add_edge_v<Graph>&& graph_has_remove_edge_v<Graph>;

template <typename Graph>
constexpr bool is_mutable_edge_property_graph_v =
    graph_has_add_edge_with_property_v<Graph>&& graph_has_remove_edge_v<Graph>;

template <typename Graph>
constexpr bool is_mutable_graph_v = is_mutable_vertex_graph_v<Graph>&& is_mutable_edge_graph_v<Graph>;

template <typename Graph>
constexpr bool is_mutable_property_graph_v =
    is_mutable_vertex_property_graph_v<Graph>&& is_mutable_edge_property_graph_v<Graph>;

template <typename Graph>
constexpr bool is_add_only_property_graph_v =
    std::is_convertible_v<graph_mutability_category_t<Graph>, add_only_property_graph_tag>;

}  // namespace bagl

#endif  // BAGL_BAGL_GRAPH_MUTABILITY_TRAITS_H_
