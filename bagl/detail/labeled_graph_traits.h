// Copyright (C) 2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_LABELED_GRAPH_TRAITS_H_
#define BAGL_BAGL_DETAIL_LABELED_GRAPH_TRAITS_H_

#include "bagl/graph_mutability_traits.h"

namespace bagl {

// Extend the graph mutability traits (and metafunctions) to include options
// for labeled graphs.

// NOTE: the label_vertex tag denotes the fact that you can basically assign
// arbitrary labels to vertices without modifying the actual graph.

// TODO: We might also overlay the uniqueness/multiplicity of labels in this
// hierarchy also. For now, we just assumed that labels are unique.

struct label_vertex_tag {};
struct labeled_add_vertex_tag : virtual label_vertex_tag {};
struct labeled_add_vertex_property_tag : virtual labeled_add_vertex_tag {};
struct labeled_remove_vertex_tag {};
struct labeled_add_edge_tag : virtual label_vertex_tag {};
struct labeled_add_edge_property_tag : virtual labeled_add_edge_tag {};
struct labeled_remove_edge_tag {};

struct labeled_mutable_vertex_graph_tag : virtual labeled_add_vertex_tag, virtual labeled_remove_vertex_tag {};
struct labeled_mutable_vertex_property_graph_tag : virtual labeled_add_vertex_property_tag,
                                                   virtual labeled_remove_vertex_tag {};
struct labeled_mutable_edge_graph_tag : virtual labeled_add_edge_tag, virtual labeled_remove_edge_tag {};
struct labeled_mutable_edge_property_graph_tag : virtual labeled_add_edge_property_tag,
                                                 virtual labeled_remove_edge_tag {};

struct labeled_graph_tag : virtual label_vertex_tag {};
struct labeled_mutable_graph_tag : virtual labeled_mutable_vertex_graph_tag, virtual labeled_mutable_edge_graph_tag {};
struct labeled_mutable_property_graph_tag : virtual labeled_mutable_vertex_property_graph_tag,
                                            virtual labeled_mutable_edge_property_graph_tag {};
struct labeled_add_only_property_graph_tag : virtual labeled_add_vertex_property_tag,
                                             virtual labeled_mutable_edge_property_graph_tag {};

// Metafunctions

template <typename Graph>
constexpr bool graph_has_add_vertex_by_label_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_add_vertex_tag>;

template <typename Graph>
constexpr bool graph_has_add_vertex_by_label_with_property_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_add_vertex_property_tag>;

template <typename Graph>
constexpr bool graph_has_remove_vertex_by_label_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_remove_vertex_tag>;

template <typename Graph>
constexpr bool graph_has_add_edge_by_label_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_add_edge_tag>;

template <typename Graph>
constexpr bool graph_has_add_edge_by_label_with_property_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_add_edge_property_tag>;

template <typename Graph>
constexpr bool graph_has_remove_edge_by_label_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_remove_edge_tag>;

template <typename Graph>
constexpr bool is_labeled_mutable_vertex_graph_v =
    graph_has_add_vertex_by_label_v<Graph>&& graph_has_remove_vertex_by_label_v<Graph>;

template <typename Graph>
constexpr bool is_labeled_mutable_vertex_property_graph_v =
    graph_has_add_vertex_by_label_v<Graph>&& graph_has_remove_vertex_by_label_v<Graph>;

template <typename Graph>
constexpr bool is_labeled_mutable_edge_graph_v =
    graph_has_add_edge_by_label_v<Graph>&& graph_has_remove_edge_by_label_v<Graph>;

template <typename Graph>
constexpr bool is_labeled_mutable_edge_property_graph_v =
    graph_has_add_edge_by_label_v<Graph>&& graph_has_remove_edge_by_label_v<Graph>;

template <typename Graph>
constexpr bool is_labeled_mutable_graph_v =
    is_labeled_mutable_vertex_graph_v<Graph>&& is_labeled_mutable_edge_graph_v<Graph>;

template <typename Graph>
constexpr bool is_labeled_mutable_property_graph_v =
    is_labeled_mutable_vertex_property_graph_v<Graph>&& is_labeled_mutable_edge_property_graph_v<Graph>;

template <typename Graph>
constexpr bool is_labeled_add_only_property_graph_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, labeled_add_only_property_graph_tag>;

template <typename Graph>
constexpr bool is_labeled_graph_v =
    std::is_convertible_v<typename graph_mutability_traits<Graph>::category, label_vertex_tag>;

template <typename>
struct graph_mutability_traits;

namespace graph_detail {
// The determine mutability metafunction computes a labeled mutability tag
// based on the mutability of the given graph type. This is used by the
// graph_mutability_traits specialization below.
template <typename Graph>
struct determine_mutability {
  using type = std::conditional_t<
      is_add_only_property_graph_v<Graph>, labeled_add_only_property_graph_tag,
      std::conditional_t<is_mutable_property_graph_v<Graph>, labeled_mutable_property_graph_tag,
                         std::conditional_t<is_mutable_graph_v<Graph>, labeled_mutable_graph_tag,
                                            std::conditional_t<is_mutable_edge_graph_v<Graph>, labeled_graph_tag,
                                                               typename graph_mutability_traits<Graph>::category>>>>;
};
}  // namespace graph_detail

// Specialize mutability traits for the labeled graph.
// This specialization depends on the mutability of the underlying graph type.
// If the underlying graph is fully mutable, this is also fully mutable.
// Otherwise, it's different.
template <typename G, typename L, typename S>
struct graph_mutability_traits<labeled_graph<G, L, S>> {
  using category = typename graph_detail::determine_mutability<typename labeled_graph<G, L, S>::graph_type>::type;
};

}  // namespace bagl

#endif  // BAGL_BAGL_DETAIL_LABELED_GRAPH_TRAITS_H_
