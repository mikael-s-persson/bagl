// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_FILTERED_GRAPH_H_
#define BAGL_BAGL_FILTERED_GRAPH_H_

#include <ranges>

#include "bagl/adjacency_range.h"
#include "bagl/bfl_d_ary_tree.h"
#include "bagl/detail/set_adaptor.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

//=========================================================================
// Some predicate classes.

struct keep_all {
  template <typename T>
  [[nodiscard]] bool operator()(const T&) const {
    return true;
  }
};

// Keep residual edges (used in maximum-flow algorithms).
template <typename ResidualCapacityEdgeMap>
struct is_residual_edge {
  is_residual_edge() = default;
  explicit is_residual_edge(ResidualCapacityEdgeMap rcap) : rcap_(rcap) {}
  template <typename Edge>
  requires concepts::ReadablePropertyMap<ResidualCapacityEdgeMap, Edge>
  [[nodiscard]] bool operator()(const Edge& e) const { return 0 < get(rcap_, e); }
  ResidualCapacityEdgeMap rcap_;
};

template <typename Set>
struct is_in_subset {
  is_in_subset() = default;
  explicit is_in_subset(const Set& s) : s_(&s) {}

  template <typename Elt>
  [[nodiscard]] bool operator()(const Elt& x) const {
    return set_contains(*s_, x);
  }
  const Set* s_ = nullptr;
};

template <typename Set>
struct is_not_in_subset {
  is_not_in_subset() = default;
  explicit is_not_in_subset(const Set& s) : s_(&s) {}

  template <typename Elt>
  [[nodiscard]] bool operator()(const Elt& x) const {
    return !set_contains(*s_, x);
  }
  const Set* s_ = nullptr;
};

namespace filtered_detail {

template <typename EdgePredicate, typename VertexPredicate, concepts::IncidenceGraph G>
struct out_edge_predicate {
  out_edge_predicate() = default;
  out_edge_predicate(EdgePredicate ep, VertexPredicate vp, const G& g) : edge_pred_(ep), vertex_pred_(vp), g_(&g) {}

  [[nodiscard]] bool operator()(graph_edge_descriptor_t<G> e) const {
    return edge_pred_(e) && vertex_pred_(target(e, *g_));
  }
  EdgePredicate edge_pred_;
  VertexPredicate vertex_pred_;
  const G* g_ = nullptr;
};

template <typename EdgePredicate, typename VertexPredicate, concepts::IncidenceGraph G>
struct in_edge_predicate {
  in_edge_predicate() = default;
  in_edge_predicate(EdgePredicate ep, VertexPredicate vp, const G& g) : edge_pred_(ep), vertex_pred_(vp), g_(&g) {}

  [[nodiscard]] bool operator()(graph_edge_descriptor_t<G> e) const {
    return edge_pred_(e) && vertex_pred_(source(e, *g_));
  }
  EdgePredicate edge_pred_;
  VertexPredicate vertex_pred_;
  const G* g_ = nullptr;
};

template <typename EdgePredicate, typename VertexPredicate, concepts::IncidenceGraph G>
struct edge_predicate {
  edge_predicate() = default;
  edge_predicate(EdgePredicate ep, VertexPredicate vp, const G& g) : edge_pred_(ep), vertex_pred_(vp), g_(&g) {}

  [[nodiscard]] bool operator()(graph_edge_descriptor_t<G> e) const {
    return edge_pred_(e) && vertex_pred_(source(e, *g_)) && vertex_pred_(target(e, *g_));
  }
  EdgePredicate edge_pred_;
  VertexPredicate vertex_pred_;
  const G* g_ = nullptr;
};

}  // namespace filtered_detail

//===========================================================================
// Filtered Graph

struct filtered_graph_tag {};

// This base class is a stupid hack to change overload resolution
// rules for the source and target functions so that they are a
// worse match than the source and target functions defined for
// pairs in graph_traits.hpp. I feel dirty. -JGS
template <class G>
struct filtered_graph_base {
  using vertex_descriptor = graph_vertex_descriptor_t<G>;
  using edge_descriptor = graph_edge_descriptor_t<G>;
  explicit filtered_graph_base(G& g) : g_(g) {}
  // protected:
  G& g_;
};

template <typename Graph, typename EdgePredicate, typename VertexPredicate = keep_all>
class filtered_graph : public filtered_graph_base<Graph> {
  using Base = filtered_graph_base<Graph>;
  using self = filtered_graph;

 public:
  using graph_type = Graph;
  using OutEdgePred = filtered_detail::out_edge_predicate<EdgePredicate, VertexPredicate, self>;
  using InEdgePred = filtered_detail::in_edge_predicate<EdgePredicate, VertexPredicate, self>;
  using EdgePred = filtered_detail::edge_predicate<EdgePredicate, VertexPredicate, self>;

  // Constructors
  filtered_graph(Graph& g, EdgePredicate ep) : Base(g), edge_pred_(ep) {}

  filtered_graph(Graph& g, EdgePredicate ep, VertexPredicate vp) : Base(g), edge_pred_(ep), vertex_pred_(vp) {}

  // Graph requirements
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using edge_descriptor = graph_edge_descriptor_t<Graph>;
  using directed_category = graph_directed_category_t<Graph>;
  using edge_parallel_category = graph_edge_parallel_category_t<Graph>;
  using traversal_category = graph_traversal_category_t<Graph>;

  // IncidenceGraph requirements
  auto out_edges(vertex_descriptor u) const {
    return out_edges(u, this->g_) | std::views::filter(OutEdgePred{edge_pred_, vertex_pred_, *this});
  }
  using out_edge_range = decltype(std::declval<self>().out_edges(std::declval<vertex_descriptor>()));
  using degree_size_type = graph_degree_size_type_t<Graph>;

  // AdjacencyGraph requirements
  using adjacency_range = decltype(adjacency_range(std::declval<out_edge_range>(), std::declval<self>()));

  // BidirectionalGraph requirements
  auto in_edges(vertex_descriptor u) const {
    return in_edges(u, this->g_) | std::views::filter(InEdgePred{edge_pred_, vertex_pred_, *this});
  }
  using in_edge_range = decltype(std::declval<self>().in_edges(std::declval<vertex_descriptor>()));

  // VertexListGraph requirements
  auto vertices() const { return vertices(this->g_) | std::views::filter(vertex_pred_); }
  using vertex_range = decltype(std::declval<self>().vertices());
  using vertices_size_type = graph_vertices_size_type_t<Graph>;

  // EdgeListGraph requirements
  auto edges() const { return edges(this->g_) | std::views::filter(EdgePred{edge_pred_, vertex_pred_, *this}); }
  using edge_range = decltype(std::declval<self>().edges());
  using edges_size_type = graph_edges_size_type_t<Graph>;

  using graph_tag = filtered_graph_tag;

  using vertex_property_type = vertex_property_type<Graph>;
  using edge_property_type = edge_property_type<Graph>;
  using graph_property_type = graph_property_type<Graph>;

  using graph_bundled = graph_bundle_type<Graph>;
  using vertex_bundled = vertex_bundle_type<Graph>;
  using edge_bundled = edge_bundle_type<Graph>;

  // Bundled properties support
  template <typename Descriptor>
  typename graph_traits_detail::bundled_result<Graph, Descriptor>::type& operator[](Descriptor x) {
    return this->g_[x];
  }

  template <typename Descriptor>
  typename graph_traits_detail::bundled_result<Graph, Descriptor>::type const& operator[](Descriptor x) const {
    return this->g_[x];
  }

  static vertex_descriptor null_vertex() { return graph_traits<Graph>::null_vertex(); }

  // private:
  EdgePredicate edge_pred_;
  VertexPredicate vertex_pred_;
};

//===========================================================================
// Non-member functions for the Filtered Edge Graph

// Helper functions
template <typename Graph, typename EdgePredicate>
auto make_filtered_graph(Graph& g, EdgePredicate ep) {
  return filtered_graph<Graph, EdgePredicate>(g, ep);
}
template <typename Graph, typename EdgePredicate, typename VertexPredicate>
auto make_filtered_graph(Graph& g, EdgePredicate ep, VertexPredicate vp) {
  return filtered_graph<Graph, EdgePredicate, VertexPredicate>(g, ep, vp);
}

template <typename Graph, typename EdgePredicate>
auto make_filtered_graph(const Graph& g, EdgePredicate ep) {
  return filtered_graph<const Graph, EdgePredicate>(g, ep);
}
template <typename Graph, typename EdgePredicate, typename VertexPredicate>
auto make_filtered_graph(const Graph& g, EdgePredicate ep, VertexPredicate vp) {
  return filtered_graph<const Graph, EdgePredicate, VertexPredicate>(g, ep, vp);
}

template <typename G, typename EP, typename VP>
auto vertices(const filtered_graph<G, EP, VP>& g) {
  return g.vertices();
}

template <typename G, typename EP, typename VP>
auto edges(const filtered_graph<G, EP, VP>& g) {
  return g.edges();
}

// An alternative for num_vertices() and num_edges() would be to
// count the number in the filtered graph. This is problematic
// because of the interaction with the vertex indices...  they would
// no longer go from 0 to num_vertices(), which would cause trouble
// for algorithms allocating property storage in an array. We could
// try to create a mapping to new recalibrated indices, but I don't
// see an efficient way to do this.
//
// However, the current solution is still unsatisfactory because
// the following semantic constraints no longer hold:
// boost::tie(vi, viend) = vertices(g);
// assert(std::distance(vi, viend) == num_vertices(g));

template <typename G, typename EP, typename VP>
auto num_vertices(const filtered_graph<G, EP, VP>& g) {
  return num_vertices(g.g_);
}

template <typename G, typename EP, typename VP>
auto num_edges(const filtered_graph<G, EP, VP>& g) {
  return num_edges(g.g_);
}

template <typename G>
auto source(graph_edge_descriptor_t<G> e, const filtered_graph_base<G>& g) {
  return source(e, g.g_);
}

template <typename G>
auto target(graph_edge_descriptor_t<G> e, const filtered_graph_base<G>& g) {
  return target(e, g.g_);
}

template <typename G, typename EP, typename VP>
auto out_edges(graph_vertex_descriptor_t<G> u, const filtered_graph<G, EP, VP>& g) {
  return g.out_edges(u);
}

template <typename G, typename EP, typename VP>
auto out_degree(graph_vertex_descriptor_t<G> u, const filtered_graph<G, EP, VP>& g) {
  return std::ranges::distance(out_edges(u, g));
}

template <typename G, typename EP, typename VP>
auto adjacent_vertices(graph_vertex_descriptor_t<G> u, const filtered_graph<G, EP, VP>& g) {
  return adjacency_range(out_edges(u, g), g);
}

template <typename G, typename EP, typename VP>
auto in_edges(graph_vertex_descriptor_t<G> u, const filtered_graph<G, EP, VP>& g) {
  return g.in_edges(u);
}

template <typename G, typename EP, typename VP>
auto in_degree(graph_vertex_descriptor_t<G> u, const filtered_graph<G, EP, VP>& g) {
  return std::ranges::distance(in_edges(u, g));
}

template <typename G, typename EP, typename VP>
auto degree(graph_vertex_descriptor_t<G> u, const filtered_graph<G, EP, VP>& g) {
  if constexpr (is_undirected_graph_v<G>) {
    return out_degree(u, g);
  } else {
    return out_degree(u, g) + in_degree(u, g);
  }
}

template <typename G, typename EP, typename VP>
auto edge(graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v, const filtered_graph<G, EP, VP>& g) {
  auto [e, exists] = edge(u, v, g.g_);
  return std::pair{e, exists && g.edge_pred_(e)};
}

//===========================================================================
// Property map

template <typename G, typename EP, typename VP, typename Property>
struct property_map<filtered_graph<G, EP, VP>, Property> : property_map<G, Property> {};

template <typename G, typename EP, typename VP, typename Property>
auto get(Property p, filtered_graph<G, EP, VP>& g) {
  return get(p, g.g_);
}

template <typename G, typename EP, typename VP, typename Property>
auto get(Property p, const filtered_graph<G, EP, VP>& g) {
  return get(p, g.g_);
}

template <typename G, typename EP, typename VP, typename Property, typename Key>
auto get(Property p, const filtered_graph<G, EP, VP>& g, const Key& k) {
  return get(p, g.g_, k);
}

template <typename G, typename EP, typename VP, typename Property, typename Key, typename Value>
void put(Property p, filtered_graph<G, EP, VP>& g, const Key& k, const Value& val) {
  put(p, g.g_, k, val);
}

//===========================================================================
// Some filtered subgraph specializations

template <typename Graph, typename Set>
using vertex_subset_filter_t = filtered_graph<Graph, keep_all, is_in_subset<Set>>;
template <typename Graph, typename Set>
auto make_vertex_subset_filter(Graph& g, const Set& s) {
  return vertex_subset_filter_t<Graph, Set>(g, keep_all(), is_in_subset{s});
}

// This is misspelled, but present for backwards compatibility; new code
// should use the version below that has the correct spelling.
template <typename Graph, typename Set>
using vertex_subset_complement_filter_t = filtered_graph<Graph, keep_all, is_not_in_subset<Set>>;
template <typename Graph, typename Set>
auto make_vertex_subset_complement_filter(Graph& g, const Set& s) {
  return vertex_subset_complement_filter_t<Graph, Set>(g, keep_all(), is_not_in_subset{s});
}

// Filter that uses a property map whose value_type is a boolean
template <typename PropertyMap>
struct property_map_filter {
  property_map_filter() = default;
  explicit property_map_filter(const PropertyMap& property_map) : property_map_(property_map) {}

  template <typename Key>
  [[nodiscard]] bool operator()(const Key& key) const {
    return get(property_map_, key);
  }

 private:
  PropertyMap property_map_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_FILTERED_GRAPH_H_
