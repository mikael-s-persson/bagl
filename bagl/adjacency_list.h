// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_ADJACENCY_LIST_H_
#define BAGL_BAGL_ADJACENCY_LIST_H_

#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <utility>

#include "bagl/adjacency_range.h"
#include "bagl/detail/adjlist_containers.h"
#include "bagl/detail/adjlist_ranges.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_mutability_traits.h"
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"
#include "bagl/tree_traits.h"
#include "bagl/zip_range.h"

namespace bagl {

template <class OutEdgeListS = vec_s, class VertexListS = vec_s, class DirectedS = directed_s>
struct adjacency_list_traits {
  using is_rand_access = typename container_detail::is_random_access<VertexListS>::type;
  using is_bidir = typename DirectedS::is_bidir_t;
  using is_directed = typename DirectedS::is_directed_t;

  using directed_category = std::conditional_t<is_bidir::value, bidirectional_tag,
                                               std::conditional_t<is_directed::value, directed_tag, undirected_tag>>;

  using edge_parallel_category = typename container_detail::parallel_edge_traits<OutEdgeListS>::type;

  using traversal_category = adjlist_detail::adjlist_traversal_tag<DirectedS>;

  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;

  // Best effort attempt at producing descriptors without knowing properties types.
  // Do not use this unless you absolutely need to.
  using vertex_descriptor_or_void =
      std::conditional_t<(std::is_same_v<VertexListS, vec_s> || std::is_same_v<VertexListS, pool_s>), std::size_t,
                         void>;
  using edge_descriptor_or_void =
      std::conditional_t<!std::is_same_v<vertex_descriptor_or_void, void> &&
                             (std::is_same_v<OutEdgeListS, vec_s> || std::is_same_v<OutEdgeListS, pool_s>),
                         container_detail::edge_desc<vertex_descriptor_or_void, std::size_t>, void>;
};

struct adj_list_tag {};

template <typename VertexListS>
struct adjacency_list_disallowed_vertex_list {
  using type = void;
};

template <>
struct adjacency_list_disallowed_vertex_list<set_s> {};
template <>
struct adjacency_list_disallowed_vertex_list<multiset_s> {};
template <>
struct adjacency_list_disallowed_vertex_list<unordered_set_s> {};
template <>
struct adjacency_list_disallowed_vertex_list<unordered_multiset_s> {};

// forward-declare:
template <typename OutEdgeListS = vec_s, typename VertexListS = vec_s, typename DirectedS = directed_s,
          typename VertexProperties = no_property, typename EdgeProperties = no_property,
          typename GraphProperties = no_property>
class adjacency_list;

template <typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename GraphProperties>
void do_graph_deep_copy(
    adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& lhs,
    const adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& rhs);

/**
 * This class implements an adjacency-list based on containers that are tailored
 * to store elements of a graph.
 * \tparam OutEdgeListS A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexListS A type tag to choose the storage policy for the vertices.
 * \tparam DirectedS A type tag to choose the directionality of the edges.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties, typename GraphProperties>
class adjacency_list {
 public:
  using check_allowed_vertex_list = typename adjacency_list_disallowed_vertex_list<VertexListS>::type;

  using self = adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties, GraphProperties>;

  using out_edge_list_selector = OutEdgeListS;
  using vertex_list_selector = VertexListS;
  using directed_selector = DirectedS;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;
  using graph_property_type = GraphProperties;

  using vertex_bundled = lookup_one_property_t<VertexProperties, vertex_bundle_t>;
  using edge_bundled = lookup_one_property_t<EdgeProperties, edge_bundle_t>;
  using graph_bundled = lookup_one_property_t<GraphProperties, graph_bundle_t>;

  using storage_type = typename adjlist_detail::adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS,
                                                                         VertexProperties, EdgeProperties>;

  using vertex_descriptor = typename storage_type::vertex_descriptor;
  using vertices_size_type = typename storage_type::vertices_size_type;

  using edge_descriptor = typename storage_type::edge_descriptor;
  using edges_size_type = typename storage_type::edges_size_type;
  using degree_size_type = edges_size_type;

  using vertex_stored_impl = typename storage_type::vertex_stored_type;
  using vertex_value_impl = typename storage_type::vertex_value_type;

  using edge_stored_impl = typename storage_type::edge_stored_type;
  using edge_value_impl = typename storage_type::edge_value_type;

  using Traits = adjacency_list_traits<OutEdgeListS, VertexListS, DirectedS>;

  using directed_category = typename Traits::directed_category;
  using edge_parallel_category = typename Traits::edge_parallel_category;
  using traversal_category = typename Traits::traversal_category;

  using mutability_category = mutable_property_graph_tag;

  using graph_tag = adj_list_tag;

  /**
   * This static member function outputs the null-vertex (invalid vertex descriptor).
   * \return A null-vertex descriptor (invalid vertex descriptor).
   */
  static vertex_descriptor null_vertex() { return container_detail::null_desc<vertex_descriptor>::value(); }

  /**
   * This static member function outputs the null-edge (invalid edge descriptor).
   * \return A null-edge descriptor (invalid edge descriptor).
   */
  static edge_descriptor null_edge() { return container_detail::null_desc<edge_descriptor>::value(); }

  // private:
  storage_type m_pack;
  graph_property_type m_graph_prop;

  /**
   * Constructs an empty adjacency-list.
   */
  adjacency_list() : m_pack(), m_graph_prop() {}

  explicit adjacency_list(std::size_t num_vertices, graph_property_type graph_prop = {})
      : m_pack(num_vertices), m_graph_prop(std::move(graph_prop)) {}

  ~adjacency_list() = default;

  adjacency_list(const self& rhs) : m_pack(), m_graph_prop(rhs.m_graph_prop) { do_graph_deep_copy(*this, rhs); }
  self& operator=(const self& rhs) {
    if (this != &rhs) {
      m_graph_prop = rhs.m_graph_prop;
      do_graph_deep_copy(*this, rhs);
    }
    return *this;
  }

  adjacency_list(self&& rhs) noexcept : m_pack(std::move(rhs.m_pack)), m_graph_prop(std::move(rhs.m_graph_prop)) {}
  self& operator=(self&& rhs) noexcept {
    m_pack = std::move(rhs.m_pack);
    m_graph_prop = std::move(rhs.m_graph_prop);
    return *this;
  }

  // Construct from a given number of vertices and an edge range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
  adjacency_list(vertices_size_type num_vertices, const EdgeRange& e_range, graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto& v : tmp_vs) {
      v = m_pack.add_vertex(vertex_property_type{});
    }

    for (auto [u, v] : e_range) {
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], edge_property_type{});
    }
  }

  // Construct from a given number of vertices and an edge and edge-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
      std::constructible_from<edge_property_type, std::ranges::range_reference_t<EdgePropRange>>
      adjacency_list(vertices_size_type num_vertices, const EdgeRange& e_range, const EdgePropRange& ep_range,
                     graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto& v : tmp_vs) {
      v = m_pack.add_vertex(vertex_property_type{});
    }

    for (auto [e, ep] : zip_range_ref(e_range, ep_range)) {
      auto [u, v] = e;
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], ep);
    }
  }

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
      std::constructible_from<vertex_property_type, std::ranges::range_reference_t<VertexPropRange>>
      adjacency_list(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                     graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto [v, vp] : zip_range_ref(tmp_vs, vp_range)) {
      v = m_pack.add_vertex(vp);
    }

    for (auto [u, v] : e_range) {
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], edge_property_type{});
    }
  }

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange,
            std::ranges::input_range EdgePropRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
      std::constructible_from<vertex_property_type, std::ranges::range_reference_t<VertexPropRange>> &&
      std::constructible_from<edge_property_type, std::ranges::range_reference_t<EdgePropRange>>
      adjacency_list(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                     const EdgePropRange& ep_range, graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto [v, vp] : zip_range_ref(tmp_vs, vp_range)) {
      v = m_pack.add_vertex(vp);
    }

    for (auto [e, ep] : zip_range_ref(e_range, ep_range)) {
      auto [u, v] = e;
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], ep);
    }
  }

  /**
   * Swaps the adjacency-list with another.
   */
  void swap(self& rhs) {
    using std::swap;
    swap(m_graph_prop, rhs.m_graph_prop);
    m_pack.swap(rhs.m_pack);
  }

  /**
   * Clears the adjacency-list of all vertices and edges.
   */
  void clear() { m_pack.clear(); }

  // Indexing operator. Returns a reference to the vertex-bundle associated to the given vertex descriptor.
  auto& operator[](vertex_descriptor v) { return get_property_value(m_pack.get_stored_vertex(v).data, vertex_bundle); }
  // Indexing operator. Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
  const auto& operator[](vertex_descriptor v) const {
    return get_property_value(m_pack.get_stored_vertex(v).data, vertex_bundle);
  }

  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  auto& operator[](const edge_descriptor& e) {
    return get_property_value(m_pack.get_stored_edge_property(e), edge_bundle);
  }
  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  const auto& operator[](const edge_descriptor& e) const {
    return get_property_value(m_pack.get_stored_edge_property(e), edge_bundle);
  }

  // Indexing operator. Returns a reference to the graph-bundle associated to the graph.
  auto& operator[](graph_bundle_t /*unused*/) { return get_property_value(m_graph_prop, graph_bundle); }
  // Indexing operator. Returns a const-reference to the graph-bundle associated to the graph.
  const auto& operator[](graph_bundle_t /*unused*/) const { return get_property_value(m_graph_prop, graph_bundle); }

  // Get a reference to the vertex-property associated to the given vertex descriptor.
  auto& get_property(vertex_descriptor v) { return m_pack.get_stored_vertex(v).data; }
  // Get a const-reference to the vertex-property associated to the given vertex descriptor.
  const auto& get_property(vertex_descriptor v) const { return m_pack.get_stored_vertex(v).data; }

  // Get a reference to the edge-property associated to the given edge descriptor.
  auto& get_property(const edge_descriptor& e) { return m_pack.get_stored_edge_property(e); }
  // Get a const-reference to the edge-property associated to the given edge descriptor.
  const auto& get_property(const edge_descriptor& e) const { return m_pack.get_stored_edge_property(e); }

  // Get a reference to the graph-property associated to the graph.
  auto& get_property(graph_all_t /*unused*/) { return m_graph_prop; }
  // Get a const-reference to the graph-property associated to the graph.
  const auto& get_property(graph_all_t /*unused*/) const { return m_graph_prop; }
};

/**
 * This is the tree-storage specifier for an adjacency-list being used as a tree storage.
 */
template <typename OutEdgeListS = vec_s, typename VertexListS = list_s, typename DirectedS = directed_s>
struct adj_list_as_tree_storage {};

template <typename VertexProperty, typename EdgeProperty, typename OutEdgeListS, typename VertexListS,
          typename DirectedS>
struct tree_storage<VertexProperty, EdgeProperty, adj_list_as_tree_storage<OutEdgeListS, VertexListS, DirectedS>> {
  using type = adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty>;
};

template <typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_storage_traits<adj_list_as_tree_storage<OutEdgeListS, VertexListS, DirectedS>>
    : adjacency_list_traits<OutEdgeListS, VertexListS, DirectedS> {};

template <typename VertexProperty, typename EdgeProperty, typename OutEdgeListS, typename VertexListS,
          typename DirectedS>
struct tree_traits<adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty>> {
  using node_descriptor =
      typename adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty>::vertex_descriptor;
};

#define BAGL_ADJACENCY_LIST_ARGS                                                                                       \
  typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, \
      typename GraphProperties
#define BAGL_ADJACENCY_LIST \
  adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties, GraphProperties>

template <BAGL_ADJACENCY_LIST_ARGS>
void swap(BAGL_ADJACENCY_LIST& lhs, BAGL_ADJACENCY_LIST& rhs) {
  lhs.swap(rhs);
}

/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
auto source(const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, const BAGL_ADJACENCY_LIST& /*unused*/) {
  return e.source;
}

template <BAGL_ADJACENCY_LIST_ARGS>
auto target(const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_stored_edge_target(e);
}

template <BAGL_ADJACENCY_LIST_ARGS>
auto out_edges(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.out_edges(v);
}

template <BAGL_ADJACENCY_LIST_ARGS>
std::size_t out_degree(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_out_degree(v);
}

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
auto in_edges(
    std::enable_if_t<!std::is_same_v<DirectedS, directed_s>, typename BAGL_ADJACENCY_LIST::vertex_descriptor> v,
    const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.in_edges(v);
}

template <BAGL_ADJACENCY_LIST_ARGS>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>, std::size_t> in_degree(
    typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_in_degree(v);
}

template <BAGL_ADJACENCY_LIST_ARGS>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>, std::size_t> degree(
    typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
}

/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
auto vertices(const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.vertices();
}

template <BAGL_ADJACENCY_LIST_ARGS>
auto num_vertices(const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.size();
}

/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
auto edges(const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.edges();
}

template <BAGL_ADJACENCY_LIST_ARGS>
auto num_edges(const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.num_edges();
}

/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
auto adjacent_vertices(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return adjacency_range(out_edges(v, g), g);
}

/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>> inv_adjacent_vertices(
    typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return inv_adjacency_range(in_edges(v, g), g);
}

/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS>
auto edge(typename BAGL_ADJACENCY_LIST::vertex_descriptor u, typename BAGL_ADJACENCY_LIST::vertex_descriptor v,
          const BAGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_edge(u, v);
}

/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

// auto add_vertex(BAGL_ADJACENCY_LIST& g); -> variadic overload

template <BAGL_ADJACENCY_LIST_ARGS>
void clear_vertex(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, BAGL_ADJACENCY_LIST& g) {
  g.m_pack.clear_vertex(v);
}

template <BAGL_ADJACENCY_LIST_ARGS>
void remove_vertex(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, BAGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_vertex(v);
}

// auto add_edge(u, v, BAGL_ADJACENCY_LIST& g); -> variadic overload

template <BAGL_ADJACENCY_LIST_ARGS>
void remove_edge(typename BAGL_ADJACENCY_LIST::vertex_descriptor u, typename BAGL_ADJACENCY_LIST::vertex_descriptor v,
                 BAGL_ADJACENCY_LIST& g) {
  auto [e, e_found] = edge(u, v, g);
  while (e_found) {
    g.m_pack.remove_edge(e);
    std::tie(e, e_found) = edge(u, v, g);
  }
}

template <BAGL_ADJACENCY_LIST_ARGS>
void remove_edge(const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, BAGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_edge(e);
}

/*******************************************************************************************
 *                    MutableIncidenceGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS, typename EdgePred>
void remove_out_edge_if(typename BAGL_ADJACENCY_LIST::vertex_descriptor u, EdgePred pred, BAGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_out_edge_if(u, pred);
}

/*******************************************************************************************
 *                    MutableBidirectionalGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS, typename EdgePred>
std::enable_if_t<!std::is_same_v<DirectedS, directed_s>> remove_in_edge_if(
    typename BAGL_ADJACENCY_LIST::vertex_descriptor v, EdgePred pred, BAGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_in_edge_if(v, pred);
}

/*******************************************************************************************
 *                    MutableEdgeListGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS, typename EdgePred>
void remove_edge_if(EdgePred pred, BAGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_edge_if(pred);
}

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_ARGS, typename... VPArgs>
auto add_vertex(BAGL_ADJACENCY_LIST& g, VPArgs&&... vp_args) {
  return g.m_pack.add_vertex(std::forward<VPArgs>(vp_args)...);
}

template <BAGL_ADJACENCY_LIST_ARGS, typename VProp>
void remove_vertex(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, BAGL_ADJACENCY_LIST& g, VProp* vp) {
  if (vp != nullptr) {
    *vp = std::move(g.get_property(v));
  }
  g.m_pack.remove_vertex(v);
}

template <BAGL_ADJACENCY_LIST_ARGS, typename... EPArgs>
auto add_edge(typename BAGL_ADJACENCY_LIST::vertex_descriptor u, typename BAGL_ADJACENCY_LIST::vertex_descriptor v,
              BAGL_ADJACENCY_LIST& g, EPArgs&&... ep_args) {
  return std::pair<typename BAGL_ADJACENCY_LIST::edge_descriptor, bool>{
      g.m_pack.add_edge(u, v, std::forward<EPArgs>(ep_args)...)};
}

template <BAGL_ADJACENCY_LIST_ARGS, typename EProp>
void remove_edge(typename BAGL_ADJACENCY_LIST::vertex_descriptor u, typename BAGL_ADJACENCY_LIST::vertex_descriptor v,
                 BAGL_ADJACENCY_LIST& g, EProp* ep) {
  auto [e, e_found] = edge(u, v, g);
  while (e_found) {
    if (ep != nullptr) {
      *ep = std::move(g.get_property(e));
    }
    g.m_pack.remove_edge(e);
    std::tie(e, e_found) = edge(u, v, g);
  }
}

template <BAGL_ADJACENCY_LIST_ARGS, typename EProp>
void remove_edge(const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, BAGL_ADJACENCY_LIST& g, EProp* ep) {
  if (ep != nullptr) {
    *ep = std::move(g.get_property(e));
  }
  g.m_pack.remove_edge(e);
}

/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

// Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
template <BAGL_ADJACENCY_LIST_ARGS>
const auto& get(const BAGL_ADJACENCY_LIST& g, typename BAGL_ADJACENCY_LIST::vertex_descriptor v) {
  return g[v];
}
template <BAGL_ADJACENCY_LIST_ARGS>
auto& get(BAGL_ADJACENCY_LIST& g, typename BAGL_ADJACENCY_LIST::vertex_descriptor v) {
  return g[v];
}

// Returns a const-reference to the edge-bundle associated to the given edge descriptor.
template <BAGL_ADJACENCY_LIST_ARGS>
const auto& get(const BAGL_ADJACENCY_LIST& g, const typename BAGL_ADJACENCY_LIST::edge_descriptor& e) {
  return g[e];
}
template <BAGL_ADJACENCY_LIST_ARGS>
auto& get(BAGL_ADJACENCY_LIST& g, const typename BAGL_ADJACENCY_LIST::edge_descriptor& e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <BAGL_ADJACENCY_LIST_ARGS>
const auto& get(const BAGL_ADJACENCY_LIST& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}
template <BAGL_ADJACENCY_LIST_ARGS>
auto& get(BAGL_ADJACENCY_LIST& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <BAGL_ADJACENCY_LIST_ARGS, typename VProp>
void put(BAGL_ADJACENCY_LIST& g, typename BAGL_ADJACENCY_LIST::vertex_descriptor v, VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-bundle associated to the given edge descriptor.
template <BAGL_ADJACENCY_LIST_ARGS, typename EProp>
void put(BAGL_ADJACENCY_LIST& g, const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Sets the graph-bundle associated to the graph.
template <BAGL_ADJACENCY_LIST_ARGS, typename GProp>
void put(BAGL_ADJACENCY_LIST& g, graph_bundle_t /*unused*/, GProp&& value) {
  g[graph_bundle] = std::forward<GProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_ADJACENCY_LIST_ARGS>
auto& get_property(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, BAGL_ADJACENCY_LIST& g) {
  return g.get_property(v);
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_ADJACENCY_LIST_ARGS>
const auto& get_property(typename BAGL_ADJACENCY_LIST::vertex_descriptor v, const BAGL_ADJACENCY_LIST& g) {
  return g.get_property(v);
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <BAGL_ADJACENCY_LIST_ARGS>
auto& get_property(const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, BAGL_ADJACENCY_LIST& g) {
  return g.get_property(e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <BAGL_ADJACENCY_LIST_ARGS>
const auto& get_property(const typename BAGL_ADJACENCY_LIST::edge_descriptor& e, const BAGL_ADJACENCY_LIST& g) {
  return g.get_property(e);
}

template <BAGL_ADJACENCY_LIST_ARGS>
auto& get_property(BAGL_ADJACENCY_LIST& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

template <BAGL_ADJACENCY_LIST_ARGS>
const auto& get_property(const BAGL_ADJACENCY_LIST& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Handle graph property tags, also handles graph_bundle_t.

template <BAGL_ADJACENCY_LIST_ARGS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>, lookup_one_property_t<GraphProperties, Tag>&>
get_property(BAGL_ADJACENCY_LIST& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_ADJACENCY_LIST_ARGS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>,
                 const lookup_one_property_t<GraphProperties, Tag>&>
get_property(const BAGL_ADJACENCY_LIST& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle>
struct property_map<BAGL_ADJACENCY_LIST, T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v =
      std::is_convertible_v<typename BAGL_ADJACENCY_LIST::vertex_bundled*, non_const_Bundle*>;
  static constexpr bool is_edge_bundle_v =
      std::is_convertible_v<typename BAGL_ADJACENCY_LIST::edge_bundled*, non_const_Bundle*>;
  using tag_type = std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                                      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<non_const_T, BAGL_ADJACENCY_LIST, tag_type>;
  using const_type = bundle_member_property_map<const non_const_T, const BAGL_ADJACENCY_LIST, tag_type>;
};

template <BAGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, BAGL_ADJACENCY_LIST& g) {
  return typename property_map<BAGL_ADJACENCY_LIST, T Bundle::*>::type{&g, p};
}

template <BAGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, const BAGL_ADJACENCY_LIST& g) {
  return typename property_map<BAGL_ADJACENCY_LIST, T Bundle::*>::const_type{&g, p};
}

template <BAGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p, const BAGL_ADJACENCY_LIST& g, const Key& k) {
  return (g[k]).*p;
}

template <BAGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_ADJACENCY_LIST& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct adj_list_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;

    using type = tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type = tagged_in_property_property_map<const value_type, const Graph, Tag>;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<adj_list_tag> {
  using type = adj_list_property_selector;
};

template <>
struct edge_property_selector<adj_list_tag> {
  using type = adj_list_property_selector;
};

template <BAGL_ADJACENCY_LIST_ARGS, typename Property>
auto get(Property p, BAGL_ADJACENCY_LIST& g) {
  using Map = property_map_t<BAGL_ADJACENCY_LIST, Property>;
  return Map(&g, p);
}

template <BAGL_ADJACENCY_LIST_ARGS, typename Property>
auto get(Property p, const BAGL_ADJACENCY_LIST& g) {
  using Map = property_map_const_t<BAGL_ADJACENCY_LIST, Property>;
  return Map(&g, p);
}

template <BAGL_ADJACENCY_LIST_ARGS, typename Property, typename Key>
decltype(auto) get(Property p, const BAGL_ADJACENCY_LIST& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_ADJACENCY_LIST_ARGS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_ADJACENCY_LIST& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_ADJACENCY_LIST_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BAGL_ADJACENCY_LIST& g, const Key& k, Value&& val) {
  get_property_value(g.get_property(k), p) = std::forward<Value>(val);
}

template <BAGL_ADJACENCY_LIST_ARGS>
void do_graph_deep_copy(BAGL_ADJACENCY_LIST& lhs, const BAGL_ADJACENCY_LIST& rhs) {
  using Vertex = typename BAGL_ADJACENCY_LIST::vertex_descriptor;

  lhs.m_pack.clear();

  // first, add all the vertices (keep unordered_map of correspondance).
  std::unordered_map<Vertex, Vertex, container_detail::desc_hasher> v_map;
  for (auto u : vertices(rhs)) {
    Vertex v = add_vertex(lhs, rhs.get_property(u));
    v_map[u] = v;
  }

  // then, go through all the edges and add them to the lhs:
  for (auto e : edges(rhs)) {
    add_edge(v_map[source(e, rhs)], v_map[target(e, rhs)], lhs, rhs.get_property(e));
  }
}

#undef BAGL_ADJACENCY_LIST_ARGS
#undef BAGL_ADJACENCY_LIST

/***********************************************************************************************
 *                             Vertex Index Map for vec_s and pool_s
 * ********************************************************************************************/

// Note that the vertex index values for a pool_s vertex list will not be contiguous or ordered
// since holes in the vertex pool retain their associated index (index into raw vector).

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
struct property_map<adjacency_list<OutEdgeListS, vec_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>,
                    vertex_index_t> {
  using type = typed_identity_property_map<std::size_t>;
  using const_type = typed_identity_property_map<std::size_t>;
};

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
auto get(vertex_index_t /*unused*/,
         adjacency_list<OutEdgeListS, vec_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/) {
  return typed_identity_property_map<std::size_t>{};
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
auto get(
    vertex_index_t /*unused*/,
    const adjacency_list<OutEdgeListS, vec_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/) {
  return typed_identity_property_map<std::size_t>{};
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
std::size_t get(
    vertex_index_t /*unused*/,
    const adjacency_list<OutEdgeListS, vec_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/,
    const std::size_t& k) {
  return k;
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
std::size_t get(
    vertex_index_t /*unused*/,
    adjacency_list<OutEdgeListS, vec_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/,
    const std::size_t& k) {
  return k;
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
struct property_map<adjacency_list<OutEdgeListS, pool_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>,
                    vertex_index_t> {
  using type = typed_identity_property_map<std::size_t>;
  using const_type = typed_identity_property_map<std::size_t>;
};

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
auto get(vertex_index_t /*unused*/,
         adjacency_list<OutEdgeListS, pool_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/) {
  return typed_identity_property_map<std::size_t>{};
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
auto get(
    vertex_index_t /*unused*/,
    const adjacency_list<OutEdgeListS, pool_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/) {
  return typed_identity_property_map<std::size_t>{};
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
std::size_t get(
    vertex_index_t /*unused*/,
    adjacency_list<OutEdgeListS, pool_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/,
    const std::size_t& k) {
  return k;
}

template <typename OutEdgeListS, typename DirectedS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
std::size_t get(
    vertex_index_t /*unused*/,
    const adjacency_list<OutEdgeListS, pool_s, DirectedS, VertexProperties, EdgeProperties, GraphProperties>& /*g*/,
    const std::size_t& k) {
  return k;
}

/**
 * This class implements an adjacency-list based on Boost.Containers that is tailored
 * to store elements of an undirected graph.
 * \note This is a partial specialization for undirected edges.
 * \tparam OutEdgeListS A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexListS A type tag to choose the storage policy for the vertices.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS, typename VertexListS, typename VertexProperties, typename EdgeProperties,
          typename GraphProperties>
class adjacency_list<OutEdgeListS, VertexListS, undirected_s, VertexProperties, EdgeProperties, GraphProperties> {
 public:
  using check_allowed_vertex_list = typename adjacency_list_disallowed_vertex_list<VertexListS>::type;

  using self =
      adjacency_list<OutEdgeListS, VertexListS, undirected_s, VertexProperties, EdgeProperties, GraphProperties>;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;
  using graph_property_type = GraphProperties;

  using vertex_bundled = lookup_one_property_t<VertexProperties, vertex_bundle_t>;
  using edge_bundled = lookup_one_property_t<EdgeProperties, edge_bundle_t>;
  using graph_bundled = lookup_one_property_t<GraphProperties, graph_bundle_t>;

  using bidir_storage_type =
      typename adjlist_detail::adjlist_vertex_container<VertexListS, OutEdgeListS, bidirectional_s, VertexProperties,
                                                        EdgeProperties>;

  using vertex_descriptor = typename bidir_storage_type::vertex_descriptor;
  using vertices_size_type = typename bidir_storage_type::vertices_size_type;

  using bidir_edge_descriptor = typename bidir_storage_type::edge_descriptor;
  using edge_descriptor = container_detail::undir_edge_desc<bidir_edge_descriptor>;
  using edges_size_type = typename bidir_storage_type::edges_size_type;
  using degree_size_type = edges_size_type;

  using vertex_stored_impl = typename bidir_storage_type::vertex_stored_type;
  using vertex_value_impl = typename bidir_storage_type::vertex_value_type;

  using edge_stored_impl = typename bidir_storage_type::edge_stored_type;
  using edge_value_impl = typename bidir_storage_type::edge_value_type;

  using Traits = adjacency_list_traits<OutEdgeListS, VertexListS, undirected_s>;

  using directed_category = typename Traits::directed_category;
  using edge_parallel_category = typename Traits::edge_parallel_category;
  using traversal_category = typename Traits::traversal_category;

  using mutability_category = mutable_property_graph_tag;

  using graph_tag = adj_list_tag;

  /**
   * This static member function outputs the null-vertex (invalid vertex descriptor).
   * \return A null-vertex descriptor (invalid vertex descriptor).
   */
  static vertex_descriptor null_vertex() { return container_detail::null_desc<vertex_descriptor>::value(); }

  /**
   * This static member function outputs the null-edge (invalid edge descriptor).
   * \return A null-edge descriptor (invalid edge descriptor).
   */
  static edge_descriptor null_edge() {
    return edge_descriptor(container_detail::null_desc<bidir_edge_descriptor>::value());
  }

  // private:
  bidir_storage_type m_pack;
  graph_property_type m_graph_prop;

  /**
   * Constructs an empty adjacency-list.
   */
  adjacency_list() : m_pack(), m_graph_prop() {}

  explicit adjacency_list(std::size_t num_vertices, graph_property_type graph_prop = {})
      : m_pack(num_vertices), m_graph_prop(std::move(graph_prop)) {}

  ~adjacency_list() = default;

  adjacency_list(const self& rhs) : m_pack(), m_graph_prop(rhs.m_graph_prop) { do_graph_deep_copy(*this, rhs); }
  self& operator=(const self& rhs) {
    if (this != &rhs) {
      m_graph_prop = rhs.m_graph_prop;
      do_graph_deep_copy(*this, rhs);
    }
    return *this;
  }

  adjacency_list(self&& rhs) noexcept : m_pack(std::move(rhs.m_pack)), m_graph_prop(rhs.m_graph_prop) {}
  self& operator=(self&& rhs) noexcept {
    m_pack = std::move(rhs.m_pack);
    m_graph_prop = std::move(rhs.m_graph_prop);
    return *this;
  }

  // Construct from a given number of vertices and an edge range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
  adjacency_list(vertices_size_type num_vertices, const EdgeRange& e_range, graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto& v : tmp_vs) {
      v = m_pack.add_vertex(vertex_property_type{});
    }

    for (auto [u, v] : e_range) {
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], edge_property_type{});
    }
  }

  // Construct from a given number of vertices and an edge and edge-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
      std::constructible_from<edge_property_type, std::ranges::range_reference_t<EdgePropRange>>
      adjacency_list(vertices_size_type num_vertices, const EdgeRange& e_range, const EdgePropRange& ep_range,
                     graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto& v : tmp_vs) {
      v = m_pack.add_vertex(vertex_property_type{});
    }

    for (auto [e, ep] : zip_range_ref(e_range, ep_range)) {
      auto [u, v] = e;
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], ep);
    }
  }

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
      std::constructible_from<vertex_property_type, std::ranges::range_reference_t<VertexPropRange>>
      adjacency_list(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                     graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto [v, vp] : zip_range_ref(tmp_vs, vp_range)) {
      v = m_pack.add_vertex(vp);
    }

    for (auto [u, v] : e_range) {
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], edge_property_type{});
    }
  }

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange,
            std::ranges::input_range EdgePropRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
      std::constructible_from<vertex_property_type, std::ranges::range_reference_t<VertexPropRange>> &&
      std::constructible_from<edge_property_type, std::ranges::range_reference_t<EdgePropRange>>
      adjacency_list(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                     const EdgePropRange& ep_range, graph_property_type graph_prop = {})
      : m_pack(), m_graph_prop(std::move(graph_prop)) {
    std::vector<vertex_descriptor> tmp_vs(num_vertices);
    for (auto [v, vp] : zip_range_ref(tmp_vs, vp_range)) {
      v = m_pack.add_vertex(vp);
    }

    for (auto [e, ep] : zip_range_ref(e_range, ep_range)) {
      auto [u, v] = e;
      m_pack.add_edge(tmp_vs[u], tmp_vs[v], ep);
    }
  }

  /**
   * Swaps the adjacency-list with another.
   */
  void swap(self& rhs) {
    using std::swap;
    m_pack.swap(rhs.m_pack);
    swap(m_graph_prop, rhs.m_graph_prop);
  }

  /**
   * Clears the adjacency-list of all vertices and edges.
   */
  void clear() { m_pack.clear(); }

  // Indexing operator. Returns a reference to the vertex-bundle associated to the given vertex descriptor.
  vertex_bundled& operator[](vertex_descriptor v) {
    return get_property_value(m_pack.get_stored_vertex(v).data, vertex_bundle);
  }
  // Indexing operator. Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
  const vertex_bundled& operator[](vertex_descriptor v) const {
    return get_property_value(m_pack.get_stored_vertex(v).data, vertex_bundle);
  }

  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  edge_bundled& operator[](const edge_descriptor& e) {
    return get_property_value(m_pack.get_stored_edge_property(bidir_edge_descriptor(e)), edge_bundle);
  }
  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  const edge_bundled& operator[](const edge_descriptor& e) const {
    return get_property_value(m_pack.get_stored_edge_property(bidir_edge_descriptor(e)), edge_bundle);
  }

  // Indexing operator. Returns a reference to the graph-bundle associated to the graph.
  edge_bundled& operator[](graph_bundle_t /*unused*/) { return get_property_value(m_graph_prop, graph_bundle); }
  // Indexing operator. Returns a const-reference to the graph-bundle associated to the graph.
  const edge_bundled& operator[](graph_bundle_t /*unused*/) const {
    return get_property_value(m_graph_prop, graph_bundle);
  }

  // Get a reference to the vertex-property associated to the given vertex descriptor.
  auto& get_property(vertex_descriptor v) { return m_pack.get_stored_vertex(v).data; }
  // Get a const-reference to the vertex-property associated to the given vertex descriptor.
  const auto& get_property(vertex_descriptor v) const { return m_pack.get_stored_vertex(v).data; }

  // Get a reference to the edge-property associated to the given edge descriptor.
  auto& get_property(const edge_descriptor& e) { return m_pack.get_stored_edge_property(bidir_edge_descriptor(e)); }
  // Get a const-reference to the edge-property associated to the given edge descriptor.
  const auto& get_property(const edge_descriptor& e) const {
    return m_pack.get_stored_edge_property(bidir_edge_descriptor(e));
  }

  // Get a reference to the graph-property associated to the graph.
  auto& get_property(graph_all_t /*unused*/) { return m_graph_prop; }
  // Get a const-reference to the graph-property associated to the graph.
  const auto& get_property(graph_all_t /*unused*/) const { return m_graph_prop; }
};

#define BAGL_ADJACENCY_LIST_UNDIR_ARGS                                                             \
  typename OutEdgeListS, typename VertexListS, typename VertexProperties, typename EdgeProperties, \
      typename GraphProperties
#define BAGL_ADJACENCY_LIST_UNDIR \
  adjacency_list<OutEdgeListS, VertexListS, undirected_s, VertexProperties, EdgeProperties, GraphProperties>

/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto source(const typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor& e, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  if (e.is_reversed) {
    return g.m_pack.get_stored_edge_target(e);
  }
  return e.source;
}

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto target(const typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor& e, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  if (e.is_reversed) {
    return e.source;
  }
  return g.m_pack.get_stored_edge_target(e);
}

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto out_edges(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return make_adjlist_undir_ioerange<typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor>(true, g.m_pack.in_edges(v),
                                                                                          g.m_pack.out_edges(v));
}

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
std::size_t out_degree(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.get_out_degree(v) + g.m_pack.get_in_degree(v);
}

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto in_edges(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return make_adjlist_undir_ioerange<typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor>(false, g.m_pack.in_edges(v),
                                                                                          g.m_pack.out_edges(v));
}

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
std::size_t in_degree(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
}

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
std::size_t degree(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
}

/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto edges(const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return adjlist_detail::adjlist_undir_eiter_range(g.m_pack.edges());
}

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto num_edges(const BAGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.num_edges();
}

/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto edge(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
          typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BAGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  auto [be, be_found] = g.m_pack.get_edge(u, v);
  if (be_found) {
    return std::pair(Edge(be), true);
  }
  std::tie(be, be_found) = g.m_pack.get_edge(v, u);
  if (be_found) {
    return std::pair(Edge(be, true), true);
  }
  return std::pair(Edge(be), false);
}

/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS>
auto add_edge(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
              typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, BAGL_ADJACENCY_LIST_UNDIR& g) {
  using EdgeProp = typename BAGL_ADJACENCY_LIST_UNDIR::edge_property_type;
  using Edge = typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  auto [be, be_found] = g.m_pack.add_edge(u, v, EdgeProp{});
  return std::pair(Edge(be), be_found);
}

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS, typename EProp>
auto add_edge(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
              typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, EProp&& ep, BAGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  auto [be, be_found] = g.m_pack.add_edge(u, v, std::forward<EProp>(ep));
  return std::pair(Edge(be), be_found);
}

/*******************************************************************************************
 *                    MutableIncidenceGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS, typename EdgePred>
void remove_out_edge_if(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u, EdgePred pred,
                        BAGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  g.m_pack.remove_out_edge_if(u, [&pred](const auto& be) { return pred(Edge(be, /*reversed=*/false)); });
  g.m_pack.remove_in_edge_if(u, [&pred](const auto& be) { return pred(Edge(be, /*reversed=*/true)); });
}

/*******************************************************************************************
 *                    MutableBidirectionalGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS, typename EdgePred>
void remove_in_edge_if(typename BAGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, EdgePred pred,
                       BAGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  g.m_pack.remove_out_edge_if(v, [&pred](const auto& be) { return pred(Edge(be, /*reversed=*/true)); });
  g.m_pack.remove_in_edge_if(v, [&pred](const auto& be) { return pred(Edge(be, /*reversed=*/false)); });
}

/*******************************************************************************************
 *                    MutableEdgeListGraph concept
 ******************************************************************************************/

template <BAGL_ADJACENCY_LIST_UNDIR_ARGS, typename EdgePred>
void remove_edge_if(EdgePred pred, BAGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BAGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  g.m_pack.remove_edge_if([&pred](const auto& be) { return pred(Edge(be)); });
}

#undef BAGL_ADJACENCY_LIST_UNDIR_ARGS
#undef BAGL_ADJACENCY_LIST_UNDIR

}  // namespace bagl

#endif  // BAGL_BAGL_ADJACENCY_LIST_H_
