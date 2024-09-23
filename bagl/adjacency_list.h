// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_ADJACENCY_LIST_H_
#define BAGL_BAGL_ADJACENCY_LIST_H_

#include <unordered_map>
#include <utility>

#include "bagl/adjacency_range.h"
#include "bagl/detail/adjlist_containers.h"
#include "bagl/detail/adjlist_ranges.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_mutability_traits.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/tree_traits.h"

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
          typename VertexProperties = no_property, typename EdgeProperties = no_property>
class adjacency_list;

template <typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
void do_graph_deep_copy(
    adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>& lhs,
    const adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>& rhs);

/**
 * This class implements an adjacency-list based on Boost.Containers that is tailored
 * to store elements of a graph.
 * \tparam OutEdgeListS A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexListS A type tag to choose the storage policy for the vertices.
 * \tparam DirectedS A type tag to choose the directionality of the edges.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties,
          typename EdgeProperties>
class adjacency_list {
 public:
  using check_allowed_vertex_list = typename adjacency_list_disallowed_vertex_list<VertexListS>::type;

  using self = adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>;

  using out_edge_list_selector = OutEdgeListS;
  using vertex_list_selector = VertexListS;
  using directed_selector = DirectedS;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;

  using vertex_bundled = VertexProperties;
  using edge_bundled = EdgeProperties;

  using storage_type = typename adjlist_detail::adjlist_vertex_container<VertexListS, OutEdgeListS, DirectedS,
                                                                         VertexProperties, EdgeProperties>;

  using vertex_descriptor = typename storage_type::vertex_descriptor;
  using vertices_size_type = typename storage_type::vertices_size_type;

  using edge_descriptor = typename storage_type::edge_descriptor;
  using edges_size_type = typename storage_type::edges_size_type;
  using degree_size_type = edges_size_type;

  using vertex_range = typename storage_type::vertex_range;

  using out_edge_range = typename storage_type::out_edge_range;
  using in_edge_range = typename storage_type::in_edge_range;
  using edge_range = typename storage_type::edge_range;

  using adjacency_range = decltype(adjacency_range(std::declval<out_edge_range>(), std::declval<self>()));
  using inv_adjacency_range = decltype(inv_adjacency_range(std::declval<in_edge_range>(), std::declval<self>()));

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

  /**
   * Constructs an empty adjacency-list.
   */
  adjacency_list() : m_pack() {}

  ~adjacency_list() = default;

  adjacency_list(const self& rhs) : m_pack() { do_graph_deep_copy(*this, rhs); }
  self& operator=(const self& rhs) {
    if (this != &rhs) {
      do_graph_deep_copy(*this, rhs);
    }
    return *this;
  }

  adjacency_list(self&& rhs) noexcept : m_pack(std::move(rhs.m_pack)) {}
  self& operator=(self&& rhs) noexcept {
    m_pack = std::move(rhs.m_pack);
    return *this;
  }

  /**
   * Swaps the adjacency-list with another.
   */
  void swap(self& rhs) { m_pack.swap(rhs.m_pack); }

  /**
   * Clears the adjacency-list of all vertices and edges.
   */
  void clear() { m_pack.clear(); }

  /**
   * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
   * \param v The vertex descriptor of the sought-after vertex-property.
   * \return The vertex-property, by reference, associated to the given vertex descriptor.
   */
  vertex_property_type& operator[](vertex_descriptor v) { return m_pack.get_stored_vertex(v).data; }
  /**
   * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
   * \param v The vertex descriptor of the sought-after vertex-property.
   * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
   */
  const vertex_property_type& operator[](vertex_descriptor v) const { return m_pack.get_stored_vertex(v).data; }
  /**
   * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
   * \param e The edge descriptor of the sought-after edge-property.
   * \return The edge-property, by reference, associated to the given edge descriptor.
   */
  edge_property_type& operator[](const edge_descriptor& e) { return m_pack.get_stored_edge(e).data; }
  /**
   * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
   * \param e The edge descriptor of the sought-after edge-property.
   * \return The edge-property, by const-reference, associated to the given edge descriptor.
   */
  const edge_property_type& operator[](const edge_descriptor& e) const { return m_pack.get_stored_edge(e).data; }
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
  using child_vertex_range =
      typename adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperty, EdgeProperty>::adjacency_range;
};

#define BGL_ADJACENCY_LIST_ARGS \
  typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties, typename EdgeProperties
#define BGL_ADJACENCY_LIST adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>

template <BGL_ADJACENCY_LIST_ARGS>
void swap(BGL_ADJACENCY_LIST& lhs, BGL_ADJACENCY_LIST& rhs) {
  lhs.swap(rhs);
}

/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto source(const typename BGL_ADJACENCY_LIST::edge_descriptor& e, const BGL_ADJACENCY_LIST& /*unused*/) {
  return e.source;
}

template <BGL_ADJACENCY_LIST_ARGS>
auto target(const typename BGL_ADJACENCY_LIST::edge_descriptor& e, const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_stored_edge(e).target;
}

template <BGL_ADJACENCY_LIST_ARGS>
auto out_edges(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.out_edges(v);
}

template <BGL_ADJACENCY_LIST_ARGS>
std::size_t out_degree(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_out_degree(v);
}

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto in_edges(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.in_edges(v);
}

template <BGL_ADJACENCY_LIST_ARGS>
std::size_t in_degree(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_in_degree(v);
}

template <BGL_ADJACENCY_LIST_ARGS>
std::size_t degree(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
}

/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto vertices(const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.vertices();
}

template <BGL_ADJACENCY_LIST_ARGS>
auto num_vertices(const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.size();
}

/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto edges(const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.edges();
}

template <BGL_ADJACENCY_LIST_ARGS>
auto num_edges(const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.num_edges();
}

/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto adjacent_vertices(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return adjacency_range(out_edges(v, g), g);
}

/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto inv_adjacent_vertices(typename BGL_ADJACENCY_LIST::vertex_descriptor v, const BGL_ADJACENCY_LIST& g) {
  return inv_adjacency_range(in_edges(v, g), g);
}

/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto edge(typename BGL_ADJACENCY_LIST::vertex_descriptor u, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
          const BGL_ADJACENCY_LIST& g) {
  return g.m_pack.get_edge(u, v);
}

/*******************************************************************************************
 *                  MutableGraph concept
 ******************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto add_vertex(BGL_ADJACENCY_LIST& g) {
  using VertexBundled = typename BGL_ADJACENCY_LIST::vertex_bundled;
  return g.m_pack.add_vertex(VertexBundled());
}

template <BGL_ADJACENCY_LIST_ARGS>
void clear_vertex(typename BGL_ADJACENCY_LIST::vertex_descriptor v, BGL_ADJACENCY_LIST& g) {
  g.m_pack.clear_vertex(v);
}

template <BGL_ADJACENCY_LIST_ARGS>
void remove_vertex(typename BGL_ADJACENCY_LIST::vertex_descriptor v, BGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_vertex(v);
}

template <BGL_ADJACENCY_LIST_ARGS>
auto add_edge(typename BGL_ADJACENCY_LIST::vertex_descriptor u, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
              BGL_ADJACENCY_LIST& g) {
  using EdgeBundled = typename BGL_ADJACENCY_LIST::edge_bundled;
  return g.m_pack.add_edge(u, v, EdgeBundled());
}

template <BGL_ADJACENCY_LIST_ARGS>
void remove_edge(typename BGL_ADJACENCY_LIST::vertex_descriptor u, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
                 BGL_ADJACENCY_LIST& g) {
  auto [e, e_found] = edge(u, v, g);
  if (e_found) {
    g.m_pack.remove_edge(e);
  }
}

template <BGL_ADJACENCY_LIST_ARGS>
void remove_edge(const typename BGL_ADJACENCY_LIST::edge_descriptor& e, BGL_ADJACENCY_LIST& g) {
  g.m_pack.remove_edge(e);
}

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <BGL_ADJACENCY_LIST_ARGS>
auto add_vertex(const typename BGL_ADJACENCY_LIST::vertex_property_type& vp, BGL_ADJACENCY_LIST& g) {
  return g.m_pack.add_vertex(vp);
}

template <BGL_ADJACENCY_LIST_ARGS>
void remove_vertex(typename BGL_ADJACENCY_LIST::vertex_descriptor v,
                   typename BGL_ADJACENCY_LIST::vertex_property_type& vp, BGL_ADJACENCY_LIST& g) {
  vp = std::move(g[v]);
  g.m_pack.remove_vertex(v);
}

template <BGL_ADJACENCY_LIST_ARGS>
auto add_edge(typename BGL_ADJACENCY_LIST::vertex_descriptor u, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
              const typename BGL_ADJACENCY_LIST::edge_property_type& ep, BGL_ADJACENCY_LIST& g) {
  return g.m_pack.add_edge(u, v, ep);
}

template <BGL_ADJACENCY_LIST_ARGS>
void remove_edge(typename BGL_ADJACENCY_LIST::vertex_descriptor u, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
                 typename BGL_ADJACENCY_LIST::edge_property_type& ep, BGL_ADJACENCY_LIST& g) {
  auto [e, e_found] = edge(u, v, g);
  if (e_found) {
    ep = std::move(g[e]);
    g.m_pack.remove_edge(e);
  }
}

template <BGL_ADJACENCY_LIST_ARGS>
void remove_edge(const typename BGL_ADJACENCY_LIST::edge_descriptor& e,
                 typename BGL_ADJACENCY_LIST::edge_property_type& ep, BGL_ADJACENCY_LIST& g) {
  ep = std::move(g[e]);
  g.m_pack.remove_edge(e);
}

template <BGL_ADJACENCY_LIST_ARGS>
auto add_vertex(typename BGL_ADJACENCY_LIST::vertex_property_type&& vp, BGL_ADJACENCY_LIST& g) {
  return g.m_pack.add_vertex(std::move(vp));
}

template <BGL_ADJACENCY_LIST_ARGS>
auto add_edge(typename BGL_ADJACENCY_LIST::vertex_descriptor u, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
              typename BGL_ADJACENCY_LIST::edge_property_type&& ep, BGL_ADJACENCY_LIST& g) {
  return g.m_pack.add_edge(u, v, std::move(ep));
}

/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

/**
 * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
 * \param g The tree from which to draw the vertex.
 * \param v The vertex descriptor of the sought-after vertex-property.
 * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
 */
template <BGL_ADJACENCY_LIST_ARGS>
const typename BGL_ADJACENCY_LIST::vertex_property_type& get(const BGL_ADJACENCY_LIST& g,
                                                             typename BGL_ADJACENCY_LIST::vertex_descriptor v) {
  return g[v];
}

/**
 * Returns a const-reference to the edge-property associated to the given edge descriptor.
 * \param g The tree from which to draw the edge.
 * \param e The edge descriptor of the sought-after edge-property.
 * \return The edge-property, by const-reference, associated to the given edge descriptor.
 */
template <BGL_ADJACENCY_LIST_ARGS>
const typename BGL_ADJACENCY_LIST::edge_property_type& get(const BGL_ADJACENCY_LIST& g,
                                                           const typename BGL_ADJACENCY_LIST::edge_descriptor& e) {
  return g[e];
}

/**
 * Sets the vertex-property associated to the given vertex descriptor.
 * \param g The tree from which the vertex is drawn.
 * \param v The vertex descriptor of the vertex-property to be set.
 * \param value The vertex-property, by const-reference, to be associated to the given vertex.
 */
template <BGL_ADJACENCY_LIST_ARGS>
void put(BGL_ADJACENCY_LIST& g, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
         const typename BGL_ADJACENCY_LIST::vertex_property_type& value) {
  g[v] = value;
}

/**
 * Sets the edge-property associated to the given edge descriptor.
 * \param g The tree from which the edge is drawn.
 * \param e The edge descriptor of the edge-property to be set.
 * \param value The edge-property, by const-reference, to be associated to the given edge.
 */
template <BGL_ADJACENCY_LIST_ARGS>
void put(BGL_ADJACENCY_LIST& g, const typename BGL_ADJACENCY_LIST::edge_descriptor& e,
         const typename BGL_ADJACENCY_LIST::edge_property_type& value) {
  g[e] = value;
}

/**
 * Sets the vertex-property associated to the given vertex descriptor.
 * \param g The tree from which the vertex is drawn.
 * \param v The vertex descriptor of the vertex-property to be set.
 * \param value The vertex-property, by rvalue-reference, to be associated to the given vertex.
 */
template <BGL_ADJACENCY_LIST_ARGS>
void put(BGL_ADJACENCY_LIST& g, typename BGL_ADJACENCY_LIST::vertex_descriptor v,
         typename BGL_ADJACENCY_LIST::vertex_property_type&& value) {
  g[v] = std::move(value);
}

/**
 * Sets the edge-property associated to the given edge descriptor.
 * \param g The tree from which the edge is drawn.
 * \param e The edge descriptor of the edge-property to be set.
 * \param value The edge-property, by rvalue-reference, to be associated to the given edge.
 */
template <BGL_ADJACENCY_LIST_ARGS>
void put(BGL_ADJACENCY_LIST& g, const typename BGL_ADJACENCY_LIST::edge_descriptor& e,
         typename BGL_ADJACENCY_LIST::edge_property_type&& value) {
  g[e] = std::move(value);
}

/**
 * Returns a reference to the vertex-property associated to the given vertex descriptor.
 * \param v The vertex descriptor of the sought-after vertex-property.
 * \param g The tree from which to draw the vertex.
 * \return The vertex-property, by reference, associated to the given vertex descriptor.
 */
template <BGL_ADJACENCY_LIST_ARGS>
typename BGL_ADJACENCY_LIST::vertex_property_type& get_property(typename BGL_ADJACENCY_LIST::vertex_descriptor v,
                                                                BGL_ADJACENCY_LIST& g) {
  return g[v];
}

/**
 * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
 * \param v The vertex descriptor of the sought-after vertex-property.
 * \param g The tree from which to draw the vertex.
 * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
 */
template <BGL_ADJACENCY_LIST_ARGS>
const typename BGL_ADJACENCY_LIST::vertex_property_type& get_property(typename BGL_ADJACENCY_LIST::vertex_descriptor v,
                                                                      const BGL_ADJACENCY_LIST& g) {
  return g[v];
}

/**
 * Returns a reference to the edge-property associated to the given edge descriptor.
 * \param e The edge descriptor of the sought-after edge-property.
 * \param g The tree from which to draw the edge.
 * \return The edge-property, by reference, associated to the given edge descriptor.
 */
template <BGL_ADJACENCY_LIST_ARGS>
typename BGL_ADJACENCY_LIST::edge_property_type& get_property(const typename BGL_ADJACENCY_LIST::edge_descriptor& e,
                                                              BGL_ADJACENCY_LIST& g) {
  return g[e];
}

/**
 * Returns a const-reference to the edge-property associated to the given edge descriptor.
 * \param e The edge descriptor of the sought-after edge-property.
 * \param g The tree from which to draw the edge.
 * \return The edge-property, by const-reference, associated to the given edge descriptor.
 */
template <BGL_ADJACENCY_LIST_ARGS>
const typename BGL_ADJACENCY_LIST::edge_property_type& get_property(
    const typename BGL_ADJACENCY_LIST::edge_descriptor& e, const BGL_ADJACENCY_LIST& g) {
  return g[e];
}

template <BGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle>
struct property_map<BGL_ADJACENCY_LIST, T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v = std::is_convertible_v<
      typename adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>::vertex_bundled*,
      non_const_Bundle*>;
  using type =
      bundle_member_property_map<non_const_T,
                                 adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>,
                                 std::conditional_t<is_vertex_bundle_v, vertex_bundle_t, edge_bundle_t>>;
  using const_type = bundle_member_property_map<
      const non_const_T, const adjacency_list<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>,
      std::conditional_t<is_vertex_bundle, vertex_bundle_t, edge_bundle_t>>;
};

template <BGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle>
typename property_map<BGL_ADJACENCY_LIST, T Bundle::*>::type get(T Bundle::*p, BGL_ADJACENCY_LIST& g) {
  return {&g, p};
}

template <BGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle>
typename property_map<BGL_ADJACENCY_LIST, T Bundle::*>::const_type get(T Bundle::*p, const BGL_ADJACENCY_LIST& g) {
  return {&g, p};
}

template <BGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p, const BGL_ADJACENCY_LIST& g, const Key& k) {
  return (g[k]).*p;
}

template <BGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BGL_ADJACENCY_LIST& g, const Key& k, const T& val) {
  (g[k]).*p = val;
}

template <BGL_ADJACENCY_LIST_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BGL_ADJACENCY_LIST& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct adj_list_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;

    using type = tagged_from_bundle_property_map<value_type, Graph, Tag>;
    using const_type = tagged_from_bundle_property_map<const value_type, const Graph, Tag>;
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

template <BGL_ADJACENCY_LIST_ARGS, typename Property>
auto get(Property p, BGL_ADJACENCY_LIST& g) {
  using Map = typename property_map<BGL_ADJACENCY_LIST, Property>::type;
  return Map(&g, p);
}

template <BGL_ADJACENCY_LIST_ARGS, typename Property>
auto get(Property p, const BGL_ADJACENCY_LIST& g) {
  using Map = typename property_map<BGL_ADJACENCY_LIST, Property>::const_type;
  return Map(&g, p);
}

template <BGL_ADJACENCY_LIST_ARGS, typename Property, typename Key>
auto get(Property p, const BGL_ADJACENCY_LIST& g, const Key& k) {
  return get_property_value(g[k], p);
}

template <BGL_ADJACENCY_LIST_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_ADJACENCY_LIST& g, const Key& k, const Value& val) {
  get_property_value(g[k], p) = val;
}

template <BGL_ADJACENCY_LIST_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BGL_ADJACENCY_LIST& g, const Key& k, Value&& val) {
  get_property_value(g[k], p) = std::forward<Value>(val);
}

template <BGL_ADJACENCY_LIST_ARGS>
void do_graph_deep_copy(BGL_ADJACENCY_LIST& lhs, const BGL_ADJACENCY_LIST& rhs) {
  using Vertex = typename BGL_ADJACENCY_LIST::vertex_descriptor;

  lhs.m_pack.clear();

  // first, add all the vertices (keep unordered_map of correspondance).
  std::unordered_map<Vertex, Vertex, container_detail::desc_hasher> v_map;
  for (auto u : vertices(rhs)) {
    Vertex v = add_vertex(rhs[u], lhs);
    v_map[u] = v;
  }

  // then, go through all the edges and add them to the lhs:
  for (auto e : edges(rhs)) {
    add_edge(v_map[source(e, rhs)], v_map[target(e, rhs)], rhs[e], lhs);
  }
}

#undef BGL_ADJACENCY_LIST_ARGS
#undef BGL_ADJACENCY_LIST

/**
 * This class implements an adjacency-list based on Boost.Containers that is tailored
 * to store elements of an undirected graph.
 * \note This is a partial specialization for undirected edges.
 * \tparam OutEdgeListS A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexListS A type tag to choose the storage policy for the vertices.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS, typename VertexListS, typename VertexProperties, typename EdgeProperties>
class adjacency_list<OutEdgeListS, VertexListS, undirected_s, VertexProperties, EdgeProperties> {
 public:
  using check_allowed_vertex_list = typename adjacency_list_disallowed_vertex_list<VertexListS>::type;

  using self = adjacency_list<OutEdgeListS, VertexListS, undirected_s, VertexProperties, EdgeProperties>;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;

  using vertex_bundled = VertexProperties;
  using edge_bundled = EdgeProperties;

  using bidir_storage_type =
      typename adjlist_detail::adjlist_vertex_container<VertexListS, OutEdgeListS, bidirectional_s, VertexProperties,
                                                        EdgeProperties>;

  using vertex_descriptor = typename bidir_storage_type::vertex_descriptor;
  using vertices_size_type = typename bidir_storage_type::vertices_size_type;

  using bidir_edge_descriptor = typename bidir_storage_type::edge_descriptor;
  using edge_descriptor = container_detail::undir_edge_desc<bidir_edge_descriptor>;
  using edges_size_type = typename bidir_storage_type::edges_size_type;
  using degree_size_type = edges_size_type;

  using vertex_range = typename bidir_storage_type::vertex_range;

  using bidir_out_edge_range = typename bidir_storage_type::out_edge_range;
  using bidir_in_edge_range = typename bidir_storage_type::in_edge_range;
  using bidir_edge_range = typename bidir_storage_type::edge_range;

  using out_edge_range =
      adjlist_detail::adjlist_undir_ioerange<edge_descriptor, bidir_in_edge_range, bidir_out_edge_range>;
  using in_edge_range = out_edge_range;

  using edge_range = decltype(adjlist_detail::adjlist_undir_eiter_range(std::declval<bidir_edge_range>()));

  using adjacency_range = decltype(adjacency_range(std::declval<out_edge_range>(), std::declval<self>()));
  using inv_adjacency_range = decltype(inv_adjacency_range(std::declval<in_edge_range>(), std::declval<self>()));

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

  /**
   * Constructs an empty adjacency-list.
   */
  adjacency_list() : m_pack() {}

  ~adjacency_list() = default;

  adjacency_list(const self& rhs) : m_pack() { do_graph_deep_copy(*this, rhs); }
  self& operator=(const self& rhs) {
    if (this != &rhs) {
      do_graph_deep_copy(*this, rhs);
    }
    return *this;
  }

  adjacency_list(self&& rhs) noexcept : m_pack(std::move(rhs.m_pack)) {}
  self& operator=(self&& rhs) noexcept {
    m_pack = std::move(rhs.m_pack);
    return *this;
  }

  /**
   * Swaps the adjacency-list with another.
   */
  void swap(self& rhs) { m_pack.swap(rhs.m_pack); }

  /**
   * Clears the adjacency-list of all vertices and edges.
   */
  void clear() { m_pack.clear(); }

  /**
   * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
   * \param v The vertex descriptor of the sought-after vertex-property.
   * \return The vertex-property, by reference, associated to the given vertex descriptor.
   */
  vertex_property_type& operator[](vertex_descriptor v) { return m_pack.get_stored_vertex(v).data; }
  /**
   * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
   * \param v The vertex descriptor of the sought-after vertex-property.
   * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
   */
  const vertex_property_type& operator[](vertex_descriptor v) const { return m_pack.get_stored_vertex(v).data; }
  /**
   * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
   * \param e The edge descriptor of the sought-after edge-property.
   * \return The edge-property, by reference, associated to the given edge descriptor.
   */
  edge_property_type& operator[](const edge_descriptor& e) {
    return m_pack.get_stored_edge(bidir_edge_descriptor(e)).data;
  }
  /**
   * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
   * \param e The edge descriptor of the sought-after edge-property.
   * \return The edge-property, by const-reference, associated to the given edge descriptor.
   */
  const edge_property_type& operator[](const edge_descriptor& e) const {
    return m_pack.get_stored_edge(bidir_edge_descriptor(e)).data;
  }
};

#define BGL_ADJACENCY_LIST_UNDIR_ARGS \
  typename OutEdgeListS, typename VertexListS, typename VertexProperties, typename EdgeProperties
#define BGL_ADJACENCY_LIST_UNDIR \
  adjacency_list<OutEdgeListS, VertexListS, undirected_s, VertexProperties, EdgeProperties>

/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto source(const typename BGL_ADJACENCY_LIST_UNDIR::edge_descriptor& e, const BGL_ADJACENCY_LIST_UNDIR& g) {
  if (e.is_reversed) {
    return g.m_pack.get_stored_edge(e).target;
  }
  return e.source;
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto target(const typename BGL_ADJACENCY_LIST_UNDIR::edge_descriptor& e, const BGL_ADJACENCY_LIST_UNDIR& g) {
  if (e.is_reversed) {
    return e.source;
  }
  return g.m_pack.get_stored_edge(e).target;
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto out_edges(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_UNDIR& g) {
  return out_edge_range(true, g.m_pack.in_edges(v), g.m_pack.out_edges(v));
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
std::size_t out_degree(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.get_out_degree(v) + g.m_pack.get_in_degree(v);
}

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto in_edges(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_UNDIR& g) {
  return in_edge_range(false, g.m_pack.in_edges(v), g.m_pack.out_edges(v));
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
std::size_t in_degree(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
std::size_t degree(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
}

/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto edges(const BGL_ADJACENCY_LIST_UNDIR& g) {
  return adjlist_detail::adjlist_undir_eiter_range(g.m_pack.edges());
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto num_edges(const BGL_ADJACENCY_LIST_UNDIR& g) {
  return g.m_pack.num_edges();
}

/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto edge(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
          typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, const BGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
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

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto add_edge(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
              typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v, BGL_ADJACENCY_LIST_UNDIR& g) {
  using EdgeBundled = typename BGL_ADJACENCY_LIST_UNDIR::edge_bundled;
  using Edge = typename BGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  auto [be, be_found] = g.m_pack.add_edge(u, v, EdgeBundled());
  return std::pair(Edge(be), be_found);
}

/*******************************************************************************************
 *                  MutablePropertyGraph concept
 ******************************************************************************************/

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto add_edge(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
              typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v,
              const typename BGL_ADJACENCY_LIST_UNDIR::edge_property_type& ep, BGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  auto [be, be_found] = g.m_pack.add_edge(u, v, ep);
  return std::pair(Edge(be), be_found);
}

template <BGL_ADJACENCY_LIST_UNDIR_ARGS>
auto add_edge(typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor u,
              typename BGL_ADJACENCY_LIST_UNDIR::vertex_descriptor v,
              typename BGL_ADJACENCY_LIST_UNDIR::edge_property_type&& ep, BGL_ADJACENCY_LIST_UNDIR& g) {
  using Edge = typename BGL_ADJACENCY_LIST_UNDIR::edge_descriptor;
  auto [be, be_found] = g.m_pack.add_edge(u, v, std::move(ep));
  return std::pair(Edge(be), be_found);
}

#undef BGL_ADJACENCY_LIST_UNDIR_ARGS
#undef BGL_ADJACENCY_LIST_UNDIR

}  // namespace bagl

#endif  // BAGL_BAGL_ADJACENCY_LIST_H_
