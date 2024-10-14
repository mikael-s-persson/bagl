// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_LINKED_TREE_H_
#define BAGL_BAGL_LINKED_TREE_H_

#include <utility>

#include "bagl/adjacency_range.h"
#include "bagl/detail/ltree_containers.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/tree_traits.h"

namespace bagl {

// This library provides a class that implements a (doubly-)linked tree structure. This is a
// classic tree implementation in which each node contain a list of edges to its children,
// and a link to its parent.

// This traits class template is used to obtain the types (and meta-values) that describe
// the basic types used in a linked-tree with the given out-edge-list and vertex-list storage
// policies. This traits class is useful to obtain type-erased (or type-agnostic) vertex and edge
// descriptors. Note, this traits class is essentially the linked-tree equivalent of the BGL adjacency_list_traits
// class.
template <typename OutEdgeListS = vec_s, typename VertexListS = vec_s, typename DirectedS = directed_s>
struct linked_tree_traits {
  /** This meta-value tells if the edges are bidirectional, or not. */
  static constexpr bool is_bidir = DirectedS::is_bidir;
  /** This meta-value tells if the edges are directional, or not. */
  static constexpr bool is_directed = DirectedS::is_directed;

  /** This tag gives the edges' directional categorization. */
  using directed_category =
      std::conditional_t<is_bidir, bidirectional_tag, std::conditional_t<is_directed, directed_tag, undirected_tag>>;

  /** This meta-value tells if the parallel edges are allowed, or not. */
  using edge_parallel_category = typename container_detail::parallel_edge_traits<OutEdgeListS>::type;

  using traversal_category = adjlist_detail::adjlist_traversal_tag<DirectedS>;

  /** This meta-value tells if the vertex storage is random-access, or not. */
  using vertex_rand_access = typename container_detail::is_random_access<VertexListS>::type;
  using is_rand_access = vertex_rand_access;

  /** This meta-value tells if the vertex storage is random-access, or not. */
  using edge_rand_access = typename container_detail::is_random_access<OutEdgeListS>::type;

  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
};

struct linked_tree_tag {};

template <typename VertexListS>
struct linked_tree_disallowed_vertex_list {
  using type = void;
};

template <>
struct linked_tree_disallowed_vertex_list<set_s> {};
template <>
struct linked_tree_disallowed_vertex_list<multiset_s> {};
template <>
struct linked_tree_disallowed_vertex_list<unordered_set_s> {};
template <>
struct linked_tree_disallowed_vertex_list<unordered_multiset_s> {};

template <typename DirectedS>
struct linked_tree_disallowed_undirected {
  using type = void;
};

template <>
struct linked_tree_disallowed_undirected<undirected_s> {};

/**
 * This class implements a D-Ary Breadth-first tree that is tailored
 * to store elements of a tree as if their were inserted in a breadth-first manner. This type
 * of tree structure is good for both breadth-first search and depth-first search because
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree
 * should be kept balanced, and this implementation assumes that. The storage pattern is
 * similar to a binary heap tree-structure.
 * \tparam OutEdgeList A type tag to choose the storage policy for the out-edge lists.
 * \tparam VertexList A type tag to choose the storage policy for the vertices.
 * \tparam DirectedS A type tag to choose the directionality of the edges.
 * \tparam VertexProperties A POD type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A POD type to be attached to each edge in the tree.
 */
template <typename OutEdgeListS = vec_s, typename VertexListS = vec_s, typename DirectedS = directed_s,
          typename VertexProperties = no_property, typename EdgeProperties = no_property,
          typename GraphProperties = no_property>
class linked_tree {
 public:
  using check_allowed_vertex_list = typename linked_tree_disallowed_vertex_list<VertexListS>::type;
  using check_allowed_directionality = typename linked_tree_disallowed_undirected<DirectedS>::type;

  using self = linked_tree<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties, GraphProperties>;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;
  using graph_property_type = GraphProperties;

  using vertex_bundled = lookup_one_property_t<VertexProperties, vertex_bundle_t>;
  using edge_bundled = lookup_one_property_t<EdgeProperties, edge_bundle_t>;
  using graph_bundled = lookup_one_property_t<GraphProperties, graph_bundle_t>;

  using storage_type = typename adjlist_detail::ltree_vertex_container<VertexListS, OutEdgeListS, DirectedS,
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

  using Traits = linked_tree_traits<OutEdgeListS, VertexListS, DirectedS>;

  using vertex_rand_access = typename Traits::vertex_rand_access;
  using edge_rand_access = typename Traits::edge_rand_access;
  using directed_category = typename Traits::directed_category;
  using edge_parallel_category = typename Traits::edge_parallel_category;
  using traversal_category = typename Traits::traversal_category;

  using graph_tag = linked_tree_tag;

  /**
   * This static member function outputs the null-vertex (invalid vertex descriptor).
   * \return A null-vertex descriptor (invalid vertex descriptor).
   */
  static vertex_descriptor null_vertex() { return container_detail::null_desc<vertex_descriptor>::value(); };

  /**
   * This static member function outputs the null-edge (invalid edge descriptor).
   * \return A null-edge descriptor (invalid edge descriptor).
   */
  static edge_descriptor null_edge() { return container_detail::null_desc<edge_descriptor>::value(); };

  // private:
  storage_type m_pack;
  graph_property_type m_graph_prop;

  /**
   * Constructs an empty linked-tree.
   */
  linked_tree() : m_pack(), m_graph_prop(){};

  ~linked_tree() = default;

 private:
  void do_deep_copy_from(const self& rhs);

 public:
  /**
   * Constructs a linked-tree as a copy of the given tree.
   */
  linked_tree(const self& rhs) : m_pack(), m_graph_prop(rhs.m_graph_prop) { do_deep_copy_from(rhs); };
  /**
   * Assigns the linked-tree as a copy of the given tree.
   */
  self& operator=(const self& rhs) {
    if (this != &rhs) {
      m_graph_prop = rhs.m_graph_prop;
      do_deep_copy_from(rhs);
    }
    return *this;
  };

  /**
   * Constructs a linked-tree by moving the given tree into it.
   */
  linked_tree(self&& rhs) noexcept : m_pack(std::move(rhs.m_pack)), m_graph_prop(std::move(rhs.m_graph_prop)){};
  /**
   * Assigns a linked-tree by moving the given tree into it.
   */
  linked_tree& operator=(self&& rhs) noexcept {
    if (this != rhs) {
      m_graph_prop = std::move(rhs.m_graph_prop);
      m_pack = std::move(rhs.m_pack);
    }
    return *this;
  };

  /**
   * Standard swap function.
   */
  void swap(self& rhs) { m_pack.swap(rhs.m_pack); };

  /**
   * Returns the size of the tree (the number of vertices it contains).
   * \return The size of the tree (the number of vertices it contains).
   */
  [[nodiscard]] std::size_t size() const { return m_pack.size(); };

  /**
   * Checks if the tree is empty.
   * \return True if the tree is empty.
   */
  [[nodiscard]] bool empty() const { return size() == 0; };

  /**
   * Returns the maximum vertex capacity of the tree (the number of vertices it can contain).
   * \return The maximum vertex capacity of the tree (the number of vertices it can contain).
   */
  [[nodiscard]] std::size_t capacity() const { return m_pack.m_vertices.capacity(); };

  /**
   * Returns the depth of the tree.
   * \note This operation must recurse through all the branches of the tree (depth-first), and is
   * thus an expensive operation (linear-time w.r.t. the number of vertices, and linear-memory (stack)
   * w.r.t. the depth of tree).
   * \return The depth of the tree.
   */
  [[nodiscard]] std::size_t depth() const {
    if (m_pack.m_root != null_vertex()) {
      return m_pack.get_depth(m_pack.m_root);
    }
    return 0;
  };

  /**
   * Clears the tree of all vertices and edges.
   */
  void clear() { m_pack.clear(); };

  // Indexing operator. Returns a reference to the vertex-bundle associated to the given vertex descriptor.
  vertex_bundled& operator[](vertex_descriptor v) {
    return get_property_value(m_pack.get_stored_vertex(v).data, vertex_bundle);
  };
  // Indexing operator. Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
  const vertex_bundled& operator[](vertex_descriptor v) const {
    return get_property_value(m_pack.get_stored_vertex(v).data, vertex_bundle);
  };
  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  edge_bundled& operator[](const edge_descriptor& e) {
    return get_property_value(m_pack.get_stored_edge_property(e), edge_bundle);
  };
  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  const edge_bundled& operator[](const edge_descriptor& e) const {
    return get_property_value(m_pack.get_stored_edge_property(e), edge_bundle);
  };
  // Indexing operator. Returns a reference to the graph-bundle associated to the graph.
  graph_bundled& operator[](graph_bundle_t /*unused*/) { return get_property_value(m_graph_prop, graph_bundle); };
  // Indexing operator. Returns a const-reference to the graph-bundle associated to the graph.
  const graph_bundled& operator[](graph_bundle_t /*unused*/) const {
    return get_property_value(m_graph_prop, graph_bundle);
  };

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
 * This is the tree-storage specifier for a linked-tree of the given edge and vertex storage policies.
 */
template <typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct linked_tree_storage {};

template <typename VertexProperties, typename EdgeProperties, typename OutEdgeListS, typename VertexListS,
          typename DirectedS>
struct tree_storage<VertexProperties, EdgeProperties, linked_tree_storage<OutEdgeListS, VertexListS, DirectedS>> {
  using type = linked_tree<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties>;
};

template <typename OutEdgeListS, typename VertexListS, typename DirectedS>
struct tree_storage_traits<linked_tree_storage<OutEdgeListS, VertexListS, DirectedS>>
    : linked_tree_traits<OutEdgeListS, VertexListS, DirectedS> {};

#define BAGL_LINKED_TREE_ARGS                                                                                          \
  typename OutEdgeListS, typename VertexListS, typename DirectedS, typename VertexProperties, typename EdgeProperties, \
      typename GraphProperties
#define BAGL_LINKED_TREE \
  linked_tree<OutEdgeListS, VertexListS, DirectedS, VertexProperties, EdgeProperties, GraphProperties>

/**
 * Standard swap function. Swaps the contents of two objects.
 * \param lhs The left-hand-side of the swap.
 * \param rhs The right-hand-side of the swap.
 */
template <BAGL_LINKED_TREE_ARGS>
void swap(BAGL_LINKED_TREE& lhs, BAGL_LINKED_TREE& rhs) {
  lhs.swap(rhs);
};

/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

/**
 * Returns the source vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The source vertex of the given edge descriptor.
 */
template <BAGL_LINKED_TREE_ARGS>
auto source(typename BAGL_LINKED_TREE::edge_descriptor e, const BAGL_LINKED_TREE& /*unused*/) {
  return e.source;
};

/**
 * Returns the target vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The target vertex of the given edge descriptor.
 */
template <BAGL_LINKED_TREE_ARGS>
auto target(typename BAGL_LINKED_TREE::edge_descriptor e, const BAGL_LINKED_TREE& g) {
  return g.m_pack.get_stored_edge_target(e);
};

/**
 * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the out-edges of a given vertex descriptor.
 */
template <BAGL_LINKED_TREE_ARGS>
auto out_edges(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.m_pack.out_edges(v);
};

/**
 * Returns the out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The out-degree of the given vertex descriptor.
 */
template <BAGL_LINKED_TREE_ARGS>
std::size_t out_degree(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.m_pack.get_out_degree(v);
};

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

/**
 * Returns the edge iterator range for the in-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the in-edges of a given vertex descriptor.
 */
template <BAGL_LINKED_TREE_ARGS>
auto in_edges(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.m_pack.in_edges(v);
};

/**
 * Returns the in-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template <BAGL_LINKED_TREE_ARGS>
std::size_t in_degree(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.m_pack.get_in_degree(v);
};

/**
 * Returns the in-degree plus out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree plus out-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template <BAGL_LINKED_TREE_ARGS>
std::size_t degree(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.m_pack.get_in_degree(v) + g.m_pack.get_out_degree(v);
};

/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for all the vertices of the tree.
 * \param g The graph.
 * \return The vertex iterator range for all the vertices of the tree.
 */
template <BAGL_LINKED_TREE_ARGS>
auto vertices(const BAGL_LINKED_TREE& g) {
  return g.m_pack.vertices();
};

/**
 * Returns the size of the tree (the number of vertices it contains).
 * \param g The graph.
 * \return The size of the tree (the number of vertices it contains).
 */
template <BAGL_LINKED_TREE_ARGS>
auto num_vertices(const BAGL_LINKED_TREE& g) {
  return g.m_pack.size();
};

/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the edge iterator range for all the edges of the tree.
 * \param g The graph.
 * \return The edge iterator range for all the edges of the tree.
 */
template <BAGL_LINKED_TREE_ARGS>
auto edges(const BAGL_LINKED_TREE& g) {
  return g.m_pack.edges();
};

/**
 * Returns the number of edges in the tree.
 * \param g The graph.
 * \return The number of edges in the tree.
 */
template <BAGL_LINKED_TREE_ARGS>
typename BAGL_LINKED_TREE::edges_size_type num_edges(const BAGL_LINKED_TREE& g) {
  std::size_t tmp = g.m_pack.size();
  if (tmp > 0) {
    return tmp - 1;
  }
  return 0;
};

/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

template <BAGL_LINKED_TREE_ARGS>
auto adjacent_vertices(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return adjacency_range(out_edges(v, g), g);
}

/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

template <BAGL_LINKED_TREE_ARGS>
auto inv_adjacent_vertices(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return inv_adjacency_range(in_edges(v, g), g);
}

/***********************************************************************************************
 *                             AdjacencyMatrixConcept
 * ********************************************************************************************/

/**
 * Returns the edge descriptor for the edge between two given vertex descriptors.
 * \param u The vertex descriptor of the source vertex.
 * \param v The vertex descriptor of the target vertex.
 * \param g The graph.
 * \return The edge descriptor for the given vertex descriptor pair.
 */
template <BAGL_LINKED_TREE_ARGS>
auto edge(typename BAGL_LINKED_TREE::vertex_descriptor u, typename BAGL_LINKED_TREE::vertex_descriptor v,
          const BAGL_LINKED_TREE& g) {
  return g.m_pack.get_edge(u, v);
};

/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

/**
 * Returns the vertex-descriptor of the root of the tree.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_LINKED_TREE_ARGS>
auto tree_root(const BAGL_LINKED_TREE& g) {
  return g.m_pack.m_root;
};

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template <BAGL_LINKED_TREE_ARGS>
auto children(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return adjacent_vertices(v, g);
};

/***********************************************************************************************
 *                             BidirectionalTreeConcept
 * ********************************************************************************************/

/**
 * Returns the parent vertex of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The parent vertex of the given vertex descriptor (will be null_vertex() if it is the root (no parent)).
 */
template <BAGL_LINKED_TREE_ARGS>
auto parent(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.m_pack.get_parent(v);
};

/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_LINKED_TREE_ARGS>
auto create_root(BAGL_LINKED_TREE& g) {
  using VProp = typename BAGL_LINKED_TREE::vertex_property_type;
  if (g.m_pack.size()) {
    g.m_pack.clear();
  }
  g.m_pack.add_root_vertex(VProp());
  return g.m_pack.m_root;
};

/**
 * Adds a child vertex to the given parent vertex, and default-initializes the properties of
 * the newly created vertex and edge.
 * \param v The parent vertex to which a child will be added.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <BAGL_LINKED_TREE_ARGS>
auto add_child(typename BAGL_LINKED_TREE::vertex_descriptor v, BAGL_LINKED_TREE& g) {
  using VProp = typename BAGL_LINKED_TREE::vertex_property_type;
  using EProp = typename BAGL_LINKED_TREE::edge_property_type;
  return g.m_pack.add_child(v, VProp(), EProp());
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex.
 * \param v The root of the sub-tree to be removed.
 * \param g The graph.
 */
template <BAGL_LINKED_TREE_ARGS>
void clear_children(typename BAGL_LINKED_TREE::vertex_descriptor v, BAGL_LINKED_TREE& g) {
  g.m_pack.clear_children_impl(v);
};

/**
 * Removes a branch (sub-tree) starting from and including the given vertex.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param g The graph.
 */
template <BAGL_LINKED_TREE_ARGS>
void remove_branch(typename BAGL_LINKED_TREE::vertex_descriptor v, BAGL_LINKED_TREE& g) {
  g.m_pack.remove_branch_impl(v);
};

/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_LINKED_TREE_ARGS, typename VProp>
auto create_root(VProp&& vp, BAGL_LINKED_TREE& g) {
  if (g.m_pack.size()) {
    g.m_pack.clear();
  }
  g.m_pack.add_root_vertex(std::forward<VProp>(vp));
  return g.m_pack.m_root;
};

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
 * vertex to the given property value.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <BAGL_LINKED_TREE_ARGS, typename VProp>
auto add_child(typename BAGL_LINKED_TREE::vertex_descriptor v, VProp&& vp, BAGL_LINKED_TREE& g) {
  using EProp = typename BAGL_LINKED_TREE::edge_property_type;
  return g.m_pack.add_child(v, std::forward<VProp>(vp), EProp());
};

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
 * vertex and edge to the given property values.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param ep The property value to be moved into the newly created edge.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <BAGL_LINKED_TREE_ARGS, typename VProp, typename EProp>
auto add_child(typename BAGL_LINKED_TREE::vertex_descriptor v, VProp&& vp, EProp&& ep, BAGL_LINKED_TREE& g) {
  return g.m_pack.add_child(v, std::forward<VProp>(vp), std::forward<EProp>(ep));
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template <BAGL_LINKED_TREE_ARGS, typename OutputIter>
OutputIter clear_children(typename BAGL_LINKED_TREE::vertex_descriptor v, OutputIter it_out, BAGL_LINKED_TREE& g) {
  return g.m_pack.clear_children_impl(v, it_out);
};

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template <BAGL_LINKED_TREE_ARGS, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children(typename BAGL_LINKED_TREE::vertex_descriptor v, VertexOIter vit_out,
                                                 EdgeOIter eit_out, BAGL_LINKED_TREE& g) {
  return g.m_pack.clear_children_impl(v, vit_out, eit_out);
};

/**
 * Removes a branch (sub-tree) starting from and including the given vertex, while
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 * \note The first vertex-property to figure in the output range is that of the vertex v.
 */
template <BAGL_LINKED_TREE_ARGS, typename OutputIter>
OutputIter remove_branch(typename BAGL_LINKED_TREE::vertex_descriptor v, OutputIter it_out, BAGL_LINKED_TREE& g) {
  return g.m_pack.remove_branch_impl(v, it_out);
};

/**
 * Removes a branch (sub-tree) starting from and including the given vertex, while
 * recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 * \note The first vertex-property to figure in the output range is that of the vertex v.
 */
template <BAGL_LINKED_TREE_ARGS, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch(typename BAGL_LINKED_TREE::vertex_descriptor v, VertexOIter vit_out,
                                                EdgeOIter eit_out, BAGL_LINKED_TREE& g) {
  return g.m_pack.remove_branch_impl(v, vit_out, eit_out);
};

/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

// Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
template <BAGL_LINKED_TREE_ARGS>
const auto& get(const BAGL_LINKED_TREE& g, typename BAGL_LINKED_TREE::vertex_descriptor v) {
  return g[v];
}
template <BAGL_LINKED_TREE_ARGS>
auto& get(BAGL_LINKED_TREE& g, typename BAGL_LINKED_TREE::vertex_descriptor v) {
  return g[v];
}

// Returns a const-reference to the edge-bundle associated to the given edge descriptor.
template <BAGL_LINKED_TREE_ARGS>
const auto& get(const BAGL_LINKED_TREE& g, typename BAGL_LINKED_TREE::edge_descriptor e) {
  return g[e];
}
template <BAGL_LINKED_TREE_ARGS>
auto& get(BAGL_LINKED_TREE& g, typename BAGL_LINKED_TREE::edge_descriptor e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <BAGL_LINKED_TREE_ARGS>
const auto& get(const BAGL_LINKED_TREE& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}
template <BAGL_LINKED_TREE_ARGS>
auto& get(BAGL_LINKED_TREE& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <BAGL_LINKED_TREE_ARGS, typename VProp>
void put(BAGL_LINKED_TREE& g, typename BAGL_LINKED_TREE::vertex_descriptor v, VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-bundle associated to the given edge descriptor.
template <BAGL_LINKED_TREE_ARGS, typename EProp>
void put(BAGL_LINKED_TREE& g, typename BAGL_LINKED_TREE::edge_descriptor e, EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Sets the graph-bundle associated to the graph.
template <BAGL_LINKED_TREE_ARGS, typename GProp>
void put(BAGL_LINKED_TREE& g, graph_bundle_t /*unused*/, GProp&& value) {
  g[graph_bundle] = std::forward<GProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_LINKED_TREE_ARGS>
auto& get_property(typename BAGL_LINKED_TREE::vertex_descriptor v, BAGL_LINKED_TREE& g) {
  return g.get_property(v);
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_LINKED_TREE_ARGS>
const auto& get_property(typename BAGL_LINKED_TREE::vertex_descriptor v, const BAGL_LINKED_TREE& g) {
  return g.get_property(v);
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <BAGL_LINKED_TREE_ARGS>
auto& get_property(typename BAGL_LINKED_TREE::edge_descriptor e, BAGL_LINKED_TREE& g) {
  return g.get_property(e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <BAGL_LINKED_TREE_ARGS>
const auto& get_property(typename BAGL_LINKED_TREE::edge_descriptor e, const BAGL_LINKED_TREE& g) {
  return g.get_property(e);
}

// Returns a reference to the graph-property associated to the graph.
template <BAGL_LINKED_TREE_ARGS>
auto& get_property(BAGL_LINKED_TREE& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Returns a const-reference to the graph-property associated to the graph.
template <BAGL_LINKED_TREE_ARGS>
const auto& get_property(const BAGL_LINKED_TREE& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Handle graph property tags, also handles graph_bundle_t.

template <BAGL_LINKED_TREE_ARGS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>, lookup_one_property_t<GraphProperties, Tag>&>
get_property(BAGL_LINKED_TREE& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_LINKED_TREE_ARGS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>,
                 const lookup_one_property_t<GraphProperties, Tag>&>
get_property(const BAGL_LINKED_TREE& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_LINKED_TREE_ARGS, typename T, typename Bundle>
struct property_map<BAGL_LINKED_TREE, T Bundle::*> {
  using non_const_bundle = std::remove_cv_t<Bundle>;
  using non_const_t = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v =
      std::is_convertible_v<typename BAGL_LINKED_TREE::vertex_bundled*, non_const_bundle*>;
  static constexpr bool is_edge_bundle_v =
      std::is_convertible_v<typename BAGL_LINKED_TREE::edge_bundled*, non_const_bundle*>;
  using tag_type = std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                                      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<non_const_t, BAGL_LINKED_TREE, tag_type>;
  using const_type = bundle_member_property_map<const non_const_t, const BAGL_LINKED_TREE, tag_type>;
};

template <BAGL_LINKED_TREE_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, BAGL_LINKED_TREE& g) {
  return typename property_map<BAGL_LINKED_TREE, T Bundle::*>::type(&g, p);
}

template <BAGL_LINKED_TREE_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, const BAGL_LINKED_TREE& g) {
  return typename property_map<BAGL_LINKED_TREE, T Bundle::*>::const_type(&g, p);
}

template <BAGL_LINKED_TREE_ARGS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p, const BAGL_LINKED_TREE& g, const Key& k) {
  return (g[k]).*p;
}

template <BAGL_LINKED_TREE_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_LINKED_TREE& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct linked_tree_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;

    using type = tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type = tagged_in_property_property_map<const value_type, const Graph, Tag>;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<linked_tree_tag> {
  using type = linked_tree_property_selector;
};

template <>
struct edge_property_selector<linked_tree_tag> {
  using type = linked_tree_property_selector;
};

template <BAGL_LINKED_TREE_ARGS, typename Property>
auto get(Property p, BAGL_LINKED_TREE& g) {
  using Map = typename property_map<BAGL_LINKED_TREE, Property>::type;
  return Map(&g, p);
}

template <BAGL_LINKED_TREE_ARGS, typename Property>
auto get(Property p, const BAGL_LINKED_TREE& g) {
  using Map = typename property_map<BAGL_LINKED_TREE, Property>::const_type;
  return Map(&g, p);
}

template <BAGL_LINKED_TREE_ARGS, typename Property, typename Key>
decltype(auto) get(Property p, const BAGL_LINKED_TREE& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_LINKED_TREE_ARGS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_LINKED_TREE& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_LINKED_TREE_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BAGL_LINKED_TREE& g, const Key& k, Value&& val) {
  get_property_value(g.get_property(k), p) = std::forward<Value>(val);
}

template <BAGL_LINKED_TREE_ARGS>
void BAGL_LINKED_TREE::do_deep_copy_from(const BAGL_LINKED_TREE& rhs) {
  using Vertex = typename BAGL_LINKED_TREE::vertex_descriptor;
  using TaskType = std::pair<Vertex, Vertex>;

  this->m_pack.clear();

  // NOTE: Do the copy in breadth-first order because it leads to the best memory traversal patterns later on.
  std::queue<TaskType> bft_queue;
  TaskType cur;
  cur.first = tree_root(rhs);
  cur.second = create_root(rhs[cur.first], *this);
  bft_queue.push(cur);

  while (!bft_queue.empty()) {
    cur = bft_queue.front();
    bft_queue.pop();
    for (auto e : out_edges(cur.first, rhs)) {
      auto [new_u, new_e, suc] = add_child(cur.second, rhs[target(e, rhs)], rhs[e], *this);
      bft_queue.push(TaskType(target(e, rhs), new_u));
    };
  };
};

#undef BAGL_LINKED_TREE_ARGS
#undef BAGL_LINKED_TREE

};  // namespace bagl

#endif  // BAGL_BAGL_LINKED_TREE_H_
