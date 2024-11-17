// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_VEBL_D_ARY_TREE_H_
#define BAGL_BAGL_VEBL_D_ARY_TREE_H_

#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/detail/container_generators.h"
#include "bagl/detail/vebl_tree_ranges.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/tree_traits.h"

namespace bagl {

// This library provides a class that implements a D-Ary von Emde Boas Layout tree that stores
// elements in a contiguous array in a cache-friendly recursive layout. This type
// of tree structure is good for both breadth-first search and depth-first search because
// of locality of reference issues. Ideally, for the least amount of wasted memory, the tree
// should be kept balanced, and this implementation assumes that.

struct vebl_d_ary_tree_tag {};

/**
 * This class implements a D-Ary von Emde Boas Layout tree that stores
 * elements in a contiguous array in a cache-friendly recursive layout. This type
 * of tree structure is good for both breadth-first search and depth-first search because
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree
 * should be kept balanced, and this implementation assumes that.
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <std::size_t Arity = 2, typename VertexProperties = no_property, typename EdgeProperties = no_property,
          typename GraphProperties = no_property>
class vebl_d_ary_tree {
 public:
  using self = vebl_d_ary_tree<Arity, VertexProperties, EdgeProperties>;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;
  using graph_property_type = GraphProperties;

  using vertex_bundled = lookup_one_property_t<VertexProperties, vertex_bundle_t>;
  using edge_bundled = lookup_one_property_t<EdgeProperties, edge_bundle_t>;
  using graph_bundled = lookup_one_property_t<GraphProperties, graph_bundle_t>;

  using value_type = bfl_detail::bfltree_value_type<vertex_property_type, edge_property_type>;

  using container_type = std::vector<value_type>;

  using vertex_descriptor = std::size_t;
  using edge_descriptor = bfl_detail::bfltree_edge_desc;

  using vertices_size_type = std::size_t;
  using edges_size_type = vertices_size_type;
  using degree_size_type = vertices_size_type;

  /**
   * This static member function outputs the null-vertex (invalid vertex descriptor).
   * \return A null-vertex descriptor (invalid vertex descriptor).
   */
  static vertex_descriptor null_vertex() { return std::numeric_limits<std::size_t>::max(); }

  static edge_descriptor null_edge() { return edge_descriptor{null_vertex()}; }

  static auto make_unsafe_vertex_range(std::size_t first, std::size_t last) { return std::views::iota(first, last); }

  using vertex_validity = bfl_detail::vebltree_vertex_validity<container_type, Arity>;
  auto make_valid_vertex_range(std::size_t first, std::size_t last) const {
    return make_unsafe_vertex_range(first, last) | std::views::filter(vertex_validity{&m_vertices, &m_depth_recs});
  }

  template <typename VertexRange>
  static auto view_vertex_range_as_edges(VertexRange v_rg) {
    return v_rg | std::views::transform([](std::size_t u) { return edge_descriptor{u}; });
  }

  using directed_category = directed_tag;
  using edge_parallel_category = disallow_parallel_edge_tag;

  struct traversal_category : virtual public incidence_graph_tag,
                              virtual public adjacency_graph_tag,
                              virtual public bidirectional_graph_tag,
                              virtual public vertex_list_graph_tag,
                              virtual public edge_list_graph_tag {};

  using graph_tag = vebl_d_ary_tree_tag;

  // private:
  container_type m_vertices;
  vertices_size_type m_vertex_count{0};
  graph_property_type m_graph_prop;

  bfl_detail::vebl_depth_records m_depth_recs;

  value_type& get_vertex_value_impl(vertex_descriptor v) {
    return m_vertices[bfl_detail::convert_bfl_to_vebl<Arity>(v, m_depth_recs)];
  }
  const value_type& get_vertex_value_impl(vertex_descriptor v) const {
    return m_vertices[bfl_detail::convert_bfl_to_vebl<Arity>(v, m_depth_recs)];
  }

  static vertex_descriptor get_child_impl(vertex_descriptor v, std::size_t edge_id) { return Arity * v + 1 + edge_id; }

  /**
   * Construct the D-ary BF-tree with a given reserved depth.
   * \param aDepth The depth of the graph to reserve space for.
   */
  explicit vebl_d_ary_tree(vertices_size_type aDepth = 0) : m_graph_prop() {
    m_depth_recs.T.push_back(0);
    m_depth_recs.B.push_back(1);
    m_depth_recs.D.push_back(0);
    for (vertices_size_type i = 0; i < aDepth; ++i) {
      bfl_detail::extend_vebl_depth_records<Arity>(m_depth_recs);
    }
    m_vertices.resize(bfl_detail::s_treesize<Arity>(aDepth));
  }

  /**
   * Checks if the tree is empty.
   * \return True if the tree is empty.
   */
  [[nodiscard]] bool empty() const { return m_vertex_count == 0; }

  /**
   * Returns the size of the tree (the number of vertices it contains).
   * \return The size of the tree (the number of vertices it contains).
   */
  [[nodiscard]] std::size_t size() const { return m_vertex_count; }

  /**
   * Returns the maximum vertex capacity of the tree (the number of vertices it can contain).
   * \return The maximum vertex capacity of the tree (the number of vertices it can contain).
   */
  [[nodiscard]] std::size_t capacity() const { return m_vertices.capacity(); }

  /**
   * Returns the depth of the tree.
   * \return The depth of the tree.
   */
  [[nodiscard]] std::size_t depth() const { return m_depth_recs.T.size(); }

  /**
   * Standard swap function.
   */
  void swap(self& rhs) {
    using std::swap;
    m_vertices.swap(rhs.m_vertices);
    swap(m_vertex_count, rhs.m_vertex_count);
    swap(m_graph_prop, rhs.m_graph_prop);
    m_depth_recs.T.swap(rhs.m_depth_recs.T);
    m_depth_recs.B.swap(rhs.m_depth_recs.B);
    m_depth_recs.D.swap(rhs.m_depth_recs.D);
  }

  /**
   * Clears the tree of all vertices and edges.
   */
  void clear() {
    m_vertices.clear();
    m_vertices.resize(1);
    m_depth_recs.T.clear();
    m_depth_recs.T.push_back(0);
    m_depth_recs.B.clear();
    m_depth_recs.B.push_back(1);
    m_depth_recs.D.clear();
    m_depth_recs.D.push_back(0);
    m_vertex_count = 0;
  }

  // Indexing operator. Returns a reference to the vertex-bundle associated to the given vertex descriptor.
  vertex_bundled& operator[](const vertex_descriptor& v) {
    return get_property_value(get_vertex_value_impl(v).vertex(), vertex_bundle);
  }
  // Indexing operator. Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
  const vertex_bundled& operator[](const vertex_descriptor& v) const {
    return get_property_value(get_vertex_value_impl(v).vertex(), vertex_bundle);
  }
  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  edge_bundled& operator[](const edge_descriptor& e) {
    return get_property_value(get_vertex_value_impl(e.target_vertex).edge(), edge_bundle);
  }
  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  const edge_bundled& operator[](const edge_descriptor& e) const {
    return get_property_value(get_vertex_value_impl(e.target_vertex).edge(), edge_bundle);
  }
  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  graph_bundled& operator[](graph_bundle_t /*unused*/) { return get_property_value(m_graph_prop, graph_bundle); }
  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  const graph_bundled& operator[](graph_bundle_t /*unused*/) const {
    return get_property_value(m_graph_prop, graph_bundle);
  }

  // Get a reference to the vertex-property associated to the given vertex descriptor.
  auto& get_property(vertex_descriptor v) { return get_vertex_value_impl(v).vertex(); }
  // Get a const-reference to the vertex-property associated to the given vertex descriptor.
  const auto& get_property(vertex_descriptor v) const { return get_vertex_value_impl(v).vertex(); }

  // Get a reference to the edge-property associated to the given edge descriptor.
  auto& get_property(const edge_descriptor& e) { return get_vertex_value_impl(e.target_vertex).edge(); }
  // Get a const-reference to the edge-property associated to the given edge descriptor.
  const auto& get_property(const edge_descriptor& e) const { return get_vertex_value_impl(e.target_vertex).edge(); }

  // Get a reference to the graph-property associated to the graph.
  auto& get_property(graph_all_t /*unused*/) { return m_graph_prop; }
  // Get a const-reference to the graph-property associated to the graph.
  const auto& get_property(graph_all_t /*unused*/) const { return m_graph_prop; }

  /**
   * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
   * vertex and edge to the given property values.
   * \param v The parent vertex to which a child will be added.
   * \param vp The property value for the newly created vertex.
   * \param ep The property value for the newly created edge.
   * \return A pair consisting of the newly created vertex and edge (descriptors).
   */
  template <typename VP, typename EP>
  std::tuple<vertex_descriptor, edge_descriptor, bool> add_child(vertex_descriptor v, VP&& vp, EP&& ep) {
    if ((v >= m_vertices.size()) || !bfl_detail::bfltree_is_vertex_valid(get_vertex_value_impl(v))) {
      return {null_vertex(), null_edge(), false};
    }
    vertex_descriptor result = get_child_impl(v, 0);
    for (; result < get_child_impl(v, Arity); ++result) {
      if ((result >= m_vertices.size()) || !bfl_detail::bfltree_is_vertex_valid(get_vertex_value_impl(result))) {
        break;
      }
    }
    if (result == get_child_impl(v, Arity)) {
      return {null_vertex(), null_edge(), false};
    }

    if (result >= m_vertices.size()) {
      bfl_detail::extend_vebl_storage<Arity>(m_vertices, m_depth_recs);
    }

    value_type& result_prop = get_vertex_value_impl(result);
    result_prop.out_degree = 0;
    result_prop.vertex() = vertex_property_type{std::forward<VP>(vp)};
    result_prop.edge() = edge_property_type{std::forward<EP>(ep)};
    ++(get_vertex_value_impl(v).out_degree);
    ++m_vertex_count;
    return {result, edge_descriptor(result), true};
  }

  template <typename VertexOIter, typename EdgeOIter>
  void clear_children_impl(vertex_descriptor v, VertexOIter& vit_out, EdgeOIter& eit_out) {
    // this traversal order is intentional (traverse pre-order depth-first, and
    // delay removal of empty tail elements as much as possible, such that it is only required once).
    for (std::size_t i = 0; i < Arity; ++i) {
      vertex_descriptor next_v = get_child_impl(v, i);
      if (next_v >= m_vertices.size()) {
        break;
      }
      value_type& next_v_prop = get_vertex_value_impl(next_v);
      if (!bfl_detail::bfltree_is_vertex_valid(next_v_prop)) {
        continue;
      }
      --m_vertex_count;
      *(vit_out++) = std::move(next_v_prop.vertex());
      *(eit_out++) = std::move(next_v_prop.edge());
      clear_children_impl(next_v, vit_out, eit_out);
    }
    get_vertex_value_impl(v).out_degree = std::numeric_limits<std::size_t>::max();
    if (v != 0) {
      // if the node is not the root one, then update the out-degree of the parent node:
      get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
    }
  }

  /**
   * Removes a branch (sub-tree) starting from but excluding the given vertex, while
   * recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
   * \param v The root of the sub-tree to be removed.
   * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
   * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
   * \return The output-iterator after the collection of all the removed vertices.
   */
  template <typename VertexOIter, typename EdgeOIter>
  std::pair<VertexOIter, EdgeOIter> clear_children(vertex_descriptor v, VertexOIter vit_out, EdgeOIter eit_out) {
    if (v >= m_vertices.size()) {
      return {vit_out, eit_out};  // vertex is already deleted.
    }
    value_type& v_prop = get_vertex_value_impl(v);
    if (!bfl_detail::bfltree_is_vertex_valid(v_prop)) {
      return {vit_out, eit_out};  // vertex is already deleted.
    }
    clear_children_impl(v, vit_out, eit_out);
    return {vit_out, eit_out};
  }

  /**
   * Removes a branch (sub-tree) starting from but excluding the given vertex, while
   * recording the vertex-properties of all the removed vertices into an output-iterator.
   * \param v The root of the sub-tree to be removed.
   * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
   * \return The output-iterator after the collection of all the removed vertices.
   */
  template <typename VertexOIter>
  VertexOIter clear_children(vertex_descriptor v, VertexOIter vit_out) {
    return clear_children(v, vit_out, container_detail::ignore_output_iter{}).first;
  }

  /**
   * Removes a branch (sub-tree) starting from but excluding the given vertex.
   * \param v The root of the sub-tree to be removed.
   */
  void clear_children(vertex_descriptor v) {
    clear_children(v, container_detail::ignore_output_iter{}, container_detail::ignore_output_iter{});
  }

  /**
   * Removes a branch (sub-tree) starting from and including the given vertex, while
   * recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
   * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
   * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
   * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
   * \return The output-iterator after the collection of all the removed vertices.
   * \note The first vertex-property to figure in the output range is that of the vertex v.
   */
  template <typename VertexOIter, typename EdgeOIter>
  std::pair<VertexOIter, EdgeOIter> remove_branch(vertex_descriptor v, VertexOIter vit_out, EdgeOIter eit_out) {
    if (v >= m_vertices.size()) {
      return {vit_out, eit_out};  // vertex is already deleted.
    }
    value_type& v_prop = get_vertex_value_impl(v);
    if (!bfl_detail::bfltree_is_vertex_valid(v_prop)) {
      return {vit_out, eit_out};  // vertex is already deleted.
    }
    if (v == 0) {
      *(eit_out++) = edge_property_type();
    } else {
      // in-edge:  u = (v - 1) / Arity;  e_id = (v - 1) % Arity;
      *(eit_out++) = std::move(v_prop.edge());
      // if the node is not the root one, then update the out-degree of the parent node:
      get_vertex_value_impl((v - 1) / Arity).out_degree -= 1;
    }
    --m_vertex_count;
    *(vit_out++) = std::move(v_prop.vertex());
    clear_children_impl(v, vit_out, eit_out);
    return {vit_out, eit_out};
  }

  /**
   * Removes a branch (sub-tree) starting from and including the given vertex, while
   * recording the vertex-properties of all the removed vertices into an output-iterator.
   * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
   * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
   * \return The output-iterator after the collection of all the removed vertices.
   * \note The first vertex-property to figure in the output range is that of the vertex v.
   */
  template <typename VertexOIter>
  VertexOIter remove_branch(vertex_descriptor v, VertexOIter vit_out) {
    return remove_branch(v, vit_out, container_detail::ignore_output_iter{}).first;
  }

  /**
   * Removes a branch (sub-tree) starting from and including the given vertex.
   * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
   */
  void remove_branch(vertex_descriptor v) {
    remove_branch(v, container_detail::ignore_output_iter{}, container_detail::ignore_output_iter{});
  }

  /**
   * Creates a root for the tree (clears it if not empty), and moves the given vertex-property into it.
   * \param vp The vertex-property to move into the newly created root vertex.
   * \return The vertex-descriptor of the root of the tree.
   */
  template <typename VProp>
  vertex_descriptor create_root_vertex(VProp&& vp) {
    if (bfl_detail::bfltree_is_vertex_valid(m_vertices[0])) {
      remove_branch(0);
    }
    m_vertices[0].out_degree = 0;
    m_vertices[0].vertex() = std::forward<VProp>(vp);
    ++m_vertex_count;
    return 0;
  }
};

/**
 * This is the tree-storage specifier for a D-Ary von Emde Boas Layout tree of a given Arity.
 */
template <std::size_t Arity = 2>
struct vebl_d_ary_tree_storage {};

template <typename VertexDescriptor, typename EdgeDescriptor, std::size_t Arity>
struct tree_storage<VertexDescriptor, EdgeDescriptor, vebl_d_ary_tree_storage<Arity>> {
  using type = vebl_d_ary_tree<Arity, VertexDescriptor, EdgeDescriptor>;
};

template <std::size_t Arity>
struct tree_storage_traits<vebl_d_ary_tree_storage<Arity>> {
  static constexpr bool is_rand_access_v = true;
  static constexpr bool is_bidir_v = true;
  static constexpr bool is_directed_v = true;

  using directed_category = std::conditional_t<is_bidir_v, bidirectional_tag,
                                               std::conditional_t<is_directed_v, directed_tag, undirected_tag>>;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = std::size_t;
  using vertex_descriptor = std::size_t;
  using edges_size_type = std::size_t;
  using edge_descriptor = bfl_detail::bfltree_edge_desc;
};

#define BAGL_VEBL_D_ARY_TREE_ARGS \
  std::size_t Arity, typename VertexProperties, typename EdgeProperties, typename GraphProperties
#define BAGL_VEBL_D_ARY_TREE vebl_d_ary_tree<Arity, VertexProperties, EdgeProperties, GraphProperties>

/**
 * Standard swap function. Swaps the contents of two objects.
 * \param lhs The left-hand-side of the swap.
 * \param rhs The right-hand-side of the swap.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
void swap(BAGL_VEBL_D_ARY_TREE& lhs, BAGL_VEBL_D_ARY_TREE& rhs) {
  lhs.swap(rhs);
}

/***********************************************************************************************
 *                             IncidenceGraphConcept
 * ********************************************************************************************/

/**
 * Returns the source vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The source vertex of the given edge descriptor.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor source(typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e,
                                                        const BAGL_VEBL_D_ARY_TREE& /*unused*/) {
  return (e.target_vertex - 1) / Arity;
}

/**
 * Returns the target vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The target vertex of the given edge descriptor.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor target(typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e,
                                                        const BAGL_VEBL_D_ARY_TREE& /*unused*/) {
  return e.target_vertex;
}

/**
 * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the out-edges of a given vertex descriptor.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto out_edges(const typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BAGL_VEBL_D_ARY_TREE& g) {
  // Arity * v + 1 + edge_index;
  return BAGL_VEBL_D_ARY_TREE::view_vertex_range_as_edges(
      g.make_valid_vertex_range(Arity * v + 1, Arity * (v + 1) + 1));
}

/**
 * Returns the out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The out-degree of the given vertex descriptor.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
std::size_t out_degree(const typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BAGL_VEBL_D_ARY_TREE& g) {
  return g.get_vertex_value_impl(v).out_degree;
}

/***********************************************************************************************
 *                             BidirectionalGraphConcept
 * ********************************************************************************************/

/**
 * Returns the edge iterator range for the in-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the in-edges of a given vertex descriptor.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto in_edges(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BAGL_VEBL_D_ARY_TREE& /*g*/) {
  return BAGL_VEBL_D_ARY_TREE::view_vertex_range_as_edges(
      BAGL_VEBL_D_ARY_TREE::make_unsafe_vertex_range(v, (v == 0 ? v : v + 1)));
}

/**
 * Returns the in-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
std::size_t in_degree(const typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BAGL_VEBL_D_ARY_TREE& /*g*/) {
  if (v == 0) {
    return 0;
  }
  return 1;
}

/**
 * Returns the in-degree plus out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree plus out-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
std::size_t degree(const typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BAGL_VEBL_D_ARY_TREE& g) {
  return in_degree(v, g) + out_degree(v, g);
}

/***********************************************************************************************
 *                             VertexListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for all the vertices of the tree.
 * \param g The graph.
 * \return The vertex iterator range for all the vertices of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto vertices(const BAGL_VEBL_D_ARY_TREE& g) {
  return g.make_valid_vertex_range(0, g.m_vertices.size());
}

/**
 * Returns the size of the tree (the number of vertices it contains).
 * \param g The graph.
 * \return The size of the tree (the number of vertices it contains).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertices_size_type num_vertices(const BAGL_VEBL_D_ARY_TREE& g) {
  return g.m_vertex_count;
}

/***********************************************************************************************
 *                             EdgeListGraphConcept
 * ********************************************************************************************/

/**
 * Returns the edge iterator range for all the edges of the tree.
 * \param g The graph.
 * \return The edge iterator range for all the edges of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto edges(const BAGL_VEBL_D_ARY_TREE& g) {
  return BAGL_VEBL_D_ARY_TREE::view_vertex_range_as_edges(g.make_valid_vertex_range(1, g.m_vertices.size()));
}

/**
 * Returns the number of edges in the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertices_size_type num_edges(const BAGL_VEBL_D_ARY_TREE& g) {
  return num_vertices(g) - 1;
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
template <BAGL_VEBL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_VEBL_D_ARY_TREE::edge_descriptor, bool> edge(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor u,
                                                                     typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                                                                     const BAGL_VEBL_D_ARY_TREE& /*unused*/) {
  using Edge = typename BAGL_VEBL_D_ARY_TREE::edge_descriptor;
  if (v == 0 || u != (v - 1) / Arity) {
    return {BAGL_VEBL_D_ARY_TREE::null_edge(), false};
  }
  return {Edge(v), true};
}

/***********************************************************************************************
 *                             AdjacencyGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto adjacent_vertices(const typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BAGL_VEBL_D_ARY_TREE& g) {
  return g.make_valid_vertex_range(Arity * v + 1, Arity * (v + 1) + 1);
}

/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

/**
 * Returns the vertex-descriptor of the root of the tree.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor tree_root(const BAGL_VEBL_D_ARY_TREE& /*g*/) {
  return 0;
}

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto children(const typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor& v, const BAGL_VEBL_D_ARY_TREE& g) {
  return adjacent_vertices(v, g);
}

/***********************************************************************************************
 *                             BidirectionalTreeConcept
 * ********************************************************************************************/

/**
 * Returns the parent vertex of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The parent vertex of the given vertex descriptor (will be null_vertex() if it is the root (no parent)).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor parent(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                                                        const BAGL_VEBL_D_ARY_TREE& /*unused*/) {
  if (v == 0) {
    return BAGL_VEBL_D_ARY_TREE::null_vertex();
  }
  return (v - 1) / Arity;
}

/***********************************************************************************************
 *                             InvAdjacencyGraphConcept
 * ********************************************************************************************/

/**
 * Returns the vertex iterator range for the parent-vertex of a given vertex of the tree.
 * \param v The vertex descriptor whose parent is sought.
 * \param g The graph.
 * \return The vertex iterator range for the parent-vertex of a given vertex of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto inv_adjacent_vertices(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BAGL_VEBL_D_ARY_TREE& /*g*/) {
  if (v == 0) {
    return BAGL_VEBL_D_ARY_TREE::make_unsafe_vertex_range(0, 0);
  }
  return BAGL_VEBL_D_ARY_TREE::make_unsafe_vertex_range((v - 1) / Arity, (v - 1) / Arity + 1);
}

/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

/**
 * Removes a branch (sub-tree) starting from and including the given vertex.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param g The graph.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
void remove_branch(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g) {
  return g.remove_branch(v);
}

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex.
 * \param v The root of the sub-tree to be removed.
 * \param g The graph.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
void clear_children(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g) {
  return g.clear_children(v);
}

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor create_root(BAGL_VEBL_D_ARY_TREE& g) {
  if (bfl_detail::bfltree_is_vertex_valid(g.m_vertices[0])) {
    remove_branch(0, g);
  }
  g.m_vertices[0].out_degree = 0;
  ++g.m_vertex_count;
  return 0;
}

/**
 * Adds a child vertex to the given parent vertex, and default-initializes the properties of
 * the newly created vertex and edge.
 * \param v The parent vertex to which a child will be added.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto add_child(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g) {
  using VProp = typename BAGL_VEBL_D_ARY_TREE::vertex_property_type;
  using EProp = typename BAGL_VEBL_D_ARY_TREE::edge_property_type;
  return g.add_child(v, VProp(), EProp());
}

/***********************************************************************************************
 *                             MutablePropertyTreeConcept
 * ********************************************************************************************/

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param vp The vertex-property to assign to the newly created root vertex.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename VProp>
typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor create_root(BAGL_VEBL_D_ARY_TREE& g, VProp&& vp) {
  if (bfl_detail::bfltree_is_vertex_valid(g.m_vertices[0])) {
    remove_branch(0, g);
  }
  g.m_vertices[0].out_degree = 0;
  using VertexProp = typename BAGL_VEBL_D_ARY_TREE::vertex_property_type;
  g.m_vertices[0].vertex() = VertexProp{std::forward<VProp>(vp)};
  ++g.m_vertex_count;
  return 0;
}

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
 * vertex to the given property value.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename VProp>
auto add_child(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g, VProp&& vp) {
  using EProp = typename BAGL_VEBL_D_ARY_TREE::edge_property_type;
  return g.add_child(v, std::forward<VProp>(vp), EProp());
}

/**
 * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
 * vertex and edge to the given property values.
 * \param v The parent vertex to which a child will be added.
 * \param vp The property value to be moved into the newly created vertex.
 * \param ep The property value to be moved into the newly created edge.
 * \param g The graph.
 * \return A pair consisting of the newly created vertex and edge (descriptors).
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename VProp, typename EProp>
auto add_child(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g, VProp&& vp, EProp&& ep) {
  return g.add_child(v, std::forward<VProp>(vp), std::forward<EProp>(ep));
}

/**
 * Removes a branch (sub-tree) starting from and including the given vertex, while
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 * \note The first vertex-property to figure in the output range is that of the vertex v.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename OutputIter>
OutputIter remove_branch(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g,
                         OutputIter it_out) {
  return g.remove_branch(v, it_out);
}

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
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                                                BAGL_VEBL_D_ARY_TREE& g, VertexOIter vit_out, EdgeOIter eit_out) {
  return g.remove_branch(v, vit_out, eit_out);
}

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param it_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename OutputIter>
OutputIter clear_children(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g,
                          OutputIter it_out) {
  return g.clear_children(v, it_out);
}

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex, while
 * recording the vertex-properties of all the removed vertices into an output-iterator.
 * \param v The root of the sub-tree to be removed.
 * \param vit_out An output iterator (with vertex-properties as value-type) that can store the removed vertices.
 * \param eit_out An output iterator (with edge-properties as value-type) that can store the removed edges.
 * \param g The graph.
 * \return The output-iterator after the collection of all the removed vertices.
 */
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v,
                                                 BAGL_VEBL_D_ARY_TREE& g, VertexOIter vit_out, EdgeOIter eit_out) {
  return g.clear_children(v, vit_out, eit_out);
}

/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
const auto& get(const BAGL_VEBL_D_ARY_TREE& g, typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v) {
  return g[v];
}
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto& get(BAGL_VEBL_D_ARY_TREE& g, typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v) {
  return g[v];
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
const auto& get(const BAGL_VEBL_D_ARY_TREE& g, typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e) {
  return g[e];
}
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto& get(BAGL_VEBL_D_ARY_TREE& g, typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
const auto& get(const BAGL_VEBL_D_ARY_TREE& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto& get(BAGL_VEBL_D_ARY_TREE& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename VProp>
void put(BAGL_VEBL_D_ARY_TREE& g, typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-property associated to the given edge descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename EProp>
void put(BAGL_VEBL_D_ARY_TREE& g, typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e, EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Sets the edge-property associated to the given edge descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS, typename GProp>
void put(BAGL_VEBL_D_ARY_TREE& g, graph_bundle_t /*unused*/, GProp&& value) {
  g[graph_bundle] = std::forward<GProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto& get_property(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, BAGL_VEBL_D_ARY_TREE& g) {
  return g.get_property(v);
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
const auto& get_property(typename BAGL_VEBL_D_ARY_TREE::vertex_descriptor v, const BAGL_VEBL_D_ARY_TREE& g) {
  return g.get_property(v);
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto& get_property(typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e, BAGL_VEBL_D_ARY_TREE& g) {
  return g.get_property(e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
const auto& get_property(typename BAGL_VEBL_D_ARY_TREE::edge_descriptor e, const BAGL_VEBL_D_ARY_TREE& g) {
  return g.get_property(e);
}

// Returns a reference to the graph-property associated to the graph.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
auto& get_property(BAGL_VEBL_D_ARY_TREE& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Returns a const-reference to the graph-property associated to the graph.
template <BAGL_VEBL_D_ARY_TREE_ARGS>
const auto& get_property(const BAGL_VEBL_D_ARY_TREE& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Handle graph property tags, also handles graph_bundle_t.

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>, lookup_one_property_t<GraphProperties, Tag>&>
get_property(BAGL_VEBL_D_ARY_TREE& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>,
                 const lookup_one_property_t<GraphProperties, Tag>&>
get_property(const BAGL_VEBL_D_ARY_TREE& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle>
struct property_map<BAGL_VEBL_D_ARY_TREE, T Bundle::*> {
  using non_const_bundle = std::remove_cv_t<Bundle>;
  using non_const_t = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v =
      std::is_convertible_v<typename BAGL_VEBL_D_ARY_TREE::vertex_bundled*, non_const_bundle*>;
  static constexpr bool is_edge_bundle_v =
      std::is_convertible_v<typename BAGL_VEBL_D_ARY_TREE::edge_bundled*, non_const_bundle*>;
  using tag_type = std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                                      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<non_const_t, BAGL_VEBL_D_ARY_TREE, tag_type>;
  using const_type = bundle_member_property_map<const non_const_t, const BAGL_VEBL_D_ARY_TREE, tag_type>;
};

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, BAGL_VEBL_D_ARY_TREE& g) {
  return typename property_map<BAGL_VEBL_D_ARY_TREE, T Bundle::*>::type(&g, p);
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, const BAGL_VEBL_D_ARY_TREE& g) {
  return typename property_map<BAGL_VEBL_D_ARY_TREE, T Bundle::*>::const_type(&g, p);
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p, const BAGL_VEBL_D_ARY_TREE& g, const Key& k) {
  return (g[k]).*p;
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_VEBL_D_ARY_TREE& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct vebl_d_ary_tree_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;
    using type = tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type = tagged_in_property_property_map<value_type, const Graph, Tag>;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<vebl_d_ary_tree_tag> {
  using type = vebl_d_ary_tree_property_selector;
};

template <>
struct edge_property_selector<vebl_d_ary_tree_tag> {
  using type = vebl_d_ary_tree_property_selector;
};

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Property>
auto get(Property p, BAGL_VEBL_D_ARY_TREE& g) {
  using Map = typename property_map<BAGL_VEBL_D_ARY_TREE, Property>::type;
  return Map(&g, p);
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Property>
auto get(Property p, const BAGL_VEBL_D_ARY_TREE& g) {
  using Map = typename property_map<BAGL_VEBL_D_ARY_TREE, Property>::const_type;
  return Map(&g, p);
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Property, typename Key>
decltype(auto) get(Property p, const BAGL_VEBL_D_ARY_TREE& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_VEBL_D_ARY_TREE& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_VEBL_D_ARY_TREE_ARGS, typename Property, typename Key, typename Value>
void put(Property p, BAGL_VEBL_D_ARY_TREE& g, const Key& k, Value&& val) {
  get_property_value(g.get_property(k), p) = std::forward<Value>(val);
}

#undef BAGL_VEBL_D_ARY_TREE_ARGS
#undef BAGL_VEBL_D_ARY_TREE

}  // namespace bagl

#endif  // BAGL_BAGL_VEBL_D_ARY_TREE_H_
