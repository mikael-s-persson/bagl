// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_BFL_D_ARY_TREE_H_
#define BAGL_BAGL_BFL_D_ARY_TREE_H_

#include <ranges>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/tree_concepts.h"
#include "bagl/tree_traits.h"
#include "bagl/detail/bfl_tree_ranges.h"

namespace bagl {

// This library provides a class that implements a Breadth-first Layout D-Ary tree that is tailored 
// to store elements of a tree as if their were inserted breadth-first in a contiguous array. This type
// of tree structure is good for both breadth-first search and depth-first search because 
// of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
// should be kept balanced, and this implementation assumes that. The storage pattern is 
// similar to a binary heap tree-structure.

struct bfl_d_ary_tree_tag {};

/**
 * This class implements a D-Ary Breadth-first Layout tree that is tailored 
 * to store elements of a tree as if their were inserted breadth-first in a contiguous array. This type
 * of tree structure is good for both breadth-first search and depth-first search because 
 * of locality of reference issues. Ideally, for the least amount of wasted memory, the tree 
 * should be kept balanced, and this implementation assumes that. The storage pattern is 
 * similar to a binary heap tree-structure.
 * \tparam Arity The arity of the tree, e.g., 2 means a binary-tree.
 * \tparam VertexProperties A type to be attached to each vertex in the tree.
 * \tparam EdgeProperties A type to be attached to each edge in the tree.
 */
template <std::size_t Arity = 2, typename VertexProperties = no_property,
          typename EdgeProperties = no_property>
class bfl_d_ary_tree {
 public:
  using self = bfl_d_ary_tree<Arity, VertexProperties, EdgeProperties>;

  using vertex_property_type = VertexProperties;
  using edge_property_type = EdgeProperties;

  using vertex_bundled = vertex_property_type;
  using edge_bundled = edge_property_type;
  using graph_bundled = void;

  using value_type = bfl_detail::bfltree_value_type<vertex_property_type,
                                                       edge_property_type>;

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
  static vertex_descriptor null_vertex() {
    return std::numeric_limits<std::size_t>::max();
  }

  using raw_vertex_range = std::ranges::iota_view<std::size_t, std::size_t>;

  using edge_validity = bfl_detail::bfltree_edge_validity<container_type>;
  using in_edge_range = std::ranges::transform_view<raw_vertex_range, bfl_detail::bfltree_edge_desc_from_target>;
  using edge_range = std::ranges::filter_view<in_edge_range, edge_validity>;
  using out_edge_range = edge_range;

  using vertex_validity = bfl_detail::bfltree_vertex_validity<container_type>;
  using vertex_range = std::ranges::filter_view<raw_vertex_range, vertex_validity>;
  using adjacency_range = vertex_range;
  using child_vertex_range = vertex_range;
  using inv_adjacency_range = raw_vertex_range;

  using directed_category = directed_tag;
  using edge_parallel_category = disallow_parallel_edge_tag;

  struct traversal_category : virtual public incidence_graph_tag,
                              virtual public adjacency_graph_tag,
                              virtual public bidirectional_graph_tag,
                              virtual public vertex_list_graph_tag,
                              virtual public edge_list_graph_tag {};

  using graph_tag = bfl_d_ary_tree_tag;

  // private:
  container_type m_vertices;
  vertices_size_type m_vertex_count{0};

  /**
     * Construct the D-ary BF-tree with a given reserved depth.
     * \param aDepth The depth of the graph to reserve space for.
     */
  explicit bfl_d_ary_tree(vertices_size_type aDepth = 0) {
    vertices_size_type vert_count = 1;
    vertices_size_type accum = 1;
    for (vertices_size_type i = 0; i < aDepth; ++i) {
      accum *= Arity;
      vert_count += accum;
    }
    m_vertices.resize(vert_count);
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
  [[nodiscard]] std::size_t capacity() const { return m_vertices.size(); }

  /**
     * Returns the depth of the tree.
     * \return The depth of the tree.
     */
  [[nodiscard]] std::size_t depth() const {
    vertices_size_type vert_count = 1;
    vertices_size_type accum = 1;
    vertices_size_type depth_count = 0;
    for (; vert_count < m_vertices.size(); ++depth_count) {
      accum *= Arity;
      vert_count += accum;
    }
    return depth_count;
  }

  /**
     * Standard swap function.
     */
  void swap(self& rhs) {
    using std::swap;
    m_vertices.swap(rhs.m_vertices);
    swap(m_vertex_count, rhs.m_vertex_count);
  }

  /**
     * Clears the tree of all vertices and edges.
     */
  void clear() {
    m_vertices.clear();
    m_vertices.resize(1);
    m_vertex_count = 0;
  }

  /**
     * Indexing operator. Returns a reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by reference, associated to the given vertex descriptor.
     */
  vertex_property_type& operator[](vertex_descriptor v) {
    return m_vertices[v].vertex();
  }
  /**
     * Indexing operator. Returns a const-reference to the vertex-property associated to the given vertex descriptor.
     * \param v_i The vertex descriptor of the sought-after vertex-property.
     * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
     */
  const vertex_property_type& operator[](vertex_descriptor v) const {
    return m_vertices[v].vertex();
  }
  /**
     * Indexing operator. Returns a reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by reference, associated to the given edge descriptor.
     */
  edge_property_type& operator[](edge_descriptor e) {
    return m_vertices[e.target_vertex].edge();
  }
  /**
     * Indexing operator. Returns a const-reference to the edge-property associated to the given edge descriptor.
     * \param e_i The edge descriptor of the sought-after edge-property.
     * \return The edge-property, by const-reference, associated to the given edge descriptor.
     */
  const edge_property_type& operator[](edge_descriptor e) const {
    return m_vertices[e.target_vertex].edge();
  }

  /**
     * Adds a child vertex to the given parent vertex, and initializes the properties of the newly created 
     * vertex and edge to the given property values.
     * \param v The parent vertex to which a child will be added.
     * \param vp The property value to be forwarded into the newly created vertex.
     * \param ep The property value to be forwarded into the newly created edge.
     * \return A pair consisting of the newly created vertex and edge (descriptors).
     */
  template <typename VP, typename EP>
  std::pair<vertex_descriptor, edge_descriptor> add_child(
      vertex_descriptor v, VP&& vp, EP&& ep) {
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return std::pair<vertex_descriptor, edge_descriptor>(
          null_vertex(), edge_descriptor(null_vertex()));
    }
    std::size_t result = Arity * v + 1;
    for (; result < Arity * (v + 1) + 1; ++result) {
      if ((result >= m_vertices.size()) ||
          !bfl_detail::bfltree_is_vertex_valid(m_vertices[result])) {
        break;
      }
    }
    if (result == Arity * (v + 1) + 1) {
      return std::pair<vertex_descriptor, edge_descriptor>(
          null_vertex(), edge_descriptor(null_vertex()));
    }
    if (result >= m_vertices.size()) {
      m_vertices.resize(result + 1);
    }
    m_vertices[result].out_degree = 0;
    m_vertices[result].vertex() = std::forward<VP>(vp);
    m_vertices[result].edge() = std::forward<EP>(ep);
    ++(m_vertices[v].out_degree);
    ++m_vertex_count;
    return std::make_pair(result, edge_descriptor(result));
  }

  template <typename VertexOIter, typename EdgeOIter>
  void clear_children_impl(vertex_descriptor v, VertexOIter& vit_out,
                           EdgeOIter& eit_out) {
    // this traversal order is intentional (traverse pre-order depth-first, and
    // delay removal of empty tail elements as much as possible, such that it is only required once).
    for (std::size_t i = 0; i < Arity; ++i) {
      vertex_descriptor next_v = Arity * v + 1 + i;
      if (next_v >= m_vertices.size()) {
        break;
      }
      if (!bfl_detail::bfltree_is_vertex_valid(m_vertices[next_v])) {
        continue;
      }
      --m_vertex_count;
      *(vit_out++) = std::move(m_vertices[next_v].vertex());
      *(eit_out++) = std::move(m_vertices[next_v].edge());
      clear_children_impl(next_v, vit_out, eit_out);
    }
    m_vertices[v].out_degree = std::numeric_limits<std::size_t>::max();
    // remove empty vertices from the end of the container:
    if (v == m_vertices.size() - 1) {
      while ((v > 0) &&
             !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
        --v;
      }
      ++v;
      m_vertices.erase(m_vertices.begin() + v, m_vertices.end());
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
  std::pair<VertexOIter, EdgeOIter> clear_children(vertex_descriptor v,
                                                   VertexOIter vit_out,
                                                   EdgeOIter eit_out) {
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return std::pair<VertexOIter, EdgeOIter>(
          vit_out, eit_out);  // vertex is already deleted.
    }
    clear_children_impl(v, vit_out, eit_out);
    return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
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
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return vit_out;  // vertex is already deleted.
    }
    container_detail::ignore_output_iter eit_out;
    clear_children_impl(v, vit_out, eit_out);
    return vit_out;
  }

  /**
     * Removes a branch (sub-tree) starting from but excluding the given vertex.
     * \param v The root of the sub-tree to be removed.
     */
  void clear_children(vertex_descriptor v) {
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return;  // vertex is already deleted.
    }
    container_detail::ignore_output_iter vit_out;
    container_detail::ignore_output_iter eit_out;
    clear_children_impl(v, vit_out, eit_out);
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
  std::pair<VertexOIter, EdgeOIter> remove_branch(vertex_descriptor v,
                                                  VertexOIter vit_out,
                                                  EdgeOIter eit_out) {
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return std::pair<VertexOIter, EdgeOIter>(
          vit_out, eit_out);  // vertex is already deleted.
    }
    if (v == 0) {
      *(eit_out++) = edge_property_type();
    } else {
      // in-edge:  u = (v - 1) / Arity;  e_id = (v - 1) % Arity;
      *(eit_out++) = std::move(m_vertices[v].edge());
      // if the node is not the root one, then update the out-degree of the parent node:
      m_vertices[(v - 1) / Arity].out_degree -= 1;
    }
    --m_vertex_count;
    *(vit_out++) = std::move(m_vertices[v].vertex());
    clear_children_impl(v, vit_out, eit_out);
    return std::pair<VertexOIter, EdgeOIter>(vit_out, eit_out);
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
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return vit_out;  // vertex is already deleted.
    }
    if (v != 0) {
      // if the node is not the root one, then update the out-degree of the parent node:
      m_vertices[(v - 1) / Arity].out_degree -= 1;
    }
    --m_vertex_count;
    *(vit_out++) = std::move(m_vertices[v].vertex());
    container_detail::ignore_output_iter eit_out;
    clear_children_impl(v, vit_out, eit_out);
    return vit_out;
  }

  /**
     * Removes a branch (sub-tree) starting from and including the given vertex.
     * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
     */
  void remove_branch(vertex_descriptor v) {
    if ((v >= m_vertices.size()) ||
        !bfl_detail::bfltree_is_vertex_valid(m_vertices[v])) {
      return;  // vertex is already deleted.
    }
    if (v != 0) {
      // if the node is not the root one, then update the out-degree of the parent node:
      m_vertices[(v - 1) / Arity].out_degree -= 1;
    }
    --m_vertex_count;
    container_detail::ignore_output_iter vit_out;
    container_detail::ignore_output_iter eit_out;
    clear_children_impl(v, vit_out, eit_out);
  }
};

/**
 * This is the tree-storage specifier for a D-Ary Breadth-first Layout tree of a given Arity.
 */
template <std::size_t Arity = 2>
struct bfl_d_ary_tree_storage {};

template <typename VertexDescriptor, typename EdgeDescriptor, std::size_t Arity>
struct tree_storage<VertexDescriptor, EdgeDescriptor,
                    bfl_d_ary_tree_storage<Arity>> {
  using type = bfl_d_ary_tree<Arity, VertexDescriptor, EdgeDescriptor>;
};

template <std::size_t Arity>
struct tree_storage_traits<bfl_d_ary_tree_storage<Arity>> {
  static constexpr bool is_rand_access_v = true;
  static constexpr bool is_bidir_v = true;
  static constexpr bool is_directed_v = true;

  using directed_category = std::conditional_t<
      is_bidir_v, bidirectional_tag,
      std::conditional_t<is_directed_v, directed_tag,
                               undirected_tag>>;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = std::size_t;
  using vertex_descriptor = std::size_t;
  using edges_size_type = std::size_t;
  using edge_descriptor = bfl_detail::bfltree_edge_desc;
};

#define BAGL_BFL_D_ARY_TREE_ARGS \
  std::size_t Arity, typename VertexProperties, typename EdgeProperties
#define BAGL_BFL_D_ARY_TREE \
  bfl_d_ary_tree<Arity, VertexProperties, EdgeProperties>

/**
 * Standard swap function. Swaps the contents of two objects.
 * \param lhs The left-hand-side of the swap.
 * \param rhs The right-hand-side of the swap.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void swap(BAGL_BFL_D_ARY_TREE& lhs, BAGL_BFL_D_ARY_TREE& rhs) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto source(
    typename BAGL_BFL_D_ARY_TREE::edge_descriptor e,
    const BAGL_BFL_D_ARY_TREE& /*unused*/) {
  return (e.target_vertex - 1) / Arity;
}

/**
 * Returns the target vertex of a given edge descriptor in the tree.
 * \param e The edge descriptor.
 * \param g The graph.
 * \return The target vertex of the given edge descriptor.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto target(
    typename BAGL_BFL_D_ARY_TREE::edge_descriptor e,
    const BAGL_BFL_D_ARY_TREE& /*unused*/) {
  return e.target_vertex;
}

/**
 * Returns the edge iterator range for the out-edges of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The edge iterator range for the out-edges of a given vertex descriptor.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::out_edge_iterator,
          typename BAGL_BFL_D_ARY_TREE::out_edge_iterator>
out_edges(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
          const BAGL_BFL_D_ARY_TREE& g) {
  using OutIter = typename bfl_d_ary_tree<Arity, VertexProperties,
                                          EdgeProperties>::out_edge_iterator;
  using RawContainer = typename bfl_d_ary_tree<Arity, VertexProperties,
                                               EdgeProperties>::container_type;
  // Arity * v + 1 + edge_index;
  bfl_detail::bfltree_eiter v_beg(
      bfl_detail::bfltree_edge_desc(Arity * v + 1));
  bfl_detail::bfltree_eiter v_end(
      bfl_detail::bfltree_edge_desc(Arity * (v + 1) + 1));
  return std::pair<OutIter, OutIter>(
      OutIter(bfl_detail::bfltree_edge_validity<RawContainer>(&g.m_vertices),
              v_beg, v_end),
      OutIter(bfl_detail::bfltree_edge_validity<RawContainer>(&g.m_vertices),
              v_end, v_end));
}

/**
 * Returns the out-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The out-degree of the given vertex descriptor.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::size_t out_degree(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                       const BAGL_BFL_D_ARY_TREE& g) {
  return g.m_vertices[v].out_degree;
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::in_edge_iterator,
          typename BAGL_BFL_D_ARY_TREE::in_edge_iterator>
in_edges(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
         const BAGL_BFL_D_ARY_TREE& g) {
  using InIter = typename bfl_d_ary_tree<Arity, VertexProperties,
                                         EdgeProperties>::in_edge_iterator;
  if (v == 0) {
    return std::make_pair(InIter(bfl_detail::bfltree_edge_desc(0)),
                          InIter(bfl_detail::bfltree_edge_desc(0)));
  }
  return std::make_pair(InIter(bfl_detail::bfltree_edge_desc(v)),
                        InIter(bfl_detail::bfltree_edge_desc(v + 1)));
}

/**
 * Returns the in-degree of a given vertex descriptor in the tree.
 * \param v The vertex descriptor.
 * \param g The graph.
 * \return The in-degree of the given vertex descriptor (will be 1 or 0 (root or invalid vertex)).
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::size_t in_degree(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                      const BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::size_t degree(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                   const BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::vertex_iterator,
          typename BAGL_BFL_D_ARY_TREE::vertex_iterator>
vertices(const BAGL_BFL_D_ARY_TREE& g) {
  using VIter = typename bfl_d_ary_tree<Arity, VertexProperties,
                                        EdgeProperties>::vertex_iterator;
  using RawContainer = typename bfl_d_ary_tree<Arity, VertexProperties,
                                               EdgeProperties>::container_type;
  bfl_detail::bfltree_viter v_beg(0);
  bfl_detail::bfltree_viter v_end(g.m_vertices.size());
  return std::pair<VIter, VIter>(
      VIter(bfl_detail::bfltree_vertex_validity<RawContainer>(&g.m_vertices),
            v_beg, v_end),
      VIter(bfl_detail::bfltree_vertex_validity<RawContainer>(&g.m_vertices),
            v_end, v_end));
}

/**
 * Returns the size of the tree (the number of vertices it contains).
 * \param g The graph.
 * \return The size of the tree (the number of vertices it contains).
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto num_vertices(const BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::edge_iterator,
          typename BAGL_BFL_D_ARY_TREE::edge_iterator>
edges(const BAGL_BFL_D_ARY_TREE& g) {
  using EIter = typename BAGL_BFL_D_ARY_TREE::edge_iterator;
  using RawContainer = typename BAGL_BFL_D_ARY_TREE::container_type;

  bfl_detail::bfltree_eiter v_beg(bfl_detail::bfltree_edge_desc(1));
  bfl_detail::bfltree_eiter v_end(
      bfl_detail::bfltree_edge_desc(g.m_vertices.size()));
  return std::pair<EIter, EIter>(
      EIter(bfl_detail::bfltree_edge_validity<RawContainer>(&g.m_vertices),
            v_beg, v_end),
      EIter(bfl_detail::bfltree_edge_validity<RawContainer>(&g.m_vertices),
            v_end, v_end));
}

/**
 * Returns the number of edges in the tree.
 * \param g The graph.
 * \return The number of edges in the tree.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto num_edges(const BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::edge_descriptor, bool> edge(
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor u,
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
    const BAGL_BFL_D_ARY_TREE& /*unused*/) {
  using Edge = typename BAGL_BFL_D_ARY_TREE::edge_descriptor;
  if (u == (v - 1) / Arity) {
    return {Edge(v), true};
  }
  return {Edge(0), false};
}

/***********************************************************************************************
 *                             TreeConcept
 * ********************************************************************************************/

/**
 * Returns the vertex-descriptor of the root of the tree.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto tree_root(const BAGL_BFL_D_ARY_TREE& g) {
  return 0;
}

/**
 * Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
 * \param v The vertex descriptor whose children are sought.
 * \param g The graph.
 * \return The vertex iterator range for all the child-vertices of a given vertex of the tree.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::vertex_iterator,
          typename BAGL_BFL_D_ARY_TREE::vertex_iterator>
children(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
               const BAGL_BFL_D_ARY_TREE& g) {
  using VIter = typename bfl_d_ary_tree<Arity, VertexProperties,
                                        EdgeProperties>::vertex_iterator;
  using RawContainer = typename bfl_d_ary_tree<Arity, VertexProperties,
                                               EdgeProperties>::container_type;
  bfl_detail::bfltree_viter v_beg(Arity * v + 1);
  bfl_detail::bfltree_viter v_end(Arity * (v + 1) + 1);
  return std::pair<VIter, VIter>(
      VIter(bfl_detail::bfltree_vertex_validity<RawContainer>(&g.m_vertices),
            v_beg, v_end),
      VIter(bfl_detail::bfltree_vertex_validity<RawContainer>(&g.m_vertices),
            v_end, v_end));
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto adjacent_vertices(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                  const BAGL_BFL_D_ARY_TREE& g) {
  return children(v, g);
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
std::pair<typename BAGL_BFL_D_ARY_TREE::inv_adjacency_iterator,
          typename BAGL_BFL_D_ARY_TREE::inv_adjacency_iterator>
inv_adjacent_vertices(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                      const BAGL_BFL_D_ARY_TREE& g) {
  using InvAdjIter =
      typename bfl_d_ary_tree<Arity, VertexProperties,
                              EdgeProperties>::inv_adjacency_iterator;
  if (v == 0) {
    return std::pair<InvAdjIter, InvAdjIter>(bfl_detail::bfltree_viter(0),
                                             bfl_detail::bfltree_viter(0));
  }
  return std::pair<InvAdjIter, InvAdjIter>(
      bfl_detail::bfltree_viter((v - 1) / Arity),
      bfl_detail::bfltree_viter((v - 1) / Arity + 1));
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto parent(
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
    const BAGL_BFL_D_ARY_TREE& /*unused*/) {
  if (v == 0) {
    return BAGL_BFL_D_ARY_TREE::null_vertex();
  }
  return (v - 1) / Arity;
}

/***********************************************************************************************
 *                             MutableTreeConcept
 * ********************************************************************************************/

/**
 * Removes a branch (sub-tree) starting from and including the given vertex.
 * \param v The vertex to remove, along with the sub-tree rooted at that vertex.
 * \param g The graph.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void remove_branch(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                   BAGL_BFL_D_ARY_TREE& g) {
  return g.remove_branch(v);
}

/**
 * Removes a branch (sub-tree) starting from but excluding the given vertex.
 * \param v The root of the sub-tree to be removed.
 * \param g The graph.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void clear_children(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                    BAGL_BFL_D_ARY_TREE& g) {
  return g.clear_children(v);
}

/**
 * Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
 * \param g The graph.
 * \return The vertex-descriptor of the root of the tree.
 */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto create_root(BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto add_child(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                 BAGL_BFL_D_ARY_TREE& g) {
  using VProp = typename BAGL_BFL_D_ARY_TREE::vertex_property_type;
  using EProp = typename BAGL_BFL_D_ARY_TREE::edge_property_type;
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename VProp>
auto create_root(VProp&& vp, BAGL_BFL_D_ARY_TREE& g) {
  if (bfl_detail::bfltree_is_vertex_valid(g.m_vertices[0])) {
    remove_branch(0, g);
  }
  g.m_vertices[0].out_degree = 0;
  g.m_vertices[0].vertex() = std::forward<VProp>(vp);
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename VProp>
auto add_child(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v, VProp&& vp,
                 BAGL_BFL_D_ARY_TREE& g) {
  using EProp = typename bfl_d_ary_tree<Arity, VertexProperties,
                                        EdgeProperties>::edge_property_type;
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename VProp, typename EProp>
auto add_child(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v, VProp&& vp,
                 EProp&& ep, BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename OutputIter>
OutputIter remove_branch(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                         OutputIter it_out, BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch(
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v, VertexOIter vit_out,
    EdgeOIter eit_out, BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename OutputIter>
OutputIter clear_children(typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
                          OutputIter it_out, BAGL_BFL_D_ARY_TREE& g) {
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
template <BAGL_BFL_D_ARY_TREE_ARGS, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children(
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v, VertexOIter vit_out,
    EdgeOIter eit_out, BAGL_BFL_D_ARY_TREE& g) {
  return g.clear_children(v, vit_out, eit_out);
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
template <BAGL_BFL_D_ARY_TREE_ARGS>
const auto& get(
    const BAGL_BFL_D_ARY_TREE& g,
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v) {
  return g[v];
}

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param g The tree from which to draw the edge.
  * \param e The edge descriptor of the sought-after edge-property.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
const auto& get(
    const BAGL_BFL_D_ARY_TREE& g,
    typename BAGL_BFL_D_ARY_TREE::edge_descriptor e) {
  return g[e];
}

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by const-reference, to be associated to the given vertex.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void put(BAGL_BFL_D_ARY_TREE& g,
         typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
         const typename BAGL_BFL_D_ARY_TREE::vertex_property_type& value) {
  g[v] = value;
}

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by const-reference, to be associated to the given edge.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void put(BAGL_BFL_D_ARY_TREE& g, typename BAGL_BFL_D_ARY_TREE::edge_descriptor e,
         const typename BAGL_BFL_D_ARY_TREE::edge_property_type& value) {
  g[e] = value;
}

/**
  * Sets the vertex-property associated to the given vertex descriptor.
  * \param g The tree from which the vertex is drawn.
  * \param v The vertex descriptor of the vertex-property to be set.
  * \param value The vertex-property, by rvalue-reference, to be associated to the given vertex.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void put(BAGL_BFL_D_ARY_TREE& g,
         typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
         typename BAGL_BFL_D_ARY_TREE::vertex_property_type&& value) {
  g[v] = std::move(value);
}

/**
  * Sets the edge-property associated to the given edge descriptor.
  * \param g The tree from which the edge is drawn.
  * \param e The edge descriptor of the edge-property to be set.
  * \param value The edge-property, by rvalue-reference, to be associated to the given edge.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
void put(BAGL_BFL_D_ARY_TREE& g, typename BAGL_BFL_D_ARY_TREE::edge_descriptor e,
         typename BAGL_BFL_D_ARY_TREE::edge_property_type&& value) {
  g[e] = std::move(value);
}

/**
  * Returns a reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by reference, associated to the given vertex descriptor.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto& get_property(
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v, BAGL_BFL_D_ARY_TREE& g) {
  return g[v];
}

/**
  * Returns a const-reference to the vertex-property associated to the given vertex descriptor.
  * \param v The vertex descriptor of the sought-after vertex-property.
  * \param g The tree from which to draw the vertex.
  * \return The vertex-property, by const-reference, associated to the given vertex descriptor.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
const auto& get_property(
    typename BAGL_BFL_D_ARY_TREE::vertex_descriptor v,
    const BAGL_BFL_D_ARY_TREE& g) {
  return g[v];
}

/**
  * Returns a reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by reference, associated to the given edge descriptor.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
auto& get_property(typename BAGL_BFL_D_ARY_TREE::edge_descriptor e, BAGL_BFL_D_ARY_TREE& g) {
  return g[e];
}

/**
  * Returns a const-reference to the edge-property associated to the given edge descriptor.
  * \param e The edge descriptor of the sought-after edge-property.
  * \param g The tree from which to draw the edge.
  * \return The edge-property, by const-reference, associated to the given edge descriptor.
  */
template <BAGL_BFL_D_ARY_TREE_ARGS>
const auto& get_property(
    typename BAGL_BFL_D_ARY_TREE::edge_descriptor e,
    const BAGL_BFL_D_ARY_TREE& g) {
  return g[e];
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename T, typename Bundle>
struct property_map<BAGL_BFL_D_ARY_TREE, T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v =
      std::is_convertible_v<typename bfl_d_ary_tree<Arity, VertexProperties,
                                             EdgeProperties>::vertex_bundled*,
                     non_const_Bundle*>;
  using type = bundle_member_property_map<
      non_const_T, bfl_d_ary_tree<Arity, VertexProperties, EdgeProperties>,
      std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                        edge_bundle_t>>;
  using const_type = bundle_member_property_map<
      const non_const_T,
      const bfl_d_ary_tree<Arity, VertexProperties, EdgeProperties>,
      std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                        edge_bundle_t>>;
};

template <BAGL_BFL_D_ARY_TREE_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, BAGL_BFL_D_ARY_TREE& g) {
  return typename property_map<BAGL_BFL_D_ARY_TREE, T Bundle::*>::type(&g, p);
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename T, typename Bundle>
auto get(T Bundle::*p, const BAGL_BFL_D_ARY_TREE& g) {
  return typename property_map<BAGL_BFL_D_ARY_TREE, T Bundle::*>::const_type(&g, p);
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p,
                                          const BAGL_BFL_D_ARY_TREE& g,
                                          const Key& k) {
  return (g[k]).*p;
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_BFL_D_ARY_TREE& g, const Key& k, const T& val) {
  (g[k]).*p = val;
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_BFL_D_ARY_TREE& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct bfl_d_ary_tree_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;

    using type = tagged_from_bundle_property_map<value_type, Graph, Tag>;
    using const_type =
        tagged_from_bundle_property_map<const value_type, const Graph, Tag>;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<bfl_d_ary_tree_tag> {
  using type = bfl_d_ary_tree_property_selector;
};

template <>
struct edge_property_selector<bfl_d_ary_tree_tag> {
  using type = bfl_d_ary_tree_property_selector;
};

template <BAGL_BFL_D_ARY_TREE_ARGS, typename Property>
auto get(Property p, BAGL_BFL_D_ARY_TREE& g) {
  using Map = typename property_map<BAGL_BFL_D_ARY_TREE, Property>::type;
  return Map(&g, p);
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename Property>
auto get(Property p, const BAGL_BFL_D_ARY_TREE& g) {
  using Map = typename property_map<BAGL_BFL_D_ARY_TREE, Property>::const_type;
  return Map(&g, p);
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename Property, typename Key>
typename property_map_value<BAGL_BFL_D_ARY_TREE, Property>::type get(
    Property p, const BAGL_BFL_D_ARY_TREE& g, const Key& k) {
  return get_property_value(g[k], p);
}

template <BAGL_BFL_D_ARY_TREE_ARGS, typename Property, typename Key,
          typename Value>
void put(Property p, BAGL_BFL_D_ARY_TREE& g, const Key& k, Value&& val) {
  get_property_value(g[k], p) = std::forward<Value>(val);
}

#undef BAGL_BFL_D_ARY_TREE_ARGS
#undef BAGL_BFL_D_ARY_TREE

}  // namespace bagl

#endif // BAGL_BAGL_BFL_D_ARY_TREE_H_
