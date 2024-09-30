// Copyright 2005-2009 The Trustees of Indiana University.
// Authors: Jeremiah Willcock, Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

// Compressed sparse row graph type

#ifndef BAGL_BAGL_COMPRESSED_SPARSE_ROW_GRAPH_H_
#define BAGL_BAGL_COMPRESSED_SPARSE_ROW_GRAPH_H_

#include <algorithm>
#include <cassert>
#include <climits>
#include <iterator>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/detail/compressed_sparse_row_struct.h"
#include "bagl/detail/indexed_properties.h"
#include "bagl/detail/is_distributed_selector.h"
#include "bagl/filtered_graph.h"  // For keep_all
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/transform_value_property_map.h"
#include "bagl/zip_range.h"

namespace bagl {

// A tag type indicating that the graph in question is a compressed
// sparse row graph. This is an internal detail of the BGL.
struct csr_graph_tag;

// A type (edges_are_sorted_t) and a value (edges_are_sorted) used to indicate
// that the edge list passed into the CSR graph is already sorted by source
// vertex.
enum edges_are_sorted_t { edges_are_sorted };

// A type (edges_are_sorted_global_t) and a value (edges_are_sorted_global)
// used to indicate that the edge list passed into the CSR graph is already
// sorted by source vertex.
enum edges_are_sorted_global_t { edges_are_sorted_global };

// A type (edges_are_unsorted_t) and a value (edges_are_unsorted) used to
// indicate that the edge list passed into the CSR graph is not sorted by
// source vertex.  This version caches the edge information in memory, and thus
// requires only a single pass over the input data.
enum edges_are_unsorted_t { edges_are_unsorted };

// A type (edges_are_unsorted_multi_pass_t) and a value
// (edges_are_unsorted_multi_pass) used to indicate that the edge list passed
// into the CSR graph is not sorted by source vertex.  This version uses less
// memory but requires multi-pass capability on the iterators.
enum edges_are_unsorted_multi_pass_t { edges_are_unsorted_multi_pass };

// A type (edges_are_unsorted_multi_pass_global_t) and a value
// (edges_are_unsorted_multi_pass_global) used to indicate that the edge list
// passed into the CSR graph is not sorted by source vertex.  This version uses
// less memory but requires multi-pass capability on the iterators.  The
// global mapping and filtering is done here because it is often faster and it
// greatly simplifies handling of edge properties.
enum edges_are_unsorted_multi_pass_global_t { edges_are_unsorted_multi_pass_global };

// A type (construct_inplace_from_sources_and_targets_t) and a value
// (construct_inplace_from_sources_and_targets) used to indicate that mutable
// vectors of sources and targets (and possibly edge properties) are being used
// to construct the CSR graph.  These vectors are sorted in-place and then the
// targets and properties are swapped into the graph data structure.
enum construct_inplace_from_sources_and_targets_t { construct_inplace_from_sources_and_targets };

// A type (construct_inplace_from_sources_and_targets_global_t) and a value
// (construct_inplace_from_sources_and_targets_global) used to indicate that
// mutable vectors of sources and targets (and possibly edge properties) are
// being used to construct the CSR graph.  These vectors are sorted in-place
// and then the targets and properties are swapped into the graph data
// structure.  It is assumed that global indices (for distributed CSR) are
// used, and a map is required to convert those to local indices.  This
// constructor is intended for internal use by the various CSR graphs
// (sequential and distributed).
enum construct_inplace_from_sources_and_targets_global_t { construct_inplace_from_sources_and_targets_global };

// A type (edges_are_unsorted_global_t) and a value (edges_are_unsorted_global)
// used to indicate that the edge list passed into the CSR graph is not sorted
// by source vertex.  The data is also stored using global vertex indices, and
// must be filtered to choose only local vertices.  This constructor caches the
// edge information in memory, and thus requires only a single pass over the
// input data.  This constructor is intended for internal use by the
// distributed CSR constructors.
enum edges_are_unsorted_global_t { edges_are_unsorted_global };

/****************************************************************************
 * Local helper macros to reduce typing and clutter later on.               *
 ****************************************************************************/
#define BAGL_CSR_GRAPH_TEMPLATE_PARMS                                                                         \
  typename Directed, typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename Vertex, \
      typename EdgeIndex
#define BAGL_CSR_GRAPH_TYPE \
  compressed_sparse_row_graph<Directed, VertexProperty, EdgeProperty, GraphProperty, Vertex, EdgeIndex>
#define BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS \
  typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename Vertex, typename EdgeIndex
#define BAGL_DIR_CSR_GRAPH_TYPE \
  compressed_sparse_row_graph<directed_s, VertexProperty, EdgeProperty, GraphProperty, Vertex, EdgeIndex>
#define BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS \
  typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename Vertex, typename EdgeIndex
#define BAGL_BIDIR_CSR_GRAPH_TYPE \
  compressed_sparse_row_graph<bidirectional_s, VertexProperty, EdgeProperty, GraphProperty, Vertex, EdgeIndex>

/** Compressed sparse row graph.
 *
 * Vertex and EdgeIndex should be unsigned integral types and should
 * specialize numeric_limits.
 */
template <typename Directed = directed_s, typename VertexProperty = no_property, typename EdgeProperty = no_property,
          typename GraphProperty = no_property, typename Vertex = std::size_t,
          typename EdgeIndex = Vertex>
class compressed_sparse_row_graph;  // Not defined

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename Vertex, typename EdgeIndex>
class compressed_sparse_row_graph<directed_s, VertexProperty, EdgeProperty, GraphProperty, Vertex, EdgeIndex>
    : public indexed_detail::indexed_vertex_properties<BAGL_DIR_CSR_GRAPH_TYPE, VertexProperty, Vertex,
                                                       typed_identity_property_map<Vertex>> {
 public:
  using inherited_vertex_properties =
      indexed_detail::indexed_vertex_properties<compressed_sparse_row_graph, VertexProperty, Vertex,
                                                typed_identity_property_map<Vertex>>;

  // Some tests to prevent use of "void" is a property type (as was done in
  // some test cases):
  static_assert(!std::is_same_v<VertexProperty, void>);
  static_assert(!std::is_same_v<EdgeProperty, void>);
  static_assert(!std::is_same_v<GraphProperty, void>);

  // For Property Graph
  using graph_property_type = GraphProperty;
  using graph_bundled = lookup_one_property_t<GraphProperty, graph_bundle_t>;
  using forward_type = csr_detail::compressed_sparse_row_structure<EdgeProperty, Vertex, EdgeIndex>;

  /* At this time, the compressed sparse row graph can only be used to
   * create directed and bidirectional graphs. In the future,
   * undirected CSR graphs will also be supported.
   */
  // static_assert(std::is_same_v<Directed, directed_s>);

  // Concept requirements:
  // For Graph
  using vertex_descriptor = Vertex;
  using edge_descriptor = csr_detail::csr_edge_descriptor<Vertex, EdgeIndex>;
  using directed_category = directed_tag;
  using edge_parallel_category = allow_parallel_edge_tag;

  struct traversal_category : incidence_graph_tag, adjacency_graph_tag, vertex_list_graph_tag, edge_list_graph_tag {};

  static vertex_descriptor null_vertex() { return std::numeric_limits<vertex_descriptor>::max(); }

  // For VertexListGraph
  using vertices_size_type = Vertex;

  // For EdgeListGraph
  using edges_size_type = EdgeIndex;

  // For IncidenceGraph
  using degree_size_type = EdgeIndex;

  // For internal use
  using graph_tag = csr_graph_tag;

  using edge_bundled = typename forward_type::inherited_edge_properties::edge_bundled;
  using edge_push_back_type = typename forward_type::inherited_edge_properties::edge_push_back_type;
  using edge_property_type = typename forward_type::inherited_edge_properties::edge_property_type;

  // Constructors

  // Default constructor: an empty graph.
  compressed_sparse_row_graph() : property_() {}

  //  With numverts vertices
  explicit compressed_sparse_row_graph(std::size_t numverts)
      : inherited_vertex_properties(numverts), forward_(numverts) {}

  //  From number of vertices and unsorted list of edges
  template <std::ranges::forward_range EdgeRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numverts, typed_identity_property_map<std::size_t>(),
                                              keep_all());
  }

  //  From number of vertices and unsorted list of edges, plus edge properties
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numverts, typed_identity_property_map<std::size_t>(),
                                              keep_all());
  }

  //  From number of vertices and unsorted list of edges, with filter and
  //  global-to-local map
  template <std::ranges::forward_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numlocalverts, global_to_local, source_pred);
  }

  //  From number of vertices and unsorted list of edges, plus edge
  //  properties, with filter and global-to-local map
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange, typename GlobalToLocal,
            typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numlocalverts, global_to_local, source_pred);
  }

  //  From number of vertices and sorted list of edges (new interface)
  template <std::ranges::input_range EdgeRange>
  compressed_sparse_row_graph(edges_are_sorted_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              std::size_t numedges = 0, const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, typed_identity_property_map<std::size_t>(), keep_all(), numverts,
                                      numedges);
    inherited_vertex_properties::resize(numverts);
  }

  //  From number of vertices and sorted list of edges (new interface)
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_sorted_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, std::size_t numedges = 0,
                              const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, ep_rg, typed_identity_property_map<std::size_t>(), keep_all(), numverts,
                                      numedges);
    inherited_vertex_properties::resize(numverts);
  }

  //  From number of vertices and sorted list of edges, filtered and global
  //  (new interface)
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_sorted_global_t /*unused*/, EdgeRange edge_rg,
                              const GlobalToLocal& global_to_local, const SourcePred& source_pred, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, global_to_local, source_pred, numverts, 0);
    inherited_vertex_properties::resize(numverts);
  }

  //  From number of vertices and sorted list of edges (new interface)
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange, typename GlobalToLocal,
            typename SourcePred>
  compressed_sparse_row_graph(edges_are_sorted_global_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              const GlobalToLocal& global_to_local, const SourcePred& source_pred, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, ep_rg, global_to_local, source_pred, numverts, 0);
    inherited_vertex_properties::resize(numverts);
  }

  //  From number of vertices and mutable vectors of sources and targets;
  //  vectors are returned with unspecified contents but are guaranteed not to
  //  share storage with the constructed graph.
  compressed_sparse_row_graph(construct_inplace_from_sources_and_targets_t /*unused*/,
                              std::vector<vertex_descriptor>& sources, std::vector<vertex_descriptor>& targets,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), property_(prop) {
    forward_.assign_sources_and_targets_global(sources, targets, numverts, typed_identity_property_map<std::size_t>());
  }

  //  From number of vertices and mutable vectors of sources and targets,
  //  expressed with global vertex indices; vectors are returned with
  //  unspecified contents but are guaranteed not to share storage with the
  //  constructed graph.  This constructor should only be used by the
  //  distributed CSR graph.
  template <typename GlobalToLocal>
  compressed_sparse_row_graph(construct_inplace_from_sources_and_targets_global_t /*unused*/,
                              std::vector<vertex_descriptor>& sources, std::vector<vertex_descriptor>& targets,
                              std::size_t numlocalverts, GlobalToLocal global_to_local,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), property_(prop) {
    forward_.assign_sources_and_targets_global(sources, targets, numlocalverts, global_to_local);
  }

  //  From number of vertices and mutable vectors of sources, targets, and
  //  edge properties; vectors are returned with unspecified contents but are
  //  guaranteed not to share storage with the constructed graph.
  compressed_sparse_row_graph(construct_inplace_from_sources_and_targets_t /*unused*/,
                              std::vector<vertex_descriptor>& sources, std::vector<vertex_descriptor>& targets,
                              std::vector<typename forward_type::inherited_edge_properties::edge_bundled>& edge_props,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), property_(prop) {
    forward_.assign_sources_and_targets_global(sources, targets, edge_props, numverts,
                                               typed_identity_property_map<std::size_t>());
  }

  //  From number of vertices and mutable vectors of sources and targets and
  //  edge properties, expressed with global vertex indices; vectors are
  //  returned with unspecified contents but are guaranteed not to share
  //  storage with the constructed graph.  This constructor should only be
  //  used by the distributed CSR graph.
  template <typename GlobalToLocal>
  compressed_sparse_row_graph(construct_inplace_from_sources_and_targets_global_t /*unused*/,
                              std::vector<vertex_descriptor>& sources, std::vector<vertex_descriptor>& targets,
                              std::vector<typename forward_type::inherited_edge_properties::edge_bundled>& edge_props,
                              std::size_t numlocalverts, GlobalToLocal global_to_local,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), property_(prop) {
    forward_.assign_sources_and_targets_global(sources, targets, edge_props, numlocalverts, global_to_local);
  }

  //  From number of vertices and single-pass range of unsorted edges.  Data
  //  is cached in coordinate form before creating the actual graph.
  template <std::ranges::input_range EdgeRange>
  compressed_sparse_row_graph(edges_are_unsorted_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), property_(prop) {
    std::vector<vertex_descriptor> sources;
    std::vector<vertex_descriptor> targets;
    histogram_detail::split_into_separate_coords(edge_rg, sources, targets);
    forward_.assign_sources_and_targets_global(sources, targets, numverts, typed_identity_property_map<std::size_t>());
  }

  //  From number of vertices and single-pass range of unsorted edges and
  //  single-pass range of edge properties.  Data is cached in coordinate form
  //  before creating the actual graph.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_unsorted_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), property_(prop) {
    std::vector<vertex_descriptor> sources;
    std::vector<vertex_descriptor> targets;
    histogram_detail::split_into_separate_coords(edge_rg, sources, targets);
    std::size_t numedges = sources.size();
    std::vector<typename forward_type::inherited_edge_properties::edge_bundled> edge_props(numedges);
    auto ep_iter = ep_rg.begin();
    for (std::size_t i = 0; i < numedges; ++i) {
      edge_props[i] = *ep_iter++;
    }
    forward_.assign_sources_and_targets_global(sources, targets, edge_props, numverts,
                                               typed_identity_property_map<std::size_t>());
  }

  //  From number of vertices and single-pass range of unsorted edges.  Data
  //  is cached in coordinate form before creating the actual graph.  Edges
  //  are filtered and transformed for use in a distributed graph.
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_global_t /*unused*/, EdgeRange edge_rg, std::size_t numlocalverts,
                              GlobalToLocal global_to_local, const SourcePred& source_pred,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), property_(prop) {
    std::vector<vertex_descriptor> sources;
    std::vector<vertex_descriptor> targets;
    histogram_detail::split_into_separate_coords_filtered(edge_rg, sources, targets, source_pred);
    forward_.assign_sources_and_targets_global(sources, targets, numlocalverts, global_to_local);
  }

  //  From number of vertices and single-pass range of unsorted edges and
  //  single-pass range of edge properties.  Data is cached in coordinate form
  //  before creating the actual graph.  Edges are filtered and transformed
  //  for use in a distributed graph.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange, typename GlobalToLocal,
            typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_global_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numlocalverts, GlobalToLocal global_to_local, const SourcePred& source_pred,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), property_(prop) {
    std::vector<vertex_descriptor> sources;
    std::vector<vertex_descriptor> targets;
    std::vector<edge_bundled> edge_props;
    histogram_detail::split_into_separate_coords_filtered(edge_rg, ep_rg, sources, targets, edge_props, source_pred);
    forward_.assign_sources_and_targets_global(sources, targets, edge_props, numlocalverts, global_to_local);
  }

  //   Requires IncidenceGraph and a vertex index map
  template <typename Graph, typename VertexIndexMap>
  compressed_sparse_row_graph(const Graph& g, const VertexIndexMap& vi, std::size_t numverts, std::size_t numedges)
      : property_() {
    assign(g, vi, numverts, numedges);
    inherited_vertex_properties::resize(numverts);
  }

  //   Requires VertexListGraph and EdgeListGraph
  template <typename Graph, typename VertexIndexMap>
  compressed_sparse_row_graph(const Graph& g, const VertexIndexMap& vi) : property_() {
    std::size_t numedges = num_edges(g);
    if constexpr (is_undirected_graph_v<Graph>) {
      numedges *= 2;  // Double each edge (actual doubling done by
                      // out_edges function)
    }
    std::size_t numverts = num_vertices(g);
    assign(g, vi, numverts, numedges);
    inherited_vertex_properties::resize(numverts);
  }

  // Requires vertex index map plus requirements of previous constructor
  template <typename Graph>
  explicit compressed_sparse_row_graph(const Graph& g) : compressed_sparse_row_graph(g, get(vertex_index, g)) {}

  // From any graph (slow and uses a lot of memory)
  //   Requires IncidenceGraph and a vertex index map
  //   Internal helper function
  //   Note that numedges must be doubled for undirected source graphs
  template <typename Graph, typename VertexIndexMap>
  void assign(const Graph& g, const VertexIndexMap& vi, std::size_t numverts, std::size_t numedges) {
    forward_.assign(g, vi, numverts, numedges);
    inherited_vertex_properties::resize(numverts);
  }

  // Requires the above, plus VertexListGraph and EdgeListGraph
  template <typename Graph, typename VertexIndexMap>
  void assign(const Graph& g, const VertexIndexMap& vi) {
    std::size_t numedges = num_edges(g);
    if constexpr (is_undirected_graph_v<Graph>) {
      numedges *= 2;  // Double each edge (actual doubling done by
                      // out_edges function)
    }
    assign(g, vi, num_vertices(g), numedges);
  }

  // Requires the above, plus a vertex_index map.
  template <typename Graph>
  void assign(const Graph& g) {
    assign(g, get(vertex_index, g));
  }

  // Add edges from a sorted (smallest sources first) range of pairs and edge
  // properties
  template <std::ranges::bidirectional_range BidirRange, std::ranges::bidirectional_range EPRange,
            typename GlobalToLocal>
  void add_edges_sorted_internal(BidirRange sorted_rg, EPRange ep_sorted_rg, const GlobalToLocal& global_to_local) {
    forward_.add_edges_sorted_internal(sorted_rg, ep_sorted_rg, global_to_local);
  }

  template <std::ranges::bidirectional_range BidirRange, std::ranges::bidirectional_range EPRange>
  void add_edges_sorted_internal(BidirRange sorted_rg, EPRange ep_sorted_rg) {
    forward_.add_edges_sorted_internal(sorted_rg, ep_sorted_rg, typed_identity_property_map<std::size_t>());
  }

  // Add edges from a sorted (smallest sources first) range of pairs
  template <std::ranges::bidirectional_range BidirRange>
  void add_edges_sorted_internal(BidirRange sorted_rg) {
    forward_.add_edges_sorted_internal(sorted_rg, std::ranges::single_view<edge_bundled>());
  }

  template <std::ranges::bidirectional_range BidirRange, typename GlobalToLocal>
  void add_edges_sorted_internal_global(BidirRange sorted_rg, const GlobalToLocal& global_to_local) {
    forward_.add_edges_sorted_internal(sorted_rg, std::ranges::single_view<edge_bundled>(), global_to_local);
  }

  template <std::ranges::bidirectional_range BidirRange, std::ranges::bidirectional_range EPRange,
            typename GlobalToLocal>
  void add_edges_sorted_internal_global(BidirRange sorted_rg, EPRange ep_sorted_rg,
                                        const GlobalToLocal& global_to_local) {
    forward_.add_edges_sorted_internal(sorted_rg, ep_sorted_rg, global_to_local);
  }

  // Add edges from a range of (source, target) pairs that are unsorted
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal>
  void add_edges_internal(EdgeRange edge_rg, const GlobalToLocal& global_to_local) {
    std::vector<std::pair<Vertex, Vertex>> new_edges(edge_rg.begin(), edge_rg.end());
    if (new_edges.empty()) {
      return;
    }
    std::sort(new_edges.begin(), new_edges.end());
    add_edges_sorted_internal_global(std::ranges::ref_view(new_edges), global_to_local);
  }

  template <std::ranges::input_range EdgeRange>
  void add_edges_internal(EdgeRange edge_rg) {
    add_edges_internal(edge_rg, typed_identity_property_map<std::size_t>());
  }

  // Add edges from a range of (source, target) pairs and edge properties that
  // are unsorted
  template <std::ranges::input_range InRange, std::ranges::input_range EPRange, typename GlobalToLocal>
  void add_edges_internal(InRange uv_range, EPRange ep_range, const GlobalToLocal& global_to_local) {
    using vertex_pair = std::pair<Vertex, Vertex>;
    using edge_vector_t = std::vector<std::tuple<vertex_pair, edge_bundled>>;
    edge_vector_t new_edges;
    if constexpr (std::ranges::random_access_range<InRange>) {
      new_edges.reserve(uv_range.size());
    }
    for (auto [uv, ep] : zip_range(uv_range, ep_range)) {
      new_edges.emplace_back(uv, ep);
    }
    if (new_edges.empty()) {
      return;
    }
    std::sort(new_edges.begin(), new_edges.end(),
              [](const auto& lhs, const auto& rhs) { return std::get<0>(lhs) < std::get<0>(rhs); });
    forward_.add_edges_sorted_internal(std::views::transform(new_edges, [](const auto& x) { return std::get<0>(x); }),
                                       std::views::transform(new_edges, [](const auto& x) { return std::get<1>(x); }),
                                       global_to_local);
  }

  // Add edges from a range of (source, target) pairs and edge properties that
  // are unsorted
  template <std::ranges::input_range InRange, std::ranges::input_range EPRange>
  void add_edges_internal(InRange uv_range, EPRange ep_range) {
    add_edges_internal(uv_range, ep_range, typed_identity_property_map<std::size_t>());
  }

  using inherited_vertex_properties::operator[];

  // Directly access a edge or edge bundle
  edge_push_back_type& operator[](const edge_descriptor& v) {
    return forward_.edge_properties_[get(edge_index, *this, v)];
  }

  const edge_push_back_type& operator[](const edge_descriptor& v) const {
    return forward_.edge_properties_[get(edge_index, *this, v)];
  }

  // Directly access a graph bundle
  graph_bundled& operator[](graph_bundle_t /*unused*/) { return get_property(*this); }

  const graph_bundled& operator[](graph_bundle_t /*unused*/) const { return get_property(*this); }

  // private: non-portable, requires friend templates
  inherited_vertex_properties& vertex_properties() { return *this; }
  const inherited_vertex_properties& vertex_properties() const { return *this; }
  auto& edge_properties() { return forward_; }
  const auto& edge_properties() const { return forward_; }

  forward_type forward_;
  GraphProperty property_;
};

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty, typename Vertex, typename EdgeIndex>
class compressed_sparse_row_graph<bidirectional_s, VertexProperty, EdgeProperty, GraphProperty, Vertex, EdgeIndex>
    : public indexed_detail::indexed_vertex_properties<BAGL_BIDIR_CSR_GRAPH_TYPE, VertexProperty, Vertex,
                                                       typed_identity_property_map<Vertex>> {
 public:
  using inherited_vertex_properties =
      indexed_detail::indexed_vertex_properties<compressed_sparse_row_graph, VertexProperty, Vertex,
                                                typed_identity_property_map<Vertex>>;

  // For Property Graph
  using graph_property_type = GraphProperty;
  using graph_bundled = lookup_one_property_t<GraphProperty, graph_bundle_t>;
  // typedef GraphProperty graph_property_type;

  using forward_type = csr_detail::compressed_sparse_row_structure<EdgeProperty, Vertex, EdgeIndex>;
  using backward_edge_property =
      EdgeIndex /* std::conditional_t<std::is_same_v<EdgeProperty, no_property>, no_property, EdgeIndex> */;
  using backward_type = csr_detail::compressed_sparse_row_structure<backward_edge_property, Vertex, EdgeIndex>;

  // Concept requirements:
  // For Graph
  using vertex_descriptor = Vertex;
  using edge_descriptor = csr_detail::csr_edge_descriptor<Vertex, EdgeIndex>;
  using directed_category = bidirectional_tag;
  using edge_parallel_category = allow_parallel_edge_tag;

  struct traversal_category : bidirectional_graph_tag,
                              adjacency_graph_tag,
                              vertex_list_graph_tag,
                              edge_list_graph_tag {};

  static vertex_descriptor null_vertex() { return std::numeric_limits<vertex_descriptor>::max(); }

  // For VertexListGraph
  using vertices_size_type = Vertex;

  // For EdgeListGraph
  using edges_size_type = EdgeIndex;

  // For IncidenceGraph
  using degree_size_type = EdgeIndex;

  // For internal use
  using graph_tag = csr_graph_tag;

  using edge_bundled = typename forward_type::inherited_edge_properties::edge_bundled;
  using edge_push_back_type = typename forward_type::inherited_edge_properties::edge_push_back_type;
  using edge_property_type = typename forward_type::inherited_edge_properties::edge_property_type;

  // Constructors

  // Default constructor: an empty graph.
  compressed_sparse_row_graph() : property_() {}

  //  With numverts vertices
  explicit compressed_sparse_row_graph(std::size_t numverts)
      : inherited_vertex_properties(numverts), forward_(numverts), backward_(numverts) {}

 private:
  void set_up_backward_property_links() {
    auto e_rg = edges(*this);
    auto uv_rg = csr_detail::make_edge_to_index_pair_range(*this, get(vertex_index, *this), e_rg);
    backward_.assign_unsorted_multi_pass_edges(uv_rg, std::views::iota(0, num_edges(*this)),
                                               forward_.rowstart_.size() - 1, typed_identity_property_map<Vertex>(),
                                               keep_all());
  }

 public:
  //  From number of vertices and unsorted list of edges
  template <std::ranges::forward_range EdgeRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numverts, typed_identity_property_map<Vertex>(), keep_all());
    set_up_backward_property_links();
  }

  //  From number of vertices and unsorted list of edges, plus edge properties
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numverts, typed_identity_property_map<Vertex>(),
                                              keep_all());
    set_up_backward_property_links();
  }

  //  From number of vertices and unsorted list of edges, with filter and
  //  global-to-local map
  template <std::ranges::forward_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numlocalverts, global_to_local, source_pred);
    set_up_backward_property_links();
  }

  //  From number of vertices and unsorted list of edges, plus edge
  //  properties, with filter and global-to-local map
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange, typename GlobalToLocal,
            typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : inherited_vertex_properties(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numlocalverts, global_to_local, source_pred);
    set_up_backward_property_links();
  }

  //   Requires IncidenceGraph and a vertex index map
  template <typename Graph, typename VertexIndexMap>
  compressed_sparse_row_graph(const Graph& g, const VertexIndexMap& vi, std::size_t numverts, std::size_t numedges)
      : property_() {
    assign(g, vi, numverts, numedges);
    inherited_vertex_properties::resize(numverts);
  }

  //   Requires VertexListGraph and EdgeListGraph
  template <typename Graph, typename VertexIndexMap>
  compressed_sparse_row_graph(const Graph& g, const VertexIndexMap& vi) : property_() {
    std::size_t numedges = num_edges(g);
    if constexpr (is_undirected_graph_v<Graph>) {
      numedges *= 2;  // Double each edge (actual doubling done by
                      // out_edges function)
    }
    std::size_t numverts = num_vertices(g);
    assign(g, vi, numverts, numedges);
    inherited_vertex_properties::resize(numverts);
  }

  // Requires vertex index map plus requirements of previous constructor
  template <typename Graph>
  explicit compressed_sparse_row_graph(const Graph& g) : compressed_sparse_row_graph(g, get(vertex_index, g)) {}

  // From any graph (slow and uses a lot of memory)
  //   Requires IncidenceGraph and a vertex index map
  //   Internal helper function
  //   Note that numedges must be doubled for undirected source graphs
  template <typename Graph, typename VertexIndexMap>
  void assign(const Graph& g, const VertexIndexMap& vi, vertices_size_type numverts, edges_size_type numedges) {
    forward_.assign(g, vi, numverts, numedges);
    inherited_vertex_properties::resize(numverts);
    set_up_backward_property_links();
  }

  // Requires the above, plus VertexListGraph and EdgeListGraph
  template <typename Graph, typename VertexIndexMap>
  void assign(const Graph& g, const VertexIndexMap& vi) {
    std::size_t numedges = num_edges(g);
    if constexpr (is_undirected_graph_v<Graph>) {
      numedges *= 2;  // Double each edge (actual doubling done by
                      // out_edges function)
    }
    assign(g, vi, num_vertices(g), numedges);
  }

  // Requires the above, plus a vertex_index map.
  template <typename Graph>
  void assign(const Graph& g) {
    assign(g, get(vertex_index, g));
  }

  using inherited_vertex_properties::operator[];

  // Directly access a edge or edge bundle
  edge_push_back_type& operator[](const edge_descriptor& v) {
    return forward_.edge_properties_[get(edge_index, *this, v)];
  }

  const edge_push_back_type& operator[](const edge_descriptor& v) const {
    return forward_.edge_properties_[get(edge_index, *this, v)];
  }

  // private: non-portable, requires friend templates
  inherited_vertex_properties& vertex_properties() { return *this; }
  const inherited_vertex_properties& vertex_properties() const { return *this; }
  auto& edge_properties() { return forward_; }
  const auto& edge_properties() const { return forward_; }

  forward_type forward_;
  backward_type backward_;
  GraphProperty property_;
};

// Construction functions
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex add_vertex(BAGL_CSR_GRAPH_TYPE& g) {
  add_vertex(g, typename BAGL_CSR_GRAPH_TYPE::vertex_bundled());
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
Vertex add_vertex(BAGL_DIR_CSR_GRAPH_TYPE& g, typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_bundled const& p) {
  Vertex old_num_verts_plus_one = g.forward_.rowstart_.size();
  g.forward_.rowstart_.push_back(g.forward_.rowstart_.back());
  g.vertex_properties().push_back(p);
  return old_num_verts_plus_one - 1;
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
Vertex add_vertex(BAGL_BIDIR_CSR_GRAPH_TYPE& g, typename BAGL_BIDIR_CSR_GRAPH_TYPE::vertex_bundled const& p) {
  Vertex old_num_verts_plus_one = g.forward_.rowstart_.size();
  g.forward_.rowstart_.push_back(g.forward_.rowstart_.back());
  g.backward_.rowstart_.push_back(g.backward_.rowstart_.back());
  g.vertex_properties().push_back(p);
  return old_num_verts_plus_one - 1;
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
Vertex add_vertices(typename BAGL_DIR_CSR_GRAPH_TYPE::vertices_size_type count, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  Vertex old_num_verts_plus_one = g.forward_.rowstart_.size();
  EdgeIndex numedges = g.forward_.rowstart_.back();
  g.forward_.rowstart_.resize(old_num_verts_plus_one + count, numedges);
  g.vertex_properties().resize(num_vertices(g));
  return old_num_verts_plus_one - 1;
}

// Add edges from a sorted (smallest sources first) range of pairs and edge
// properties
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename BidirectionalIteratorOrig, typename EPIterOrig>
void add_edges_sorted(BidirectionalIteratorOrig first_sorted, BidirectionalIteratorOrig last_sorted,
                      EPIterOrig ep_iter_sorted, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_sorted_internal(first_sorted, last_sorted, ep_iter_sorted);
}

// Add edges from a sorted (smallest sources first) range of pairs
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename BidirectionalIteratorOrig>
void add_edges_sorted(BidirectionalIteratorOrig first_sorted, BidirectionalIteratorOrig last_sorted,
                      BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_sorted_internal(first_sorted, last_sorted);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename BidirectionalIteratorOrig, typename EPIterOrig,
          typename GlobalToLocal>
void add_edges_sorted_global(BidirectionalIteratorOrig first_sorted, BidirectionalIteratorOrig last_sorted,
                             EPIterOrig ep_iter_sorted, const GlobalToLocal& global_to_local,
                             BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_sorted_internal_global(first_sorted, last_sorted, ep_iter_sorted, global_to_local);
}

// Add edges from a sorted (smallest sources first) range of pairs
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename BidirectionalIteratorOrig, typename GlobalToLocal>
void add_edges_sorted_global(BidirectionalIteratorOrig first_sorted, BidirectionalIteratorOrig last_sorted,
                             const GlobalToLocal& global_to_local, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_sorted_internal_global(first_sorted, last_sorted, global_to_local);
}

// Add edges from a range of (source, target) pairs that are unsorted
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename InputIterator, typename GlobalToLocal>
void add_edges_global(InputIterator first, InputIterator last, const GlobalToLocal& global_to_local,
                      BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_internal(first, last, global_to_local);
}

// Add edges from a range of (source, target) pairs that are unsorted
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename InputIterator>
void add_edges(InputIterator first, InputIterator last, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_internal(first, last);
}

// Add edges from a range of (source, target) pairs and edge properties that
// are unsorted
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename InputIterator, typename EPIterator>
void add_edges(InputIterator first, InputIterator last, EPIterator ep_iter, EPIterator ep_iter_end,
               BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_internal(first, last, ep_iter, ep_iter_end);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename InputIterator, typename EPIterator, typename GlobalToLocal>
void add_edges_global(InputIterator first, InputIterator last, EPIterator ep_iter, EPIterator ep_iter_end,
                      const GlobalToLocal& global_to_local, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  g.add_edges_internal(first, last, ep_iter, ep_iter_end, global_to_local);
}

// From VertexListGraph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex num_vertices(const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.rowstart_.size() - 1;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto vertices(const BAGL_CSR_GRAPH_TYPE& g) {
  return std::views::iota(Vertex{0}, Vertex{num_vertices(g)});
}

// From IncidenceGraph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex source(typename BAGL_CSR_GRAPH_TYPE::edge_descriptor e, const BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return e.src;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex target(typename BAGL_CSR_GRAPH_TYPE::edge_descriptor e, const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.column_[e.idx];
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto out_edges(Vertex v, const BAGL_CSR_GRAPH_TYPE& g) {
  return csr_detail::make_csr_out_edge_range<BAGL_CSR_GRAPH_TYPE>(g.forward_.rowstart_[v], g.forward_.rowstart_[v + 1],
                                                                  v);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
EdgeIndex out_degree(Vertex v, const BAGL_CSR_GRAPH_TYPE& g) {
  EdgeIndex v_row_start = g.forward_.rowstart_[v];
  EdgeIndex next_row_start = g.forward_.rowstart_[v + 1];
  return next_row_start - v_row_start;
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
auto in_edges(Vertex v, const BAGL_BIDIR_CSR_GRAPH_TYPE& g) {
  return csr_detail::make_csr_in_edge_range(g.backward_.rowstart_[v], g.backward_.rowstart_[v + 1], g);
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
EdgeIndex in_degree(Vertex v, const BAGL_BIDIR_CSR_GRAPH_TYPE& g) {
  EdgeIndex v_row_start = g.backward_.rowstart_[v];
  EdgeIndex next_row_start = g.backward_.rowstart_[v + 1];
  return next_row_start - v_row_start;
}

// From AdjacencyGraph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto adjacent_vertices(Vertex v, const BAGL_CSR_GRAPH_TYPE& g) {
  EdgeIndex v_row_start = g.forward_.rowstart_[v];
  EdgeIndex next_row_start = g.forward_.rowstart_[v + 1];
  return std::ranges::subrange(g.forward_.column_.begin() + v_row_start, g.forward_.column_.begin() + next_row_start);
}

// Extra, common functions
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex vertex(Vertex i, const BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return i;
}

// edge() can be provided in linear time for the new interface

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
std::pair<typename BAGL_CSR_GRAPH_TYPE::edge_descriptor, bool> edge(Vertex i, Vertex j, const BAGL_CSR_GRAPH_TYPE& g) {
  for (auto e : out_edges(i, g)) {
    if (target(e, g) == j) {
      return {e, true};
    }
  }
  return {{}, false};
}

// Find an edge given its index in the graph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto edge_from_index(EdgeIndex idx, const BAGL_CSR_GRAPH_TYPE& g) {
  assert(idx < num_edges(g));
  auto src_plus_1 = std::upper_bound(g.forward_.rowstart_.begin(), g.forward_.rowstart_.end(), idx);
  // Get last source whose rowstart is at most idx
  // upper_bound returns this position plus 1
  Vertex src = (src_plus_1 - g.forward_.rowstart_.begin()) - 1;
  return typename BAGL_CSR_GRAPH_TYPE::edge_descriptor(src, idx);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
EdgeIndex num_edges(const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.column_.size();
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto edges(const BAGL_CSR_GRAPH_TYPE& g) {
  return edges_from_out_edges(g);
}

// For Property Graph

// Graph properties
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, class Tag, class Value>
void set_property(BAGL_CSR_GRAPH_TYPE& g, Tag tag, const Value& value) {
  get_property_value(g.property_, tag) = value;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, class Tag>
auto& get_property(BAGL_CSR_GRAPH_TYPE& g, Tag tag) {
  return get_property_value(g.property_, tag);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, class Tag>
const auto& get_property(const BAGL_CSR_GRAPH_TYPE& g, Tag tag) {
  return get_property_value(g.property_, tag);
}

template <typename G, typename Tag, typename Kind>
struct csr_property_map_helper {};
// Kind == void for invalid property tags, so we can use that to SFINAE out

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
struct csr_property_map_helper<BAGL_CSR_GRAPH_TYPE, Tag, vertex_property_tag> {
  using all_tag = vertex_all_t;
  using key_type = property_traits_key_t<property_map_t<BAGL_CSR_GRAPH_TYPE, vertex_all_t>>;
  using plist_type = VertexProperty;
  using all_type = property_map_t<BAGL_CSR_GRAPH_TYPE, vertex_all_t>;
  using all_const_type = property_map_const_t<BAGL_CSR_GRAPH_TYPE, vertex_all_t>;
  using type = transform_value_property_map<property_detail::lookup_one_property_f<plist_type, Tag>, all_type>;
  using const_type =
      transform_value_property_map<property_detail::lookup_one_property_f<const plist_type, Tag>, all_const_type>;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
struct csr_property_map_helper<BAGL_CSR_GRAPH_TYPE, Tag, edge_property_tag> {
  using all_tag = edge_all_t;
  using key_type = property_traits_key_t<property_map_t<BAGL_CSR_GRAPH_TYPE, edge_all_t>>;
  using plist_type = EdgeProperty;
  using all_type = property_map_t<BAGL_CSR_GRAPH_TYPE, edge_all_t>;
  using all_const_type = property_map_const_t<BAGL_CSR_GRAPH_TYPE, edge_all_t>;
  using type = transform_value_property_map<property_detail::lookup_one_property_f<plist_type, Tag>, all_type>;
  using const_type =
      transform_value_property_map<property_detail::lookup_one_property_f<const plist_type, Tag>, all_const_type>;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
struct csr_property_map_helper<BAGL_CSR_GRAPH_TYPE, Tag, graph_property_tag> {
  using all_tag = graph_all_t;
  using key_type = BAGL_CSR_GRAPH_TYPE*;
  using plist_type = GraphProperty;
  using all_type = property_map_t<BAGL_CSR_GRAPH_TYPE, graph_all_t>;
  using all_const_type = property_map_const_t<BAGL_CSR_GRAPH_TYPE, graph_all_t>;
  using type = transform_value_property_map<property_detail::lookup_one_property_f<plist_type, Tag>, all_type>;
  using const_type =
      transform_value_property_map<property_detail::lookup_one_property_f<const plist_type, Tag>, all_const_type>;
};

// disable_if isn't truly necessary but required to avoid ambiguity with
// specializations below
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
struct property_map<BAGL_CSR_GRAPH_TYPE, Tag, std::enable_if_t<!distributed_detail::is_distributed_selector_v<Vertex>>>
    : csr_property_map_helper<BAGL_CSR_GRAPH_TYPE, Tag,
                              typename properties_detail::property_kind_from_graph<BAGL_CSR_GRAPH_TYPE, Tag>::type> {};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
decltype(auto) get(Tag tag, BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_t<BAGL_CSR_GRAPH_TYPE, Tag>(tag,
                                                  get(typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::all_tag(), g));
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
decltype(auto) get(Tag tag, const BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_const_t<BAGL_CSR_GRAPH_TYPE, Tag>(
      tag, get(typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::all_tag(), g));
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
decltype(auto) get(Tag tag, BAGL_CSR_GRAPH_TYPE& g, typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::key_type k) {
  using all_tag = typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::all_tag;
  using outer_pm = typename property_map<BAGL_CSR_GRAPH_TYPE, all_tag>::type;
  return lookup_one_property<property_traits_value_t<outer_pm>, Tag>::lookup(get(all_tag(), g, k), tag);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
decltype(auto) get(Tag tag, const BAGL_CSR_GRAPH_TYPE& g, typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::key_type k) {
  using all_tag = typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::all_tag;
  using outer_pm = typename property_map<BAGL_CSR_GRAPH_TYPE, all_tag>::type;
  return lookup_one_property<const property_traits_value_t<outer_pm>, Tag>::lookup(get(all_tag(), g, k), tag);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
void put(Tag tag, BAGL_CSR_GRAPH_TYPE& g, typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::key_type k,
         lookup_one_property_t<typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::plist_type, Tag> val) {
  using all_tag = typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::all_tag;
  lookup_one_property<typename property_map<BAGL_CSR_GRAPH_TYPE, Tag>::plist_type, Tag>::lookup(get(all_tag(), g, k),
                                                                                                tag) = val;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, vertex_index_t,
                    std::enable_if_t<!distributed_detail::is_distributed_selector_v<Vertex>>> {
  using type = typed_identity_property_map<Vertex>;
  using const_type = type;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, edge_index_t,
                    std::enable_if_t<!distributed_detail::is_distributed_selector_v<Vertex>>> {
  using type = csr_detail::csr_edge_index_map<Vertex, EdgeIndex>;
  using const_type = type;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, vertex_all_t,
                    std::enable_if_t<!distributed_detail::is_distributed_selector_v<Vertex>>> {
  using type = typename BAGL_CSR_GRAPH_TYPE::inherited_vertex_properties::vertex_map_type;
  using const_type = typename BAGL_CSR_GRAPH_TYPE::inherited_vertex_properties::const_vertex_map_type;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, edge_all_t,
                    std::enable_if_t<!distributed_detail::is_distributed_selector_v<Vertex>>> {
  using type = typename BAGL_CSR_GRAPH_TYPE::forward_type::inherited_edge_properties::edge_map_type;
  using const_type = typename BAGL_CSR_GRAPH_TYPE::forward_type::inherited_edge_properties::const_edge_map_type;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, graph_all_t,
                    std::enable_if_t<!distributed_detail::is_distributed_selector_v<Vertex>>> {
  using type = ref_property_map<BAGL_CSR_GRAPH_TYPE*, typename BAGL_CSR_GRAPH_TYPE::graph_property_type>;
  using const_type = ref_property_map<BAGL_CSR_GRAPH_TYPE*, const typename BAGL_CSR_GRAPH_TYPE::graph_property_type>;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return typed_identity_property_map<Vertex>();
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex get(vertex_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*unused*/, Vertex v) {
  return v;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return typed_identity_property_map<Vertex>();
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
Vertex get(vertex_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*unused*/, Vertex v) {
  return v;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(edge_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return property_map_const_t<BAGL_CSR_GRAPH_TYPE, edge_index_t>{};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
EdgeIndex get(edge_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*unused*/,
              typename BAGL_CSR_GRAPH_TYPE::edge_descriptor e) {
  return e.idx;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(edge_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return property_map_const_t<BAGL_CSR_GRAPH_TYPE, edge_index_t>{};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
EdgeIndex get(edge_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*unused*/,
              typename BAGL_CSR_GRAPH_TYPE::edge_descriptor e) {
  return e.idx;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g) {
  return g.get_vertex_bundle(get(vertex_index, g));
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_all_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g) {
  return g.get_vertex_bundle(get(vertex_index, g));
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
VertexProperty& get(vertex_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g, Vertex v) {
  return get(vertex_all, g)[v];
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const VertexProperty& get(vertex_all_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g, Vertex v) {
  return get(vertex_all, g)[v];
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
void put(vertex_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g, Vertex v, const VertexProperty& val) {
  put(get(vertex_all, g), v, val);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(edge_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.get_edge_bundle(get(edge_index, g));
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(edge_all_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.get_edge_bundle(get(edge_index, g));
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
EdgeProperty& get(edge_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g,
                  const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e) {
  return get(edge_all, g)[e];
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const EdgeProperty& get(edge_all_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g,
                        const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e) {
  return get(edge_all, g)[e];
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
void put(edge_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g, const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e,
         const EdgeProperty& val) {
  put(get(edge_all, g), e, val);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(graph_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_t<BAGL_CSR_GRAPH_TYPE, graph_all_t>(g.property_);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(graph_all_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_const_t<BAGL_CSR_GRAPH_TYPE, graph_all_t>(g.property_);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
GraphProperty& get(graph_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g, BAGL_CSR_GRAPH_TYPE* /*unused*/) {
  return g.property_;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const GraphProperty& get(graph_all_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g, BAGL_CSR_GRAPH_TYPE* /*unused*/) {
  return g.property_;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
void put(graph_all_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g, BAGL_CSR_GRAPH_TYPE* /*unused*/, const GraphProperty& val) {
  g.property_ = val;
}

#undef BAGL_CSR_GRAPH_TYPE
#undef BAGL_CSR_GRAPH_TEMPLATE_PARMS
#undef BAGL_DIR_CSR_GRAPH_TYPE
#undef BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS
#undef BAGL_BIDIR_CSR_GRAPH_TYPE
#undef BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS

}  // namespace bagl

#endif  // BAGL_BAGL_COMPRESSED_SPARSE_ROW_GRAPH_H_
