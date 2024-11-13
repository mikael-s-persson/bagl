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
#include "bagl/detail/container_generators.h"
#include "bagl/detail/is_distributed_selector.h"
#include "bagl/filtered_graph.h"  // For keep_all
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
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
#define BAGL_CSR_GRAPH_TEMPLATE_PARMS \
  typename Directed, typename VertexProperty, typename EdgeProperty, typename GraphProperty
#define BAGL_CSR_GRAPH_TYPE compressed_sparse_row_graph<Directed, VertexProperty, EdgeProperty, GraphProperty>
#define BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS typename VertexProperty, typename EdgeProperty, typename GraphProperty
#define BAGL_DIR_CSR_GRAPH_TYPE compressed_sparse_row_graph<directed_s, VertexProperty, EdgeProperty, GraphProperty>
#define BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS typename VertexProperty, typename EdgeProperty, typename GraphProperty
#define BAGL_BIDIR_CSR_GRAPH_TYPE \
  compressed_sparse_row_graph<bidirectional_s, VertexProperty, EdgeProperty, GraphProperty>

/** Compressed sparse row graph.
 */
template <typename Directed = directed_s, typename VertexProperty = no_property, typename EdgeProperty = no_property,
          typename GraphProperty = no_property>
class compressed_sparse_row_graph;  // Not defined

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty>
class compressed_sparse_row_graph<directed_s, VertexProperty, EdgeProperty, GraphProperty> {
 public:
  // Some tests to prevent use of "void" is a property type (as was done in
  // some test cases):
  static_assert(!std::is_same_v<VertexProperty, void>);
  static_assert(!std::is_same_v<EdgeProperty, void>);
  static_assert(!std::is_same_v<GraphProperty, void>);

  // For Property Graph
  using graph_property_type = GraphProperty;
  using graph_bundled = lookup_one_property_t<GraphProperty, graph_bundle_t>;
  using forward_type = csr_detail::compressed_sparse_row_structure<EdgeProperty>;

  /* At this time, the compressed sparse row graph can only be used to
   * create directed and bidirectional graphs. In the future,
   * undirected CSR graphs will also be supported.
   */
  // static_assert(std::is_same_v<Directed, directed_s>);

  // Concept requirements:
  // For Graph
  using vertex_descriptor = std::size_t;
  using edge_descriptor = csr_detail::csr_edge_descriptor;
  using directed_category = directed_tag;
  using edge_parallel_category = allow_parallel_edge_tag;

  struct traversal_category : incidence_graph_tag, adjacency_graph_tag, vertex_list_graph_tag, edge_list_graph_tag {};

  static vertex_descriptor null_vertex() { return std::numeric_limits<vertex_descriptor>::max(); }

  // For VertexListGraph
  using vertices_size_type = std::size_t;

  // For EdgeListGraph
  using edges_size_type = std::size_t;

  // For IncidenceGraph
  using degree_size_type = std::size_t;

  // For internal use
  using graph_tag = csr_graph_tag;

  using edge_bundled = typename forward_type::edge_bundled;
  using edge_property_type = typename forward_type::edge_property_type;

  using vertex_property_type = VertexProperty;
  using vertex_bundled = lookup_one_property_t<VertexProperty, vertex_bundle_t>;

  // Constructors

  // Default constructor: an empty graph.
  compressed_sparse_row_graph() : property_() {}

  //  With numverts vertices
  explicit compressed_sparse_row_graph(std::size_t numverts) : vertex_properties_(numverts), forward_(numverts) {}

  //  From number of vertices and unsorted list of edges
  template <std::ranges::forward_range EdgeRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numverts, typed_identity_property_map<std::size_t>(),
                                              keep_all());
  }

  //  From number of vertices and unsorted list of edges, plus edge properties
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numverts, typed_identity_property_map<std::size_t>(),
                                              keep_all());
  }

  //  From number of vertices and unsorted list of edges, with filter and
  //  global-to-local map
  template <std::ranges::forward_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numlocalverts, global_to_local, source_pred);
  }

  //  From number of vertices and unsorted list of edges, plus edge
  //  properties, with filter and global-to-local map
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange, typename GlobalToLocal,
            typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numlocalverts, global_to_local, source_pred);
  }

  //  From number of vertices and sorted list of edges (new interface)
  template <std::ranges::input_range EdgeRange>
  compressed_sparse_row_graph(edges_are_sorted_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              std::size_t numedges = 0, const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, typed_identity_property_map<std::size_t>(), keep_all(), numverts,
                                      numedges);
    vertex_properties_.resize(numverts);
  }

  //  From number of vertices and sorted list of edges (new interface)
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_sorted_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, std::size_t numedges = 0,
                              const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, ep_rg, typed_identity_property_map<std::size_t>(), keep_all(), numverts,
                                      numedges);
    vertex_properties_.resize(numverts);
  }

  //  From number of vertices and sorted list of edges, filtered and global
  //  (new interface)
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_sorted_global_t /*unused*/, EdgeRange edge_rg,
                              const GlobalToLocal& global_to_local, const SourcePred& source_pred, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, global_to_local, source_pred, numverts, 0);
    vertex_properties_.resize(numverts);
  }

  //  From number of vertices and sorted list of edges (new interface)
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange, typename GlobalToLocal,
            typename SourcePred>
  compressed_sparse_row_graph(edges_are_sorted_global_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              const GlobalToLocal& global_to_local, const SourcePred& source_pred, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : property_(prop) {
    forward_.assign_from_sorted_edges(edge_rg, ep_rg, global_to_local, source_pred, numverts, 0);
    vertex_properties_.resize(numverts);
  }

  //  From number of vertices and mutable vectors of sources and targets;
  //  vectors are returned with unspecified contents but are guaranteed not to
  //  share storage with the constructed graph.
  compressed_sparse_row_graph(construct_inplace_from_sources_and_targets_t /*unused*/,
                              std::vector<vertex_descriptor>& sources, std::vector<vertex_descriptor>& targets,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), property_(prop) {
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
      : vertex_properties_(numlocalverts), property_(prop) {
    forward_.assign_sources_and_targets_global(sources, targets, numlocalverts, global_to_local);
  }

  //  From number of vertices and mutable vectors of sources, targets, and
  //  edge properties; vectors are returned with unspecified contents but are
  //  guaranteed not to share storage with the constructed graph.
  compressed_sparse_row_graph(construct_inplace_from_sources_and_targets_t /*unused*/,
                              std::vector<vertex_descriptor>& sources, std::vector<vertex_descriptor>& targets,
                              std::vector<edge_property_type>& edge_props, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), property_(prop) {
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
                              std::vector<edge_property_type>& edge_props, std::size_t numlocalverts,
                              GlobalToLocal global_to_local, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numlocalverts), property_(prop) {
    forward_.assign_sources_and_targets_global(sources, targets, edge_props, numlocalverts, global_to_local);
  }

  //  From number of vertices and single-pass range of unsorted edges.  Data
  //  is cached in coordinate form before creating the actual graph.
  template <std::ranges::input_range EdgeRange>
  compressed_sparse_row_graph(edges_are_unsorted_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), property_(prop) {
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
      : vertex_properties_(numverts), property_(prop) {
    std::vector<vertex_descriptor> sources;
    std::vector<vertex_descriptor> targets;
    histogram_detail::split_into_separate_coords(edge_rg, sources, targets);
    std::size_t numedges = sources.size();
    std::vector<edge_property_type> edge_props(numedges);
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
      : vertex_properties_(numlocalverts), property_(prop) {
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
      : vertex_properties_(numlocalverts), property_(prop) {
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
    vertex_properties_.resize(numverts);
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
    vertex_properties_.resize(numverts);
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
    vertex_properties_.resize(numverts);
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
            typename GlobalToLocal, std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<BidirRange>> &&
      std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EPRange>>
  void add_edges_sorted_internal_global(const BidirRange& sorted_rg, const EPRange& ep_sorted_rg,
                                        const GlobalToLocal& global_to_local, OutputRange& e_out_rg) {
    forward_.add_edges_sorted_internal(sorted_rg, ep_sorted_rg, global_to_local, e_out_rg);
  }

  template <std::ranges::bidirectional_range BidirRange, std::ranges::bidirectional_range EPRange,
            std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<BidirRange>> &&
      std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EPRange>>
  void add_edges_sorted_internal(const BidirRange& sorted_rg, const EPRange& ep_sorted_rg, OutputRange& e_out_rg) {
    forward_.add_edges_sorted_internal(sorted_rg, ep_sorted_rg, identity_property_map(), e_out_rg);
  }

  // Add edges from a sorted (smallest sources first) range of pairs
  template <std::ranges::bidirectional_range BidirRange, std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<BidirRange>>
  void add_edges_sorted_internal(const BidirRange& sorted_rg, OutputRange& e_out_rg) {
    forward_.add_edges_sorted_internal(
        sorted_rg,
        std::ranges::iota_view(std::size_t{0}, static_cast<std::size_t>(std::ranges::distance(sorted_rg))) |
            std::views::transform([](auto /*unused*/) { return edge_property_type{}; }),
        identity_property_map(), e_out_rg);
  }

  template <std::ranges::bidirectional_range BidirRange, typename GlobalToLocal,
            std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<BidirRange>>
  void add_edges_sorted_internal_global(const BidirRange& sorted_rg, const GlobalToLocal& global_to_local,
                                        OutputRange& e_out_rg) {
    forward_.add_edges_sorted_internal(
        sorted_rg,
        std::ranges::iota_view(std::size_t{0}, static_cast<std::size_t>(std::ranges::distance(sorted_rg))) |
            std::views::transform([](auto /*unused*/) { return edge_property_type{}; }),
        global_to_local, e_out_rg);
  }

  // Add edges from a range of (source, target) pairs that are unsorted
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal, std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
  void add_edges_internal_global(const EdgeRange& edge_rg, const GlobalToLocal& global_to_local,
                                 OutputRange& e_out_rg) {
    std::vector<std::pair<std::size_t, std::size_t>> new_edges(edge_rg.begin(), edge_rg.end());
    if (new_edges.empty()) {
      return;
    }
    std::sort(new_edges.begin(), new_edges.end());
    add_edges_sorted_internal_global(new_edges, global_to_local, e_out_rg);
  }

  template <std::ranges::input_range EdgeRange, std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
  void add_edges_internal(const EdgeRange& edge_rg, OutputRange& e_out_rg) {
    add_edges_internal_global(edge_rg, identity_property_map(), e_out_rg);
  }

  // Add edges from a range of (source, target) pairs and edge properties that
  // are unsorted
  template <std::ranges::input_range InRange, std::ranges::input_range EPRange, typename GlobalToLocal,
            std::ranges::bidirectional_range OutputRange>
  void add_edges_internal_global(InRange uv_range, EPRange ep_range, const GlobalToLocal& global_to_local,
                                 OutputRange& e_out_rg) {
    using vertex_pair = std::pair<std::size_t, std::size_t>;
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
    forward_.add_edges_sorted_internal(
        std::views::transform(std::ranges::ref_view(new_edges), [](const auto& x) { return std::get<0>(x); }),
        std::views::transform(std::ranges::ref_view(new_edges), [](const auto& x) { return std::get<1>(x); }),
        global_to_local, e_out_rg);
  }

  // Add edges from a range of (source, target) pairs and edge properties that
  // are unsorted
  template <std::ranges::input_range InRange, std::ranges::input_range EPRange,
            std::ranges::bidirectional_range OutputRange>
  void add_edges_internal(InRange uv_range, EPRange ep_range, OutputRange& e_out_rg) {
    add_edges_internal_global(uv_range, ep_range, identity_property_map(), e_out_rg);
  }

  // Remove edges from a sorted (smallest sources first) range of edges
  template <std::ranges::forward_range ERange, typename GlobalToLocal,
            std::ranges::output_range<EdgeProperty> OutputRange>
  requires std::constructible_from<edge_descriptor, std::ranges::range_value_t<ERange>>
  void remove_edges_sorted_internal_global(const ERange& sorted_rg, const GlobalToLocal& global_to_local,
                                           OutputRange& e_out_rg) {
    forward_.remove_edges_sorted_internal(sorted_rg, global_to_local, e_out_rg);
  }

  template <std::ranges::forward_range ERange, std::ranges::output_range<EdgeProperty> OutputRange>
  requires std::constructible_from<edge_descriptor, std::ranges::range_value_t<ERange>>
  void remove_edges_sorted_internal(const ERange& sorted_rg, OutputRange& e_out_rg) {
    forward_.add_edges_sorted_internal(sorted_rg, identity_property_map(), e_out_rg);
  }

  // Add edges from a range of (source, target) pairs that are unsorted
  template <std::ranges::input_range ERange, typename GlobalToLocal,
            std::ranges::output_range<EdgeProperty> OutputRange>
  requires std::constructible_from<edge_descriptor, std::ranges::range_value_t<ERange>>
  void remove_edges_internal_global(const ERange& edge_rg, const GlobalToLocal& global_to_local,
                                    OutputRange& e_out_rg) {
    std::vector<edge_descriptor> sorted_rg(edge_rg.begin(), edge_rg.end());
    if (sorted_rg.empty()) {
      return;
    }
    std::sort(sorted_rg.begin(), sorted_rg.end());
    remove_edges_sorted_internal_global(sorted_rg, global_to_local, e_out_rg);
  }

  template <std::ranges::input_range ERange, std::ranges::output_range<EdgeProperty> OutputRange>
  requires std::constructible_from<edge_descriptor, std::ranges::range_value_t<ERange>>
  void remove_edges_internal(const ERange& edge_rg, OutputRange& e_out_rg) {
    remove_edges_internal_global(edge_rg, identity_property_map(), e_out_rg);
  }

  // Directly access a vertex or edge bundle
  vertex_bundled& operator[](std::size_t v) { return get_property_value(vertex_properties_[v], vertex_bundle); }
  const vertex_bundled& operator[](std::size_t v) const {
    return get_property_value(vertex_properties_[v], vertex_bundle);
  }
  vertex_property_type& get_property(std::size_t v) { return vertex_properties_[v]; }
  const vertex_property_type& get_property(std::size_t v) const { return vertex_properties_[v]; }

  // Directly access a edge or edge bundle
  edge_bundled& operator[](const edge_descriptor& e) { return forward_[e]; }
  const edge_bundled& operator[](const edge_descriptor& e) const { return forward_[e]; }
  edge_property_type& get_property(const edge_descriptor& e) { return forward_.get_property(e); }
  const edge_property_type& get_property(const edge_descriptor& e) const { return forward_.get_property(e); }

  // Indexing operator. Returns a reference to the graph-bundle associated to the graph.
  graph_bundled& operator[](graph_bundle_t /*unused*/) { return get_property_value(property_, graph_bundle); }
  const graph_bundled& operator[](graph_bundle_t /*unused*/) const {
    return get_property_value(property_, graph_bundle);
  }

  // Get a reference to the graph-property associated to the graph.
  GraphProperty& get_property(graph_all_t /*unused*/) { return property_; }
  const GraphProperty& get_property(graph_all_t /*unused*/) const { return property_; }

  std::vector<VertexProperty> vertex_properties_;
  forward_type forward_;
  GraphProperty property_;
};

template <typename VertexProperty, typename EdgeProperty, typename GraphProperty>
class compressed_sparse_row_graph<bidirectional_s, VertexProperty, EdgeProperty, GraphProperty> {
 public:
  // For Property Graph
  using graph_property_type = GraphProperty;
  using graph_bundled = lookup_one_property_t<GraphProperty, graph_bundle_t>;
  // typedef GraphProperty graph_property_type;

  using forward_type = csr_detail::compressed_sparse_row_structure<EdgeProperty>;
  using backward_edge_property = std::size_t;
  using backward_type = csr_detail::compressed_sparse_row_structure<backward_edge_property>;

  // Concept requirements:
  // For Graph
  using vertex_descriptor = std::size_t;
  using edge_descriptor = csr_detail::csr_edge_descriptor;
  using directed_category = bidirectional_tag;
  using edge_parallel_category = allow_parallel_edge_tag;

  struct traversal_category : bidirectional_graph_tag,
                              adjacency_graph_tag,
                              vertex_list_graph_tag,
                              edge_list_graph_tag {};

  static vertex_descriptor null_vertex() { return std::numeric_limits<vertex_descriptor>::max(); }

  // For VertexListGraph
  using vertices_size_type = std::size_t;

  // For EdgeListGraph
  using edges_size_type = std::size_t;

  // For IncidenceGraph
  using degree_size_type = std::size_t;

  // For internal use
  using graph_tag = csr_graph_tag;

  using edge_bundled = typename forward_type::edge_bundled;
  using edge_property_type = typename forward_type::edge_property_type;

  using vertex_property_type = VertexProperty;
  using vertex_bundled = lookup_one_property_t<VertexProperty, vertex_bundle_t>;

  // Constructors

  // Default constructor: an empty graph.
  compressed_sparse_row_graph() : property_() {}

  //  With numverts vertices
  explicit compressed_sparse_row_graph(std::size_t numverts)
      : vertex_properties_(numverts), forward_(numverts), backward_(numverts) {}

 private:
  void set_up_backward_property_links() {
    auto e_rg = edges(*this);
    auto uv_rg = csr_detail::make_edge_to_index_pair_range(*this, get(vertex_index, *this), e_rg);
    backward_.assign_unsorted_multi_pass_edges(uv_rg, std::views::iota(std::size_t{0}, num_edges(*this)),
                                               forward_.rowstart_.size() - 1,
                                               typed_identity_property_map<std::size_t>(), keep_all());
  }

 public:
  //  From number of vertices and unsorted list of edges
  template <std::ranges::forward_range EdgeRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, std::size_t numverts,
                              const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, numverts, typed_identity_property_map<std::size_t>(),
                                              keep_all());
    set_up_backward_property_links();
  }

  //  From number of vertices and unsorted list of edges, plus edge properties
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropRange>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_t /*unused*/, EdgeRange edge_rg, EdgePropRange ep_rg,
                              std::size_t numverts, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numverts, typed_identity_property_map<std::size_t>(),
                                              keep_all());
    set_up_backward_property_links();
  }

  //  From number of vertices and unsorted list of edges, with filter and
  //  global-to-local map
  template <std::ranges::forward_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  compressed_sparse_row_graph(edges_are_unsorted_multi_pass_global_t /*unused*/, EdgeRange edge_rg,
                              std::size_t numlocalverts, const GlobalToLocal& global_to_local,
                              const SourcePred& source_pred, const GraphProperty& prop = GraphProperty())
      : vertex_properties_(numlocalverts), forward_(), property_(prop) {
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
      : vertex_properties_(numlocalverts), forward_(), property_(prop) {
    forward_.assign_unsorted_multi_pass_edges(edge_rg, ep_rg, numlocalverts, global_to_local, source_pred);
    set_up_backward_property_links();
  }

  //   Requires IncidenceGraph and a vertex index map
  template <typename Graph, typename VertexIndexMap>
  compressed_sparse_row_graph(const Graph& g, const VertexIndexMap& vi, std::size_t numverts, std::size_t numedges)
      : property_() {
    assign(g, vi, numverts, numedges);
    vertex_properties_.resize(numverts);
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
    vertex_properties_.resize(numverts);
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
    vertex_properties_.resize(numverts);
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

  // Directly access a vertex or edge bundle
  vertex_bundled& operator[](std::size_t v) { return get_property_value(vertex_properties_[v], vertex_bundle); }
  const vertex_bundled& operator[](std::size_t v) const {
    return get_property_value(vertex_properties_[v], vertex_bundle);
  }
  vertex_property_type& get_property(std::size_t v) { return vertex_properties_[v]; }
  const vertex_property_type& get_property(std::size_t v) const { return vertex_properties_[v]; }

  // Directly access a edge or edge bundle
  edge_bundled& operator[](const edge_descriptor& e) { return forward_[e]; }
  const edge_bundled& operator[](const edge_descriptor& e) const { return forward_[e]; }
  edge_property_type& get_property(const edge_descriptor& e) { return forward_.get_property(e); }
  const edge_property_type& get_property(const edge_descriptor& e) const { return forward_.get_property(e); }

  // Indexing operator. Returns a reference to the graph-bundle associated to the graph.
  graph_bundled& operator[](graph_bundle_t /*unused*/) { return get_property_value(property_, graph_bundle); }
  const graph_bundled& operator[](graph_bundle_t /*unused*/) const {
    return get_property_value(property_, graph_bundle);
  }

  // Get a reference to the graph-property associated to the graph.
  GraphProperty& get_property(graph_all_t /*unused*/) { return property_; }
  const GraphProperty& get_property(graph_all_t /*unused*/) const { return property_; }

  std::vector<VertexProperty> vertex_properties_;
  forward_type forward_;
  backward_type backward_;
  GraphProperty property_;
};

// Construction functions
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto add_vertex(BAGL_CSR_GRAPH_TYPE& g) {
  return add_vertex(g, typename BAGL_CSR_GRAPH_TYPE::vertex_bundled());
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename... VPArgs>
auto add_vertex(BAGL_DIR_CSR_GRAPH_TYPE& g, VPArgs&&... vp_args) {
  auto old_num_verts_plus_one = g.forward_.rowstart_.size();
  g.forward_.rowstart_.push_back(g.forward_.rowstart_.back());
  g.vertex_properties_.emplace_back(std::forward<VPArgs>(vp_args)...);
  return old_num_verts_plus_one - 1;
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS, typename... VPArgs>
auto add_vertex(BAGL_BIDIR_CSR_GRAPH_TYPE& g, VPArgs&&... vp_args) {
  auto old_num_verts_plus_one = g.forward_.rowstart_.size();
  g.forward_.rowstart_.push_back(g.forward_.rowstart_.back());
  g.backward_.rowstart_.push_back(g.backward_.rowstart_.back());
  g.vertex_properties_.push_back(std::forward<VPArgs>(vp_args)...);
  return old_num_verts_plus_one - 1;
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
auto add_vertices(std::size_t count, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto old_num_verts_plus_one = g.forward_.rowstart_.size();
  std::size_t numedges = g.forward_.rowstart_.back();
  g.forward_.rowstart_.resize(old_num_verts_plus_one + count, numedges);
  g.vertex_properties_.resize(num_vertices(g));
  return old_num_verts_plus_one - 1;
}

// Add edges from a sorted (smallest sources first) range of pairs and edge
// properties
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::bidirectional_range EdgeRange,
          std::ranges::bidirectional_range EdgePropRange>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
    std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EdgePropRange>>
void add_edges_sorted(const EdgeRange& e_range, const EdgePropRange& ep_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_sorted_internal(e_range, ep_range, e_out);
}

// Add edges from a sorted (smallest sources first) range of pairs
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::bidirectional_range EdgeRange>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
void add_edges_sorted(const EdgeRange& e_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_sorted_internal(e_range, e_out);
}

// Remove edges from a sorted (smallest sources first) range of edges
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::forward_range EdgeRange>
requires std::constructible_from<graph_edge_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>,
                                 std::ranges::range_value_t<EdgeRange>>
void remove_edges_sorted(const EdgeRange& e_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.remove_edges_sorted_internal(e_range, e_out);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, typename... EPArgs>
std::pair<typename BAGL_DIR_CSR_GRAPH_TYPE::edge_descriptor, bool> add_edge(
    typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_descriptor u, typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_descriptor v,
    BAGL_DIR_CSR_GRAPH_TYPE& g, EPArgs&&... ep_args) {
  auto result_view = std::ranges::single_view<typename BAGL_DIR_CSR_GRAPH_TYPE::edge_descriptor>();
  g.add_edges_sorted_internal(
      std::ranges::single_view{std::pair{u, v}},
      std::ranges::single_view{typename BAGL_DIR_CSR_GRAPH_TYPE::edge_property_type{std::forward<EPArgs>(ep_args)...}},
      result_view);
  return {*result_view.data(), true};
}

// Remove edge
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
void remove_edge(typename BAGL_DIR_CSR_GRAPH_TYPE::edge_descriptor e, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.remove_edges_sorted_internal(std::ranges::single_view{e}, e_out);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
void remove_edge(typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_descriptor u,
                 typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_descriptor v, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  std::vector<typename BAGL_DIR_CSR_GRAPH_TYPE::edge_descriptor> e_rg;
  for (auto e : out_edges(u, g)) {
    if (target(e, g) == v) {
      e_rg.push_back(e);
    }
  }
  if (e_rg.empty()) {
    return;
  }
  g.remove_edges_sorted_internal(e_rg, e_out);
}

// Clear vertices from a sorted (smallest first) range of vertices
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::forward_range VRange>
requires std::constructible_from<graph_vertex_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>, std::ranges::range_value_t<VRange>>
void clear_vertices_sorted(const VRange& v_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  if (v_range.empty()) {
    return;
  }
  // Gather a sorted list of all edges to remove.
  std::vector<typename BAGL_DIR_CSR_GRAPH_TYPE::edge_descriptor> e_rg;
  auto u = v_range.begin();
  for (auto v : vertices(g)) {
    const bool v_is_u = (u != v_range.end() && v == *u);
    if (v_is_u) {
      ++u;
    }
    for (auto e : out_edges(v, g)) {
      if (v_is_u) {
        e_rg.push_back(e);
      } else {
        if (!std::ranges::equal_range(v_range, target(e, g)).empty()) {
          e_rg.push_back(e);
        }
      }
    }
  }
  if (e_rg.empty()) {
    return;
  }
  auto e_out = container_detail::ignore_output_range();
  g.remove_edges_sorted_internal(e_rg, e_out);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::input_range VRange>
requires std::constructible_from<graph_vertex_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>, std::ranges::range_value_t<VRange>>
void clear_vertices(const VRange& v_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  if (v_range.empty()) {
    return;
  }
  std::vector<graph_vertex_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>> sorted_rg(v_range.begin(), v_range.end());
  std::sort(sorted_rg.begin(), sorted_rg.end());
  clear_vertices_sorted(sorted_rg, g);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
void clear_vertex(typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_descriptor u, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  // Gather a sorted list of all edges to remove.
  std::vector<typename BAGL_DIR_CSR_GRAPH_TYPE::edge_descriptor> e_rg;
  for (auto v : vertices(g)) {
    for (auto e : out_edges(v, g)) {
      if (v == u || target(e, g) == u) {
        e_rg.push_back(e);
      }
    }
  }
  if (e_rg.empty()) {
    return;
  }
  auto e_out = container_detail::ignore_output_range();
  g.remove_edges_sorted_internal(e_rg, e_out);
}

// Remove vertices from a sorted (smallest first) range of vertices
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::forward_range VRange>
requires std::constructible_from<graph_vertex_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>, std::ranges::range_value_t<VRange>>
void remove_vertices_sorted(const VRange& v_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  if (v_range.empty()) {
    return;
  }
  clear_vertices_sorted(v_range, g);
  std::size_t base_i = 0;
  for (auto v_cur = v_range.begin(), v_end = v_range.end(); v_cur != v_end; ++v_cur) {
    auto old_i = *v_cur;
    auto v_next = std::next(v_cur);
    auto old_i_1 = (v_next != v_range.end() ? *v_next : g.vertex_properties_.size());
    std::move(g.vertex_properties_.begin() + old_i + 1, g.vertex_properties_.begin() + old_i_1,
              g.vertex_properties_.begin() + base_i);
    base_i += old_i_1 - old_i;
  }
  g.vertex_properties_.resize(base_i);
  g.forward_.erase_vertices(v_range);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::input_range VRange>
requires std::constructible_from<graph_vertex_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>, std::ranges::range_value_t<VRange>>
void remove_vertices(const VRange& v_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  if (v_range.empty()) {
    return;
  }
  std::vector<graph_vertex_descriptor_t<BAGL_DIR_CSR_GRAPH_TYPE>> sorted_rg(v_range.begin(), v_range.end());
  std::sort(sorted_rg.begin(), sorted_rg.end());
  remove_vertices_sorted(sorted_rg, g);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS>
void remove_vertex(typename BAGL_DIR_CSR_GRAPH_TYPE::vertex_descriptor u, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  if (u >= g.vertex_properties_.size()) {
    return;
  }
  clear_vertex(u, g);
  std::move(g.vertex_properties_.begin() + u + 1, g.vertex_properties_.end(), g.vertex_properties_.begin() + u);
  g.forward_.erase_vertices(std::ranges::single_view{u});
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::bidirectional_range EdgeRange,
          std::ranges::bidirectional_range EdgePropRange, typename GlobalToLocal>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
    std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EdgePropRange>>
void add_edges_sorted_global(const EdgeRange& e_range, const EdgePropRange& ep_range,
                             const GlobalToLocal& global_to_local, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_sorted_internal_global(e_range, ep_range, global_to_local, e_out);
}

// Add edges from a sorted (smallest sources first) range of pairs
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::bidirectional_range EdgeRange, typename GlobalToLocal>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
void add_edges_sorted_global(const EdgeRange& e_range, const GlobalToLocal& global_to_local,
                             BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_sorted_internal_global(e_range, global_to_local, e_out);
}

// Add edges from a range of (source, target) pairs that are unsorted
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::input_range EdgeRange, typename GlobalToLocal>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
void add_edges_global(const EdgeRange& e_range, const GlobalToLocal& global_to_local, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_internal_global(e_range, global_to_local, e_out);
}

// Add edges from a range of (source, target) pairs that are unsorted
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::input_range EdgeRange>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>>
void add_edges(const EdgeRange& e_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_internal(e_range, e_out);
}

// Add edges from a range of (source, target) pairs and edge properties that
// are unsorted
template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
    std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EdgePropRange>>
void add_edges(const EdgeRange& e_range, const EdgePropRange& ep_range, BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_internal(e_range, ep_range, e_out);
}

template <BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS, std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange,
          typename GlobalToLocal>
requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<EdgeRange>> &&
    std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EdgePropRange>>
void add_edges_global(const EdgeRange& e_range, const EdgePropRange& ep_range, const GlobalToLocal& global_to_local,
                      BAGL_DIR_CSR_GRAPH_TYPE& g) {
  auto e_out = container_detail::ignore_output_range();
  g.add_edges_internal_global(e_range, ep_range, global_to_local, e_out);
}

// From VertexListGraph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto num_vertices(const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.rowstart_.size() - 1;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto vertices(const BAGL_CSR_GRAPH_TYPE& g) {
  return std::views::iota(std::size_t{0}, std::size_t{num_vertices(g)});
}

// From IncidenceGraph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto source(typename BAGL_CSR_GRAPH_TYPE::edge_descriptor e, const BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return e.src;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto target(typename BAGL_CSR_GRAPH_TYPE::edge_descriptor e, const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.column_[e.idx];
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto out_edges(std::size_t v, const BAGL_CSR_GRAPH_TYPE& g) {
  return csr_detail::make_csr_out_edge_range<BAGL_CSR_GRAPH_TYPE>(g.forward_.rowstart_[v], g.forward_.rowstart_[v + 1],
                                                                  v);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto out_degree(std::size_t v, const BAGL_CSR_GRAPH_TYPE& g) {
  auto v_row_start = g.forward_.rowstart_[v];
  auto next_row_start = g.forward_.rowstart_[v + 1];
  return next_row_start - v_row_start;
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
auto in_edges(std::size_t v, const BAGL_BIDIR_CSR_GRAPH_TYPE& g) {
  return csr_detail::make_csr_in_edge_range(g.backward_.rowstart_[v], g.backward_.rowstart_[v + 1], g);
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
auto in_degree(std::size_t v, const BAGL_BIDIR_CSR_GRAPH_TYPE& g) {
  auto v_row_start = g.backward_.rowstart_[v];
  auto next_row_start = g.backward_.rowstart_[v + 1];
  return next_row_start - v_row_start;
}

// From AdjacencyGraph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto adjacent_vertices(std::size_t v, const BAGL_CSR_GRAPH_TYPE& g) {
  std::size_t v_row_start = g.forward_.rowstart_[v];
  std::size_t next_row_start = g.forward_.rowstart_[v + 1];
  return std::ranges::subrange(g.forward_.column_.begin() + v_row_start, g.forward_.column_.begin() + next_row_start);
}

// Extra, common functions
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto vertex(std::size_t i, const BAGL_CSR_GRAPH_TYPE& /*unused*/) {
  return i;
}

// edge() can be provided in linear time for the new interface

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
std::pair<typename BAGL_CSR_GRAPH_TYPE::edge_descriptor, bool> edge(std::size_t i, std::size_t j,
                                                                    const BAGL_CSR_GRAPH_TYPE& g) {
  for (auto e : out_edges(i, g)) {
    if (target(e, g) == j) {
      return {e, true};
    }
  }
  return {{}, false};
}

// Find an edge given its index in the graph
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto edge_from_index(std::size_t idx, const BAGL_CSR_GRAPH_TYPE& g) {
  assert(idx < num_edges(g));
  auto src_plus_1 = std::upper_bound(g.forward_.rowstart_.begin(), g.forward_.rowstart_.end(), idx);
  // Get last source whose rowstart is at most idx
  // upper_bound returns this position plus 1
  auto src = (src_plus_1 - g.forward_.rowstart_.begin()) - 1;
  return typename BAGL_CSR_GRAPH_TYPE::edge_descriptor(src, idx);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto num_edges(const BAGL_CSR_GRAPH_TYPE& g) {
  return g.forward_.column_.size();
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto edges(const BAGL_CSR_GRAPH_TYPE& g) {
  return edges_from_out_edges(g);
}

/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

// Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const auto& get(const BAGL_CSR_GRAPH_TYPE& g, typename BAGL_CSR_GRAPH_TYPE::vertex_descriptor v) {
  return g[v];
}
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto& get(BAGL_CSR_GRAPH_TYPE& g, typename BAGL_CSR_GRAPH_TYPE::vertex_descriptor v) {
  return g[v];
}

// Returns a const-reference to the edge-bundle associated to the given edge descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const auto& get(const BAGL_CSR_GRAPH_TYPE& g, const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e) {
  return g[e];
}
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto& get(BAGL_CSR_GRAPH_TYPE& g, const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const auto& get(const BAGL_CSR_GRAPH_TYPE& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto& get(BAGL_CSR_GRAPH_TYPE& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename VProp>
void put(BAGL_CSR_GRAPH_TYPE& g, typename BAGL_CSR_GRAPH_TYPE::vertex_descriptor v, VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-bundle associated to the given edge descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename EProp>
void put(BAGL_CSR_GRAPH_TYPE& g, const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e, EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Sets the graph-bundle associated to the graph.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename GProp>
void put(BAGL_CSR_GRAPH_TYPE& g, graph_bundle_t /*unused*/, GProp&& value) {
  g[graph_bundle] = std::forward<GProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto& get_property(typename BAGL_CSR_GRAPH_TYPE::vertex_descriptor v, BAGL_CSR_GRAPH_TYPE& g) {
  return g.get_property(v);
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const auto& get_property(typename BAGL_CSR_GRAPH_TYPE::vertex_descriptor v, const BAGL_CSR_GRAPH_TYPE& g) {
  return g.get_property(v);
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto& get_property(const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e, BAGL_CSR_GRAPH_TYPE& g) {
  return g.get_property(e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const auto& get_property(const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e, const BAGL_CSR_GRAPH_TYPE& g) {
  return g.get_property(e);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto& get_property(BAGL_CSR_GRAPH_TYPE& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
const auto& get_property(const BAGL_CSR_GRAPH_TYPE& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Handle graph property tags, also handles graph_bundle_t.

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>, lookup_one_property_t<GraphProperty, Tag>&>
get_property(BAGL_CSR_GRAPH_TYPE& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>,
                 const lookup_one_property_t<GraphProperty, Tag>&>
get_property(const BAGL_CSR_GRAPH_TYPE& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Tag, typename T>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>> set_property(BAGL_CSR_GRAPH_TYPE& g,
                                                                                        Tag /*unused*/, T&& value) {
  get_property_value(g.get_property(graph_all), Tag{}) = std::forward<T>(value);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename T, typename Bundle>
struct property_map<BAGL_CSR_GRAPH_TYPE, T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v =
      std::is_convertible_v<typename BAGL_CSR_GRAPH_TYPE::vertex_bundled*, non_const_Bundle*>;
  static constexpr bool is_edge_bundle_v =
      std::is_convertible_v<typename BAGL_CSR_GRAPH_TYPE::edge_bundled*, non_const_Bundle*>;
  using tag_type = std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                                      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<non_const_T, BAGL_CSR_GRAPH_TYPE, tag_type>;
  using const_type = bundle_member_property_map<const non_const_T, const BAGL_CSR_GRAPH_TYPE, tag_type>;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename T, typename Bundle>
auto get(T Bundle::*p, BAGL_CSR_GRAPH_TYPE& g) {
  return typename property_map<BAGL_CSR_GRAPH_TYPE, T Bundle::*>::type{&g, p};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename T, typename Bundle>
auto get(T Bundle::*p, const BAGL_CSR_GRAPH_TYPE& g) {
  return typename property_map<BAGL_CSR_GRAPH_TYPE, T Bundle::*>::const_type{&g, p};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p, const BAGL_CSR_GRAPH_TYPE& g, const Key& k) {
  return (g[k]).*p;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_CSR_GRAPH_TYPE& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, vertex_degree_t> {
  using type = degree_property_map<BAGL_CSR_GRAPH_TYPE>;
  using const_type = degree_property_map<BAGL_CSR_GRAPH_TYPE>;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_degree_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_t<BAGL_CSR_GRAPH_TYPE, vertex_degree_t>{g};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_degree_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_const_t<BAGL_CSR_GRAPH_TYPE, vertex_degree_t>{g};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Key>
auto get(vertex_degree_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g, const Key& k) {
  return degree(k, g);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, vertex_out_degree_t> {
  using type = degree_property_map<BAGL_CSR_GRAPH_TYPE>;
  using const_type = degree_property_map<BAGL_CSR_GRAPH_TYPE>;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_out_degree_t /*unused*/, BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_t<BAGL_CSR_GRAPH_TYPE, vertex_out_degree_t>{g};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_out_degree_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g) {
  return property_map_const_t<BAGL_CSR_GRAPH_TYPE, vertex_out_degree_t>{g};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Key>
auto get(vertex_out_degree_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& g, const Key& k) {
  return out_degree(k, g);
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_BIDIR_CSR_GRAPH_TYPE, vertex_in_degree_t> {
  using type = degree_property_map<BAGL_BIDIR_CSR_GRAPH_TYPE>;
  using const_type = degree_property_map<BAGL_BIDIR_CSR_GRAPH_TYPE>;
};

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_in_degree_t /*unused*/, BAGL_BIDIR_CSR_GRAPH_TYPE& g) {
  return property_map_t<BAGL_BIDIR_CSR_GRAPH_TYPE, vertex_in_degree_t>{g};
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_in_degree_t /*unused*/, const BAGL_BIDIR_CSR_GRAPH_TYPE& g) {
  return property_map_const_t<BAGL_BIDIR_CSR_GRAPH_TYPE, vertex_in_degree_t>{g};
}

template <BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS, typename Key>
auto get(vertex_in_degree_t /*unused*/, const BAGL_BIDIR_CSR_GRAPH_TYPE& g, const Key& k) {
  return in_degree(k, g);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct csr_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;

    using type = tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type = tagged_in_property_property_map<value_type, const Graph, Tag>;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<csr_graph_tag> {
  using type = csr_property_selector;
};

template <>
struct edge_property_selector<csr_graph_tag> {
  using type = csr_property_selector;
};

template <>
struct graph_property_selector<csr_graph_tag> {
  using type = csr_property_selector;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Property>
auto get(Property p, BAGL_CSR_GRAPH_TYPE& g) {
  using Map = property_map_t<BAGL_CSR_GRAPH_TYPE, Property>;
  return Map(&g, p);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Property>
auto get(Property p, const BAGL_CSR_GRAPH_TYPE& g) {
  using Map = property_map_const_t<BAGL_CSR_GRAPH_TYPE, Property>;
  return Map(&g, p);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Property, typename Key>
decltype(auto) get(Property p, const BAGL_CSR_GRAPH_TYPE& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_CSR_GRAPH_TYPE& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS, typename Property, typename Key, typename Value>
void put(Property p, BAGL_CSR_GRAPH_TYPE& g, const Key& k, Value&& val) {
  get_property_value(g.get_property(k), p) = std::forward<Value>(val);
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, vertex_index_t> {
  using type = identity_property_map;
  using const_type = type;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*g*/) {
  return identity_property_map{};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(vertex_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*g*/) {
  return identity_property_map{};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
std::size_t get(vertex_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*g*/, const std::size_t& k) {
  return k;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
std::size_t get(vertex_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*g*/, const std::size_t& k) {
  return k;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
struct property_map<BAGL_CSR_GRAPH_TYPE, edge_index_t> {
  using type = csr_detail::csr_edge_index_map;
  using const_type = type;
};

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(edge_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*g*/) {
  return csr_detail::csr_edge_index_map{};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
auto get(edge_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*g*/) {
  return csr_detail::csr_edge_index_map{};
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
std::size_t get(edge_index_t /*unused*/, const BAGL_CSR_GRAPH_TYPE& /*g*/,
                const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e) {
  return e.idx;
}

template <BAGL_CSR_GRAPH_TEMPLATE_PARMS>
std::size_t get(edge_index_t /*unused*/, BAGL_CSR_GRAPH_TYPE& /*g*/,
                const typename BAGL_CSR_GRAPH_TYPE::edge_descriptor& e) {
  return e.idx;
}

#undef BAGL_CSR_GRAPH_TYPE
#undef BAGL_CSR_GRAPH_TEMPLATE_PARMS
#undef BAGL_DIR_CSR_GRAPH_TYPE
#undef BAGL_DIR_CSR_GRAPH_TEMPLATE_PARMS
#undef BAGL_BIDIR_CSR_GRAPH_TYPE
#undef BAGL_BIDIR_CSR_GRAPH_TEMPLATE_PARMS

}  // namespace bagl

#endif  // BAGL_BAGL_COMPRESSED_SPARSE_ROW_GRAPH_H_
