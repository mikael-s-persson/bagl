// Copyright 2005-2009 The Trustees of Indiana University.
// Authors: Jeremiah Willcock, Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

// Compressed sparse row graph type internal structure

#ifndef BAGL_BAGL_DETAIL_COMPRESSED_SPARSE_ROW_STRUCT_H_
#define BAGL_BAGL_DETAIL_COMPRESSED_SPARSE_ROW_STRUCT_H_

#include <limits>
#ifndef BAGL_BAGL_COMPRESSED_SPARSE_ROW_GRAPH_H_
#error This file should only be included from bagl/compressed_sparse_row_graph.h
#endif

#include <algorithm>
#include <cassert>
#include <climits>
#include <iterator>
#include <ranges>
#include <tuple>
#include <utility>
#include <vector>

#include "bagl/detail/histogram_sort.h"
#include "bagl/detail/indexed_properties.h"
#include "bagl/filtered_graph.h"  // For keep_all
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl::csr_detail {
// Forward declaration of CSR edge descriptor type, needed to pass to
// indexed_edge_properties.
template <typename Vertex, typename EdgeIndex>
struct csr_edge_descriptor;

// Add edge_index property map
template <typename Vertex, typename EdgeIndex>
struct csr_edge_index_map {
  using value_type = EdgeIndex;
  using reference = EdgeIndex;
  using key_type = csr_edge_descriptor<Vertex, EdgeIndex>;
  using category = readable_property_map_tag;
};

template <typename Vertex, typename EdgeIndex>
EdgeIndex get(const csr_edge_index_map<Vertex, EdgeIndex>&, const csr_edge_descriptor<Vertex, EdgeIndex>& key) {
  return key.idx;
}

/** Compressed sparse row graph internal structure.
 *
 * Vertex and EdgeIndex should be unsigned integral types and should
 * specialize numeric_limits.
 */
template <typename EdgeProperty, typename Vertex = std::size_t, typename EdgeIndex = Vertex>
class compressed_sparse_row_structure
    : public indexed_detail::indexed_edge_properties<compressed_sparse_row_structure<EdgeProperty, Vertex, EdgeIndex>,
                                                     EdgeProperty, csr_edge_descriptor<Vertex, EdgeIndex>,
                                                     csr_edge_index_map<Vertex, EdgeIndex>> {
 public:
  using inherited_edge_properties =
      indexed_detail::indexed_edge_properties<compressed_sparse_row_structure<EdgeProperty, Vertex, EdgeIndex>,
                                              EdgeProperty, csr_edge_descriptor<Vertex, EdgeIndex>,
                                              csr_edge_index_map<Vertex, EdgeIndex>>;

  using vertices_size_type = Vertex;
  using vertex_descriptor = Vertex;
  using edges_size_type = EdgeIndex;

  static vertex_descriptor null_vertex() { return std::numeric_limits<vertex_descriptor>::max(); }

  std::vector<EdgeIndex> rowstart_;
  std::vector<Vertex> column_;

  explicit compressed_sparse_row_structure(Vertex numverts = 0) : rowstart_(numverts + 1, EdgeIndex(0)), column_() {}

  //  Rebuild graph from number of vertices and multi-pass unsorted list
  //  of edges (filtered using source_pred and mapped using
  //  global_to_local)
  template <std::ranges::forward_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  void assign_unsorted_multi_pass_edges(EdgeRange edges_rg, std::size_t numlocalverts,
                                        const GlobalToLocal& global_to_local, const SourcePred& source_pred) {
    rowstart_.clear();
    rowstart_.resize(numlocalverts + 1, 0);
    auto sources = edges_rg | std::views::transform([](const auto& p) { return p.first; });
    auto targets = edges_rg | std::views::transform([](const auto& p) { return p.second; });

    histogram_detail::count_starts(sources, std::ranges::ref_view(rowstart_), numlocalverts, source_pred,
                                   make_property_map_function(global_to_local));

    column_.resize(rowstart_.back());
    inherited_edge_properties::resize(rowstart_.back());

    histogram_detail::histogram_sort(sources, std::ranges::ref_view(rowstart_), numlocalverts, targets,
                                     std::ranges::ref_view(column_), source_pred,
                                     make_property_map_function(global_to_local));
  }

  //  Rebuild graph from number of vertices and multi-pass unsorted list
  //  of edges and their properties (filtered using source_pred and mapped
  //  using global_to_local)
  template <std::ranges::forward_range EdgeRange, std::ranges::input_range EdgePropertyRange, typename GlobalToLocal,
            typename SourcePred>
  void assign_unsorted_multi_pass_edges(EdgeRange edges_rg, EdgePropertyRange ep_rg, std::size_t numlocalverts,
                                        const GlobalToLocal& global_to_local, const SourcePred& source_pred) {
    rowstart_.clear();
    rowstart_.resize(numlocalverts + 1, 0);
    auto sources = edges_rg | std::views::transform([](const auto& p) { return p.first; });
    auto targets = edges_rg | std::views::transform([](const auto& p) { return p.second; });

    histogram_detail::count_starts(sources, std::ranges::ref_view(rowstart_), numlocalverts, source_pred,
                                   make_property_map_function(global_to_local));

    column_.resize(rowstart_.back());
    inherited_edge_properties::resize(rowstart_.back());

    histogram_detail::histogram_sort(sources, std::ranges::ref_view(rowstart_), numlocalverts, targets,
                                     std::ranges::ref_view(column_), ep_rg, inherited_edge_properties::as_range(),
                                     source_pred, make_property_map_function(global_to_local));
  }

  //  Assign from number of vertices and sorted list of edges
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  void assign_from_sorted_edges(EdgeRange edges_rg, const GlobalToLocal& global_to_local, const SourcePred& source_pred,
                                std::size_t numlocalverts, std::size_t numedges_or_zero) {
    column_.clear();
    column_.reserve(numedges_or_zero);
    rowstart_.resize(numlocalverts + 1);
    EdgeIndex current_edge = 0;
    Vertex current_vertex_plus_one = 1;
    rowstart_[0] = 0;
    for (auto vp : edges_rg) {
      if (!source_pred(vp.first)) {
        continue;
      }
      Vertex src = get(global_to_local, vp.first);
      Vertex tgt = vp.second;
      for (; current_vertex_plus_one != src + 1; ++current_vertex_plus_one) {
        rowstart_[current_vertex_plus_one] = current_edge;
      }
      column_.push_back(tgt);
      ++current_edge;
    }

    // The remaining vertices have no edges
    for (; current_vertex_plus_one != numlocalverts + 1; ++current_vertex_plus_one) {
      rowstart_[current_vertex_plus_one] = current_edge;
    }

    // Default-construct properties for edges
    inherited_edge_properties::resize(column_.size());
  }

  //  Assign from number of vertices and sorted list of edges
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropertyRange, typename GlobalToLocal,
            typename SourcePred>
  void assign_from_sorted_edges(EdgeRange edges_rg, EdgePropertyRange ep_rg, const GlobalToLocal& global_to_local,
                                const SourcePred& source_pred, std::size_t numlocalverts,
                                std::size_t numedges_or_zero) {
    // Reserving storage in advance can save us lots of time and
    // memory, but it can only be done if we have forward iterators or
    // the user has supplied the number of edges.
    std::size_t numedges = numedges_or_zero;
    if (numedges == 0) {
      numedges = histogram_detail::reserve_count_for_single_pass(edges_rg);
    }
    column_.clear();
    column_.reserve(numedges_or_zero);
    inherited_edge_properties::clear();
    inherited_edge_properties::reserve(numedges_or_zero);
    rowstart_.resize(numlocalverts + 1);
    EdgeIndex current_edge = 0;
    Vertex current_vertex_plus_one = 1;
    rowstart_[0] = 0;
    for (const auto& [vp, ep] : zip_range(edges_rg, ep_rg)) {
      if (!source_pred(vp.first)) {
        continue;
      }
      Vertex src = get(global_to_local, vp.first);
      Vertex tgt = vp.second;
      for (; current_vertex_plus_one != src + 1; ++current_vertex_plus_one) {
        rowstart_[current_vertex_plus_one] = current_edge;
      }
      column_.push_back(tgt);
      inherited_edge_properties::push_back(ep);
      ++current_edge;
    }

    // The remaining vertices have no edges
    for (; current_vertex_plus_one != numlocalverts + 1; ++current_vertex_plus_one) {
      rowstart_[current_vertex_plus_one] = current_edge;
    }
  }

  // Replace graph with sources and targets given, sorting them in-place,
  // and using the given global-to-local property map to get local indices
  // from global ones in the two arrays.
  template <typename GlobalToLocal>
  void assign_sources_and_targets_global(std::vector<vertex_descriptor>& sources,
                                         std::vector<vertex_descriptor>& targets, std::size_t numverts,
                                         GlobalToLocal global_to_local) {
    assert(sources.size() == targets.size());
    // Do an in-place histogram sort (at least that's what I think it
    // is) to sort sources and targets
    rowstart_.clear();
    rowstart_.resize(numverts + 1);
    histogram_detail::count_starts(std::ranges::ref_view(sources), std::ranges::ref_view(rowstart_), numverts,
                                   keep_all(), make_property_map_function(global_to_local));
    histogram_detail::histogram_sort_inplace(std::ranges::ref_view(sources), std::ranges::ref_view(rowstart_), numverts,
                                             std::ranges::ref_view(targets),
                                             make_property_map_function(global_to_local));
    // Now targets is the correct vector (properly sorted by source) for
    // column_
    column_.swap(targets);
    inherited_edge_properties::resize(rowstart_.back());
  }

  // Replace graph with sources and targets and edge properties given,
  // sorting them in-place, and using the given global-to-local property
  // map to get local indices from global ones in the two arrays.
  template <typename GlobalToLocal>
  void assign_sources_and_targets_global(std::vector<vertex_descriptor>& sources,
                                         std::vector<vertex_descriptor>& targets,
                                         std::vector<typename inherited_edge_properties::edge_bundled>& edge_props,
                                         std::size_t numverts, GlobalToLocal global_to_local) {
    assert(sources.size() == targets.size());
    assert(sources.size() == edge_props.size());
    // Do an in-place histogram sort (at least that's what I think it
    // is) to sort sources and targets
    rowstart_.clear();
    rowstart_.resize(numverts + 1);
    histogram_detail::count_starts(std::ranges::ref_view(sources), std::ranges::ref_view(rowstart_), numverts,
                                   keep_all(), make_property_map_function(global_to_local));
    histogram_detail::histogram_sort_inplace(std::ranges::ref_view(sources), std::ranges::ref_view(rowstart_), numverts,
                                             std::ranges::ref_view(targets), std::ranges::ref_view(edge_props),
                                             make_property_map_function(global_to_local));
    // Now targets is the correct vector (properly sorted by source) for
    // column_, and edge_props for edge_properties_
    column_.swap(targets);
    this->edge_properties_.swap(edge_props);
  }

  // From any graph (slow and uses a lot of memory)
  //   Requires IncidenceGraph and a vertex index map
  //   Internal helper function
  //   Note that numedges must be doubled for undirected source graphs
  template <typename Graph, typename VertexIndexMap>
  void assign(const Graph& g, const VertexIndexMap& vi, std::size_t numverts, std::size_t numedges) {
    rowstart_.resize(numverts + 1);
    column_.resize(numedges);
    inherited_edge_properties::resize(numedges);
    EdgeIndex current_edge = 0;
    using g_vertex = graph_vertex_descriptor_t<Graph>;

    std::vector<g_vertex> ordered_verts_of_g(numverts);
    for (auto v : vertices(g)) {
      ordered_verts_of_g[get(vertex_index, g, v)] = v;
    }
    for (Vertex i = 0; i != numverts; ++i) {
      rowstart_[i] = current_edge;
      g_vertex v = ordered_verts_of_g[i];
      for (auto e : out_edges(v, g)) {
        column_[current_edge++] = get(vi, target(e, g));
      }
    }
    rowstart_[numverts] = current_edge;
  }

  // Add edges from a sorted (smallest sources first) range of pairs and
  // edge properties
  template <std::ranges::bidirectional_range BidirRange, std::ranges::bidirectional_range EPRange,
            typename GlobalToLocal>
  void add_edges_sorted_internal(BidirRange sorted_rg, EPRange ep_sorted_rg, const GlobalToLocal& global_to_local) {
    // Flip sequence
    auto reversed_rg = sorted_rg | std::views::reverse;
    EdgeIndex new_edge_count = std::ranges::distance(reversed_rg);
    auto ep_reversed_rg = ep_sorted_rg | std::views::reverse;
    EdgeIndex edges_added_before_i = new_edge_count;  // Count increment to add to rowstarts
    column_.resize(column_.size() + new_edge_count);
    inherited_edge_properties::resize(inherited_edge_properties::size() + new_edge_count);
    auto current_new_edge = reversed_rg.begin();
    auto prev_new_edge = reversed_rg.begin();
    auto current_new_edge_prop = ep_reversed_rg.begin();
    for (Vertex i_plus_1 = rowstart_.size() - 1; i_plus_1 > 0; --i_plus_1) {
      Vertex i = i_plus_1 - 1;
      prev_new_edge = current_new_edge;
      // edges_added_to_this_vertex = #mbrs of new_edges with reversed_rg.begin() == i
      EdgeIndex edges_added_to_this_vertex = 0;
      while (current_new_edge != reversed_rg.end()) {
        if (get(global_to_local, current_new_edge->first) != i) {
          break;
        }
        ++current_new_edge;
        ++current_new_edge_prop;
        ++edges_added_to_this_vertex;
      }
      edges_added_before_i -= edges_added_to_this_vertex;
      // Invariant: edges_added_before_i = #mbrs of new_edges with first < i
      EdgeIndex old_rowstart = rowstart_[i];
      EdgeIndex new_rowstart = rowstart_[i] + edges_added_before_i;
      EdgeIndex old_degree = rowstart_[i + 1] - rowstart_[i];
      EdgeIndex new_degree = old_degree + edges_added_to_this_vertex;
      // Move old edges forward (by #new_edges before this i) to make
      // room new_rowstart > old_rowstart, so use copy_backwards
      if (old_rowstart != new_rowstart) {
        std::copy_backward(column_.begin() + old_rowstart, column_.begin() + old_rowstart + old_degree,
                           column_.begin() + new_rowstart + old_degree);
        inherited_edge_properties::move_range(old_rowstart, old_rowstart + old_degree, new_rowstart);
      }
      // Add new edges (reversed because current_new_edge is a const_reverse_iterator)
      auto temp = current_new_edge;
      auto temp_prop = current_new_edge_prop;
      for (; temp != prev_new_edge; ++old_degree) {
        --temp;
        --temp_prop;
        column_[new_rowstart + old_degree] = temp->second;
        inherited_edge_properties::write_by_index(new_rowstart + old_degree, *temp_prop);
      }
      rowstart_[i + 1] = new_rowstart + new_degree;
      if (edges_added_before_i == 0) {
        // No more edges inserted before this point
        // rowstart_[i] will be fixed up on the next iteration (to
        // avoid changing the degree of vertex i - 1); the last
        // iteration never changes it (either because of the condition
        // of the break or because rowstart_[0] is always 0)
        break;
      }
    }
  }
};

template <typename Vertex, typename EdgeIndex>
struct csr_edge_descriptor {
  Vertex src{};
  EdgeIndex idx{};

  csr_edge_descriptor(Vertex a_src, EdgeIndex a_idx) : src(a_src), idx(a_idx) {}
  csr_edge_descriptor() = default;

  auto operator<=>(const csr_edge_descriptor& e) const { return idx <=> e.idx; }
};

// Common out edge and edge iterators
template <typename CSRGraph>
auto make_csr_out_edge_range(int f, int l, typename CSRGraph::vertex_descriptor u) {
  return std::views::iota(f, l) |
         std::views::transform([u](int i) { return typename CSRGraph::edge_descriptor(u, i); });
}

// Only for bidirectional graphs
template <typename CSRGraph>
auto make_csr_in_edge_range(std::size_t f, std::size_t l, const CSRGraph& g) {
  return std::views::iota(f, l) | std::views::transform([&g](std::size_t i) {
           return typename CSRGraph::edge_descriptor(g.backward_.column_[i], g.backward_.edge_properties_[i]);
         });
}

template <typename Range>
auto transpose_edges(Range&& r) {
  return std::ranges::transform_view(std::forward<Range>(r), [](const auto& p) {
    return std::pair{p.second, p.first};
  });
}

template <typename GraphT, typename VertexIndexMap>
class edge_to_index_pair {
 public:
  edge_to_index_pair() = default;
  edge_to_index_pair(const GraphT& g, const VertexIndexMap& index) : g_(&g), index_(index) {}

  auto operator()(graph_edge_descriptor_t<GraphT> e) const {
    return std::pair{get(index_, source(e, *g_)), get(index_, target(e, *g_))};
  }

 private:
  const GraphT* g_ = nullptr;
  VertexIndexMap index_;
};

template <typename GraphT, typename VertexIndexMap>
auto make_edge_to_index_pair(const GraphT& g, const VertexIndexMap& index) {
  return edge_to_index_pair(g, index);
}

template <typename GraphT>
auto make_edge_to_index_pair(const GraphT& g) {
  return edge_to_index_pair(g, get(vertex_index, g));
}

template <typename GraphT, typename VertexIndexMap, typename Range>
auto make_edge_to_index_pair_range(const GraphT& g, const VertexIndexMap& index, Range&& r) {
  return std::ranges::transform_view(std::forward<Range>(r), edge_to_index_pair(g, index));
}

}  // namespace bagl::csr_detail

template <typename Vertex, typename EdgeIndex>
struct std::hash<bagl::csr_detail::csr_edge_descriptor<Vertex, EdgeIndex>> {
  std::size_t operator()(bagl::csr_detail::csr_edge_descriptor<Vertex, EdgeIndex> const& x) const {
    return std::hash<Vertex>()(x.src) ^ std::hash<EdgeIndex>()(x.idx);
  }
};

#endif  // BAGL_BAGL_DETAIL_COMPRESSED_SPARSE_ROW_STRUCT_H_
