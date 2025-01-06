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
#include "bagl/filtered_graph.h"  // For keep_all
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl::csr_detail {

struct csr_edge_descriptor {
  std::size_t src = 0;
  std::size_t idx = 0;

  csr_edge_descriptor(std::size_t a_src, std::size_t a_idx) : src(a_src), idx(a_idx) {}
  csr_edge_descriptor() = default;

  bool operator==(const csr_edge_descriptor& e) const { return idx == e.idx; }
  bool operator!=(const csr_edge_descriptor& e) const { return idx != e.idx; }
  auto operator<=>(const csr_edge_descriptor& e) const { return idx <=> e.idx; }
};

// Add edge_index property map
struct csr_edge_index_map {
  using value_type = std::size_t;
};

inline std::size_t get(const csr_edge_index_map&, const csr_edge_descriptor& key) { return key.idx; }

/** Compressed sparse row graph internal structure.
 */
template <typename EdgeProperty>
class compressed_sparse_row_structure {
 public:
  using vertex_descriptor = std::size_t;

  static vertex_descriptor null_vertex() { return std::numeric_limits<vertex_descriptor>::max(); }

  std::vector<std::size_t> rowstart_;
  std::vector<std::size_t> column_;
  std::vector<EdgeProperty> edge_properties_;

  using edge_property_type = EdgeProperty;
  using edge_bundled = lookup_one_property_t<EdgeProperty, edge_bundle_t>;

  // Directly access an edge property.
  edge_property_type& get_property(const csr_edge_descriptor& e) { return edge_properties_[e.idx]; }
  const edge_property_type& get_property(const csr_edge_descriptor& e) const { return edge_properties_[e.idx]; }

  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  edge_bundled& operator[](const csr_edge_descriptor& e) {
    return get_property_value(edge_properties_[e.idx], edge_bundle);
  }
  const edge_bundled& operator[](const csr_edge_descriptor& e) const {
    return get_property_value(edge_properties_[e.idx], edge_bundle);
  }

  explicit compressed_sparse_row_structure(std::size_t numverts = 0)
      : rowstart_(numverts + 1, std::size_t(0)), column_() {}

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
    edge_properties_.resize(rowstart_.back());

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
    edge_properties_.resize(rowstart_.back());

    histogram_detail::histogram_sort(sources, std::ranges::ref_view(rowstart_), numlocalverts, targets,
                                     std::ranges::ref_view(column_), ep_rg, std::ranges::ref_view(edge_properties_),
                                     source_pred, make_property_map_function(global_to_local));
  }

  //  Assign from number of vertices and sorted list of edges
  template <std::ranges::input_range EdgeRange, typename GlobalToLocal, typename SourcePred>
  void assign_from_sorted_edges(EdgeRange edges_rg, const GlobalToLocal& global_to_local, const SourcePred& source_pred,
                                std::size_t numlocalverts, std::size_t numedges_or_zero) {
    column_.clear();
    column_.reserve(numedges_or_zero);
    rowstart_.resize(numlocalverts + 1);
    std::size_t current_edge = 0;
    std::size_t current_vertex_plus_one = 1;
    rowstart_[0] = 0;
    for (auto vp : edges_rg) {
      if (!source_pred(vp.first)) {
        continue;
      }
      std::size_t src = get(global_to_local, vp.first);
      std::size_t tgt = vp.second;
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
    edge_properties_.resize(column_.size());
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
    edge_properties_.clear();
    edge_properties_.reserve(numedges_or_zero);
    rowstart_.resize(numlocalverts + 1);
    std::size_t current_edge = 0;
    std::size_t current_vertex_plus_one = 1;
    rowstart_[0] = 0;
    for (const auto& [vp, ep] : zip_range(edges_rg, ep_rg)) {
      if (!source_pred(vp.first)) {
        continue;
      }
      std::size_t src = get(global_to_local, vp.first);
      std::size_t tgt = vp.second;
      for (; current_vertex_plus_one != src + 1; ++current_vertex_plus_one) {
        rowstart_[current_vertex_plus_one] = current_edge;
      }
      column_.push_back(tgt);
      edge_properties_.push_back(ep);
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
    edge_properties_.resize(rowstart_.back());
  }

  // Replace graph with sources and targets and edge properties given,
  // sorting them in-place, and using the given global-to-local property
  // map to get local indices from global ones in the two arrays.
  template <typename GlobalToLocal>
  void assign_sources_and_targets_global(std::vector<vertex_descriptor>& sources,
                                         std::vector<vertex_descriptor>& targets, std::vector<EdgeProperty>& edge_props,
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
    edge_properties_.resize(numedges);
    std::size_t current_edge = 0;
    using g_vertex = graph_vertex_descriptor_t<Graph>;

    std::vector<g_vertex> ordered_verts_of_g(numverts);
    for (auto v : vertices(g)) {
      ordered_verts_of_g[get(vertex_index, g, v)] = v;
    }
    for (std::size_t i = 0; i != numverts; ++i) {
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
            typename GlobalToLocal, std::ranges::bidirectional_range OutputRange>
  requires std::constructible_from<std::pair<std::size_t, std::size_t>, std::ranges::range_value_t<BidirRange>> &&
      std::constructible_from<EdgeProperty, std::ranges::range_reference_t<EPRange>> &&
      std::ranges::output_range<OutputRange, csr_edge_descriptor>
  void add_edges_sorted_internal(const BidirRange& sorted_rg, const EPRange& ep_sorted_rg,
                                 const GlobalToLocal& global_to_local, OutputRange& e_out_rg) {
    // Flip sequence
    auto reversed_rg = std::ranges::ref_view(sorted_rg) | std::views::reverse;
    std::size_t new_edge_count = std::ranges::distance(reversed_rg);
    auto ep_reversed_rg = std::ranges::ref_view(ep_sorted_rg) | std::views::reverse;
    auto e_out_reversed_rg = std::ranges::ref_view(e_out_rg) | std::views::reverse;
    std::size_t edges_added_before_i = new_edge_count;  // Count increment to add to rowstarts
    column_.resize(column_.size() + new_edge_count);
    edge_properties_.resize(edge_properties_.size() + new_edge_count);
    auto current_new_edge = reversed_rg.begin();
    auto prev_new_edge = reversed_rg.begin();
    auto current_new_edge_prop = ep_reversed_rg.begin();
    auto current_new_e_out = e_out_reversed_rg.begin();
    for (std::size_t i_plus_1 = rowstart_.size() - 1; i_plus_1 > 0; --i_plus_1) {
      std::size_t i = i_plus_1 - 1;
      prev_new_edge = current_new_edge;
      // edges_added_to_this_vertex = #mbrs of new_edges with reversed_rg.begin() == i
      std::size_t edges_added_to_this_vertex = 0;
      while (current_new_edge != reversed_rg.end()) {
        if (get(global_to_local, (*current_new_edge).first) != i) {
          break;
        }
        ++current_new_edge;
        ++current_new_edge_prop;
        ++current_new_e_out;
        ++edges_added_to_this_vertex;
      }
      edges_added_before_i -= edges_added_to_this_vertex;
      // Invariant: edges_added_before_i = #mbrs of new_edges with first < i
      std::size_t old_rowstart = rowstart_[i];
      std::size_t new_rowstart = rowstart_[i] + edges_added_before_i;
      std::size_t old_degree = rowstart_[i + 1] - rowstart_[i];
      std::size_t new_degree = old_degree + edges_added_to_this_vertex;
      // Move old edges forward (by #new_edges before this i) to make
      // room new_rowstart > old_rowstart, so use move_backward
      if (old_rowstart != new_rowstart) {
        std::move_backward(column_.begin() + old_rowstart, column_.begin() + old_rowstart + old_degree,
                           column_.begin() + new_rowstart + old_degree);
        std::move_backward(edge_properties_.begin() + old_rowstart,
                           edge_properties_.begin() + old_rowstart + old_degree,
                           edge_properties_.begin() + new_rowstart + old_degree);
      }
      // Add new edges (reversed because current_new_edge is a const_reverse_iterator)
      auto temp = current_new_edge;
      auto temp_prop = current_new_edge_prop;
      auto temp_e_out = current_new_e_out;
      for (; temp != prev_new_edge; ++old_degree) {
        --temp;
        --temp_prop;
        --temp_e_out;
        column_[new_rowstart + old_degree] = (*temp).second;
        edge_properties_[new_rowstart + old_degree] = *temp_prop;
        *temp_e_out = csr_edge_descriptor(i, new_rowstart + old_degree);
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

  template <std::ranges::forward_range ERange, typename GlobalToLocal,
            std::ranges::output_range<EdgeProperty> OutputRange>
  requires std::constructible_from<csr_edge_descriptor, std::ranges::range_value_t<ERange>>
  void remove_edges_sorted_internal(const ERange& sorted_rg, const GlobalToLocal& global_to_local,
                                    OutputRange& ep_out_rg) {
    auto current_e_rm = sorted_rg.begin();
    auto current_ep_out = ep_out_rg.begin();
    std::size_t old_rowstart = rowstart_[0];
    std::size_t new_rowstart = old_rowstart;
    for (std::size_t i = 0; i < rowstart_.size() - 1; ++i) {
      // Move old edges back (by #removed_edges so far) to compress
      // when new_rowstart < old_rowstart, so use move
      std::size_t old_degree = rowstart_[i + 1] - rowstart_[i];
      if (old_rowstart != new_rowstart) {
        std::move(column_.begin() + old_rowstart, column_.begin() + old_rowstart + old_degree,
                  column_.begin() + new_rowstart);
        std::move(edge_properties_.begin() + old_rowstart, edge_properties_.begin() + old_rowstart + old_degree,
                  edge_properties_.begin() + new_rowstart);
        rowstart_[i] = new_rowstart;
      }
      // Find next e_rm that doesn't have 'i' as source.
      auto next_e_rm = current_e_rm;
      while (next_e_rm != sorted_rg.end()) {
        csr_edge_descriptor e = *next_e_rm;
        if (get(global_to_local, e.src) != i) {
          break;
        }
        ++next_e_rm;
      }
      // A simple "remove" loop.
      std::size_t old_outedge = new_rowstart;
      std::size_t new_outedge = old_outedge;
      while (old_outedge != new_rowstart + old_degree) {
        if (current_e_rm == next_e_rm) {
          // Skip to the end of remaining out edges.
          new_outedge += (new_rowstart + old_degree) - old_outedge;
          old_outedge = new_rowstart + old_degree;
          break;
        }
        csr_edge_descriptor e = *current_e_rm;
        std::size_t old_outedge_idx = (old_outedge - new_rowstart) + old_rowstart;
        if (e.idx == old_outedge_idx) {
          // Remove this edge.
          *current_ep_out = std::move(edge_properties_[new_outedge]);
          ++old_outedge;
          ++current_e_rm;
          ++current_ep_out;
        } else {
          // Keep this edge.
          if (old_outedge != new_outedge) {
            column_[new_outedge] = column_[old_outedge];
            edge_properties_[new_outedge] = std::move(edge_properties_[old_outedge]);
          }
          ++old_outedge;
          ++new_outedge;
        }
      }
      // Now, new_outedge is at the end of this row.
      // If we had not-found edges in current e_rm range, just skip them.
      current_e_rm = next_e_rm;
      old_rowstart = rowstart_[i + 1];
      new_rowstart = new_outedge;
      // If current_e_rm reached the end, we continue to move edges and update row starts.
    }
    // We should end up with new_rowstart pointing to the last "rowstart".
    rowstart_.back() = new_rowstart;
    column_.resize(new_rowstart);
    edge_properties_.resize(new_rowstart);
  }

  // Pre: The vertices must be cleared already (no in or out edges).
  template <std::ranges::forward_range VRange>
  void erase_vertices(const VRange& sorted_rg) {
    // If a vertex is cleared, rowstart_[u] == rowstart_[u+1], and u should not appear in any column_ entry.
    // We just need to shift down all vertex indices greater than u.
    for (auto& v : column_) {
      std::size_t id_shift = 0;
      for (auto u : sorted_rg) {
        if (v <= u) {
          break;
        }
        ++id_shift;
      }
      v -= id_shift;
    }
  }
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

template <>
struct std::hash<bagl::csr_detail::csr_edge_descriptor> {
  std::size_t operator()(bagl::csr_detail::csr_edge_descriptor const& x) const {
    return std::hash<std::size_t>()(x.src) ^ std::hash<std::size_t>()(x.idx);
  }
};

#endif  // BAGL_BAGL_DETAIL_COMPRESSED_SPARSE_ROW_STRUCT_H_
