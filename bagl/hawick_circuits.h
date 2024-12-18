// Copyright Louis Dionne 2013
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_HAWICK_CIRCUITS_H_
#define BAGL_BAGL_HAWICK_CIRCUITS_H_

#include <algorithm>
#include <cassert>
#include <ranges>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/one_bit_color_map.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/tiernan_all_cycles.h"

namespace bagl {
namespace hawick_circuits_detail {
//! @internal Functor returning all the vertices adjacent to a vertex.
struct get_all_adjacent_vertices {
  template <typename Vertex, typename Graph>
  auto operator()(Vertex&& v, Graph&& g) const {
    return adjacent_vertices(std::forward<Vertex>(v), std::forward<Graph>(g));
  }
};

//! @internal Functor returning a set of the vertices adjacent to a vertex.
struct get_unique_adjacent_vertices {
  template <typename Vertex, typename Graph>
  auto operator()(Vertex v, const Graph& g) const {
    auto adj_rg = adjacent_vertices(v, g);
    std::vector<Vertex> uniq_verts;
    for (Vertex u : adj_rg) {
      if (std::find(uniq_verts.begin(), uniq_verts.end(), u) == uniq_verts.end()) {
        uniq_verts.push_back(u);
      }
    }
    return uniq_verts;
  }
};

//! @internal
//! Return whether a container contains a given value.
//! This is not meant as a general purpose membership testing function; it
//! would have to be more clever about possible optimizations.
template <std::ranges::input_range InRange, typename Value>
bool contains(const InRange& c, const Value& v) {
  return std::ranges::find(c, v) != c.end();
}

/*!
 * @internal
 * Algorithm finding all the cycles starting from a given vertex.
 *
 * The search is only done in the subgraph induced by the starting vertex
 * and the vertices with an index higher than the starting vertex.
 */
template <concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap, typename Stack, typename ClosedMatrix,
          typename GetAdjacentVertices>
requires concepts::AdjacencyGraph<Graph>
struct hawick_circuits_from {
 private:
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // The one_bit_color_map starts all white, i.e. not blocked.
  // Since we make that assumption (I looked at the implementation, but
  // I can't find anything that documents this behavior), we're gonna
  // assert it in the constructor.
  using BlockedMap = one_bit_color_map<VertexIndexMap>;
  using BlockedColor = property_traits_value_t<BlockedMap>;

  static BlockedColor blocked_false_color() { return color_traits<BlockedColor>::white(); }
  static BlockedColor blocked_true_color() { return color_traits<BlockedColor>::black(); }

  // This is used by the constructor to secure the assumption
  // documented above.
  [[nodiscard]] bool blocked_map_starts_all_unblocked() const {
    for (auto v : vertices(graph_)) {
      if (is_blocked(v)) {
        return false;
      }
    }
    return true;
  }

  // This is only used in the constructor to make sure the optimization of
  // sharing data structures between iterations does not break the code.
  [[nodiscard]] bool all_closed_rows_are_empty() const {
    for (const auto& row : closed_) {
      if (!row.empty()) {
        return false;
      }
    }
    return true;
  }

 public:
  hawick_circuits_from(const Graph& graph, Visitor& visitor, VertexIndexMap vim, Stack& stack, ClosedMatrix& closed,
                       std::size_t n_vertices, std::size_t max_length, GetAdjacentVertices /*for deduction only*/)
      : graph_(graph),
        visitor_(visitor),
        vim_(vim),
        stack_(stack),
        closed_(closed),
        blocked_(n_vertices, vim_),
        blocked_ref_(blocked_.ref()),
        max_length_(max_length) {
    assert(blocked_map_starts_all_unblocked());

    // Since sharing the data structures between iterations is
    // just an optimization, it must always be equivalent to
    // constructing new ones in this constructor.
    assert(stack_.empty());
    assert(closed_.size() == n_vertices);
    assert(all_closed_rows_are_empty());
  }

 private:
  //! @internal Return the index of a given vertex.
  auto index_of(Vertex v) const { return get(vim_, v); }

  //! @internal Return whether a vertex `v` is closed to a vertex `u`.
  bool is_closed_to(Vertex u, Vertex v) const {
    const auto& closed_to_u = closed_[index_of(u)];
    return contains(closed_to_u, v);
  }

  //! @internal Close a vertex `v` to a vertex `u`.
  void close_to(Vertex u, Vertex v) {
    assert(!is_closed_to(u, v));
    closed_[index_of(u)].push_back(v);
  }

  //! @internal Return whether a given vertex is blocked.
  bool is_blocked(Vertex v) const { return get(blocked_ref_, v) == blocked_true_color(); }

  //! @internal Block a given vertex.
  void block(Vertex v) { put(blocked_ref_, v, blocked_true_color()); }

  //! @internal Unblock a given vertex.
  void unblock(Vertex u) {
    put(blocked_ref_, u, blocked_false_color());
    auto& closed_to_u = closed_[index_of(u)];

    while (!closed_to_u.empty()) {
      const Vertex w = closed_to_u.back();
      closed_to_u.pop_back();
      if (is_blocked(w)) {
        unblock(w);
      }
    }
    assert(closed_to_u.empty());
  }

  //! @internal Main procedure as described in the paper.
  bool circuit(Vertex start, Vertex v) {
    bool found_circuit = false;
    stack_.push_back(v);
    block(v);

    // Truncate the search if any circuits would exceed max_length_.
    const bool truncate_search = (max_length_ > 0 && stack_.size() >= max_length_);

    // Cache some values that are used more than once in the function.
    const auto index_of_start = index_of(start);
    auto adj_vertices = GetAdjacentVertices()(v, graph_);
    for (Vertex w : adj_vertices) {
      // Since we're only looking in the subgraph induced by `start`
      // and the vertices with an index higher than `start`, we skip
      // any vertex that does not satisfy that.
      if (index_of(w) < index_of_start) {
        continue;
      }

      // If the last vertex is equal to `start`, we have a circuit.
      if (w == start) {
        // const to ensure the visitor does not modify the stack
        visitor_.cycle(std::as_const(stack_), graph_);
        found_circuit = true;
      } else if (truncate_search) {
        // If required, truncate the search before the subsequent
        // recursive call to circuit().
        continue;

      } else if (!is_blocked(w) && circuit(start, w)) {
        // If `w` is not blocked, we continue searching further down the
        // same path for a cycle with `w` in it.
        found_circuit = true;
      }
    }

    const bool finish_circuit = (found_circuit || truncate_search);
    if (finish_circuit) {
      unblock(v);
    } else {
      for (Vertex w : adj_vertices) {
        // Like above, we skip vertices that are not in the subgraph
        // we're considering.
        if (index_of(w) < index_of_start) {
          continue;
        }

        // If `v` is not closed to `w`, we make it so.
        if (!is_closed_to(w, v)) {
          close_to(w, v);
        }
      }
    }

    assert(v == stack_.back());
    stack_.pop_back();
    return finish_circuit;
  }

 public:
  void operator()(Vertex start) { circuit(start, start); }

 private:
  const Graph& graph_;
  Visitor& visitor_;
  VertexIndexMap vim_;
  Stack& stack_;
  ClosedMatrix& closed_;
  BlockedMap blocked_;
  property_map_ref<BlockedMap> blocked_ref_;
  std::size_t max_length_;
};

template <typename GetAdjacentVertices, concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor,
          concepts::ReadableVertexPropertyMap<Graph> VertexIndexMap>
requires concepts::AdjacencyGraph<Graph>
void call_hawick_circuits(const Graph& graph, Visitor /* by value */ visitor, VertexIndexMap vertex_index_map,
                          std::size_t max_length) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  using Stack = std::vector<Vertex>;
  using ClosedMatrix = std::vector<std::vector<Vertex>>;

  const std::size_t n_vertices = num_vertices(graph);
  Stack stack;
  stack.reserve(n_vertices);
  ClosedMatrix closed(n_vertices);

  for (auto v : vertices(graph)) {
    // Note1: The sub algorithm may NOT be reused once it has been
    // called.

    // Note2: We reuse the Stack and the ClosedMatrix (after clearing
    // them) in each iteration to avoid redundant destruction and
    // construction. It would be strictly equivalent to have these as
    // member variables of the sub algorithm.
    hawick_circuits_from sub_algo{graph,  visitor,    vertex_index_map, stack,
                                  closed, n_vertices, max_length,       GetAdjacentVertices{}};
    sub_algo(v);
    stack.clear();
    for (auto& row : closed) {
      row.clear();
    }
  }
}

template <typename GetAdjacentVertices, concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor>
requires concepts::AdjacencyGraph<Graph>
void call_hawick_circuits(Graph const& graph, Visitor&& visitor, std::size_t max_length) {
  call_hawick_circuits<GetAdjacentVertices>(graph, std::forward<Visitor>(visitor), get(vertex_index, graph),
                                            max_length);
}
}  // end namespace hawick_circuits_detail

//! Enumerate all the elementary circuits in a directed multigraph.
template <concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor,
          concepts::ReadableVertexPropertyMap<Graph> VertexIndexMap>
requires concepts::AdjacencyGraph<Graph>
void hawick_circuits(Graph&& graph, Visitor&& visitor, VertexIndexMap&& vertex_index_map, std::size_t max_length = 0) {
  hawick_circuits_detail::call_hawick_circuits<hawick_circuits_detail::get_all_adjacent_vertices>(
      std::forward<Graph>(graph), std::forward<Visitor>(visitor), std::forward<VertexIndexMap>(vertex_index_map),
      max_length);
}

template <concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor>
requires concepts::AdjacencyGraph<Graph>
void hawick_circuits(Graph&& graph, Visitor&& visitor, std::size_t max_length = 0) {
  hawick_circuits_detail::call_hawick_circuits<hawick_circuits_detail::get_all_adjacent_vertices>(
      std::forward<Graph>(graph), std::forward<Visitor>(visitor), max_length);
}

/*!
 * Same as `boost::hawick_circuits`, but duplicate circuits caused by parallel
 * edges will not be considered. Each circuit will be considered only once.
 */
template <concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor,
          concepts::ReadableVertexPropertyMap<Graph> VertexIndexMap>
requires concepts::AdjacencyGraph<Graph>
void hawick_unique_circuits(Graph&& graph, Visitor&& visitor, VertexIndexMap&& vertex_index_map,
                            std::size_t max_length = 0) {
  hawick_circuits_detail::call_hawick_circuits<hawick_circuits_detail::get_unique_adjacent_vertices>(
      std::forward<Graph>(graph), std::forward<Visitor>(visitor), std::forward<VertexIndexMap>(vertex_index_map),
      max_length);
}

template <concepts::VertexListGraph Graph, concepts::CycleVisitor<Graph> Visitor>
requires concepts::AdjacencyGraph<Graph>
void hawick_unique_circuits(Graph&& graph, Visitor&& visitor, std::size_t max_length = 0) {
  hawick_circuits_detail::call_hawick_circuits<hawick_circuits_detail::get_unique_adjacent_vertices>(
      std::forward<Graph>(graph), std::forward<Visitor>(visitor), max_length);
}

}  // namespace bagl

#endif  // BAGL_BAGL_HAWICK_CIRCUITS_H_
