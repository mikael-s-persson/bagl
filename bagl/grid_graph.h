// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRID_GRAPH_H_
#define BAGL_BAGL_GRID_GRAPH_H_

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <ranges>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/transform_value_property_map.h"

namespace bagl {

// Grid Graph

template <std::size_t Dimensions>
class grid_graph {
 private:
  using WrapDimensionArray = std::array<bool, Dimensions>;
  grid_graph() = default;

 public:
  using self = grid_graph<Dimensions>;

  // sizes
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  using degree_size_type = std::size_t;

  // descriptors
  using vertex_descriptor = std::array<std::size_t, Dimensions>;
  // Could be a pair, but safer to be a unique type.
  struct edge_descriptor {
    vertex_descriptor source;
    vertex_descriptor target;

    auto operator<=>(const edge_descriptor& rhs) const {
      return std::pair{source, target} <=> std::pair{rhs.source, rhs.target};
    }
  };

  // categories
  using directed_category = directed_tag;
  using edge_parallel_category = disallow_parallel_edge_tag;
  struct traversal_category : virtual public incidence_graph_tag,
                              virtual public adjacency_graph_tag,
                              virtual public vertex_list_graph_tag,
                              virtual public edge_list_graph_tag,
                              virtual public bidirectional_graph_tag,
                              virtual public adjacency_matrix_tag {};

  static vertex_descriptor null_vertex() {
    vertex_descriptor maxed_out_vertex;
    std::ranges::fill(maxed_out_vertex, std::numeric_limits<vertices_size_type>::max());
    return maxed_out_vertex;
  }

  // Constructor that defaults to no wrapping for all dimensions.
  explicit grid_graph(vertex_descriptor dimension_lengths) : dimension_lengths_(dimension_lengths) {
    std::ranges::fill(wrap_dimension_, false);

    precalculate();
  }

  // Constructor that allows for wrapping to be specified for all
  // dimensions at once.
  grid_graph(vertex_descriptor dimension_lengths, bool wrap_all_dimensions) : dimension_lengths_(dimension_lengths) {
    std::ranges::fill(wrap_dimension_, wrap_all_dimensions);

    precalculate();
  }

  // Constructor that allows for individual dimension wrapping to be
  // specified.
  grid_graph(vertex_descriptor dimension_lengths, WrapDimensionArray wrap_dimension)
      : dimension_lengths_(dimension_lengths), wrap_dimension_(wrap_dimension) {
    precalculate();
  }

  // Returns the number of dimensions in the graph
  [[nodiscard]] std::size_t dimensions() const { return Dimensions; }

  // Returns the length of dimension [dimension_index]
  [[nodiscard]] std::size_t length(std::size_t dimension) const { return dimension_lengths_[dimension]; }

  // Returns a value indicating if dimension [dimension_index] wraps
  [[nodiscard]] bool wrapped(std::size_t dimension) const { return wrap_dimension_[dimension]; }

  // Gets the vertex that is [distance] units ahead of [vertex] in
  // dimension [dimension_index].
  vertex_descriptor next(vertex_descriptor vertex, std::size_t dimension_index, std::size_t distance = 1) const {
    std::size_t new_position = vertex[dimension_index] + distance;

    if (wrapped(dimension_index)) {
      new_position %= length(dimension_index);
    } else {
      // Stop at the end of this dimension if necessary.
      new_position = std::min(new_position, length(dimension_index) - 1);
    }

    vertex[dimension_index] = new_position;

    return (vertex);
  }

  // Gets the vertex that is [distance] units behind [vertex] in
  // dimension [dimension_index].
  vertex_descriptor previous(vertex_descriptor vertex, std::size_t dimension_index, std::size_t distance = 1) const {
    // We're assuming that vertices_size_type is unsigned, so we
    // need to be careful about the math.
    vertex[dimension_index] =
        (distance > vertex[dimension_index])
            ? (wrapped(dimension_index) ? (length(dimension_index) - (distance % length(dimension_index))) : 0)
            : vertex[dimension_index] - distance;

    return (vertex);
  }

 protected:
  // Returns the number of vertices in the graph
  [[nodiscard]] std::size_t num_vertices() const { return num_vertices_; }

  // Returns the number of edges in the graph
  [[nodiscard]] std::size_t num_edges() const { return num_edges_; }

  // Returns the number of edges in dimension [dimension_index]
  [[nodiscard]] std::size_t num_edges(std::size_t dimension_index) const { return edge_count_[dimension_index]; }

  // Returns the index of [vertex] (See also vertex_at)
  std::size_t index_of(vertex_descriptor vertex) const {
    std::size_t vertex_id = 0;
    std::size_t index_multiplier = 1;

    for (std::size_t dimension_index = 0; dimension_index < Dimensions; ++dimension_index) {
      vertex_id += (vertex[dimension_index] * index_multiplier);
      index_multiplier *= length(dimension_index);
    }

    return vertex_id;
  }

  // Returns the vertex whose index is [vertex_id] (See also
  // index_of(vertex_descriptor))
  vertex_descriptor vertex_at(std::size_t vertex_id) const {
    vertex_descriptor vertex{};
    std::size_t index_divider = 1;

    for (std::size_t dimension_index = 0; dimension_index < Dimensions; ++dimension_index) {
      vertex[dimension_index] = (vertex_id / index_divider) % length(dimension_index);

      index_divider *= length(dimension_index);
    }

    return vertex;
  }

  // Returns the edge whose index is [edge_id] (See also
  // index_of(edge_descriptor)).  NOTE: The index mapping is
  // dependent upon dimension wrapping.
  edge_descriptor edge_at(std::size_t edge_id) const {
    // Edge indices are sorted into bins by dimension
    std::size_t dimension_index = 0;
    std::size_t dimension_edges = num_edges(0);

    while (edge_id >= dimension_edges) {
      edge_id -= dimension_edges;
      ++dimension_index;
      dimension_edges = num_edges(dimension_index);
    }

    vertex_descriptor vertex_source{};
    vertex_descriptor vertex_target{};
    bool is_forward = ((edge_id / (num_edges(dimension_index) / 2)) == 0);

    if (wrapped(dimension_index)) {
      vertex_source = vertex_at(edge_id % num_vertices());
      vertex_target = is_forward ? next(vertex_source, dimension_index) : previous(vertex_source, dimension_index);
      return {vertex_source, vertex_target};
    }

    // Dimensions can wrap arbitrarily, so an index needs to be
    // computed in a more complex manner.  This is done by
    // grouping the edges for each dimension together into "bins"
    // and considering [edge_id] as an offset into the bin.
    // Each bin consists of two parts: the "forward" looking edges
    // and the "backward" looking edges for the dimension.

    std::size_t vertex_offset = edge_id % num_edges(dimension_index);

    // Consider vertex_offset an index into the graph's vertex
    // space but with the dimension [dimension_index] reduced in
    // size by one.
    std::size_t index_divider = 1;

    for (std::size_t dimension_index_iter = 0; dimension_index_iter < Dimensions; ++dimension_index_iter) {
      std::size_t dimension_length =
          (dimension_index_iter == dimension_index) ? length(dimension_index_iter) - 1 : length(dimension_index_iter);

      vertex_source[dimension_index_iter] = (vertex_offset / index_divider) % dimension_length;

      index_divider *= dimension_length;
    }

    if (is_forward) {
      vertex_target = next(vertex_source, dimension_index);
    } else {
      // Shift forward one more unit in the dimension for backward
      // edges since the algorithm above will leave us one behind.
      vertex_target = vertex_source;
      ++vertex_source[dimension_index];
    }

    return {vertex_source, vertex_target};
  }

  // Returns the index for [edge] (See also edge_at)
  std::size_t index_of(edge_descriptor edge) const {
    assert(edge.source != edge.target);

    // Determine the dimension where the source and target vertices
    // differ (should only be one if this is a valid edge).
    std::size_t different_dimension_index = 0;

    while (edge.source[different_dimension_index] == edge.target[different_dimension_index]) {
      ++different_dimension_index;
    }

    std::size_t edge_id = 0;

    // Offset the edge index into the appropriate "bin" (see edge_at
    // for a more in-depth description).
    for (std::size_t dimension_index = 0; dimension_index < different_dimension_index; ++dimension_index) {
      edge_id += num_edges(dimension_index);
    }

    // Get the position of both vertices in the differing dimension.
    std::size_t source_position = edge.source[different_dimension_index];
    std::size_t target_position = edge.target[different_dimension_index];

    // Determine if edge is forward or backward
    bool is_forward = true;

    if (wrapped(different_dimension_index)) {
      // If the dimension is wrapped, an edge is going backward if
      // either A: its target precedes the source in the differing
      // dimension and the vertices are adjacent or B: its source
      // precedes the target and they're not adjacent.
      if (((target_position < source_position) && (source_position == 1 + target_position)) ||
          ((source_position < target_position) && (target_position > 1 + source_position))) {
        is_forward = false;
      }
    } else if (target_position < source_position) {
      is_forward = false;
    }

    // "Backward" edges are in the second half of the bin.
    if (!is_forward) {
      edge_id += num_edges(different_dimension_index) / 2;
    }

    // Finally, apply the vertex offset
    if (wrapped(different_dimension_index)) {
      edge_id += index_of(edge.source);
    } else {
      std::size_t index_multiplier = 1;

      if (!is_forward) {
        --edge.source[different_dimension_index];
      }

      for (std::size_t dimension_index = 0; dimension_index < Dimensions; ++dimension_index) {
        edge_id += edge.source[dimension_index] * index_multiplier;
        index_multiplier *=
            (dimension_index == different_dimension_index ? length(dimension_index) - 1 : length(dimension_index));
      }
    }

    return edge_id;
  }

  // Returns the number of out-edges for [vertex]
  std::size_t out_degree(vertex_descriptor vertex) const {
    std::size_t out_edge_count = 0;

    for (std::size_t dimension_index = 0; dimension_index < Dimensions; ++dimension_index) {
      // If the vertex is on the edge of this dimension, then its
      // number of out edges is dependent upon whether the dimension
      // wraps or not.
      if ((vertex[dimension_index] == 0) || (vertex[dimension_index] + 1 == length(dimension_index))) {
        out_edge_count += (wrapped(dimension_index) ? 2 : 1);
      } else {
        // Next and previous edges, regardless or wrapping
        out_edge_count += 2;
      }
    }

    return out_edge_count;
  }

  // Returns an out-edge for [vertex] by index. Indices are in the
  // range [0, out_degree(vertex)).
  edge_descriptor out_edge_at(vertex_descriptor vertex, std::size_t out_edge_id) const {
    std::size_t edges_left = out_edge_id + 1;
    std::size_t dimension_index = 0;
    bool is_forward = false;

    // Walks the out edges of [vertex] and accommodates for dimension
    // wrapping.
    while (edges_left > 0) {
      if (!wrapped(dimension_index)) {
        if (!is_forward && (vertex[dimension_index] == 0)) {
          is_forward = true;
          continue;
        }
        if (is_forward && (vertex[dimension_index] == (length(dimension_index) - 1))) {
          is_forward = false;
          ++dimension_index;
          continue;
        }
      }

      --edges_left;

      if (edges_left > 0) {
        is_forward = !is_forward;

        if (!is_forward) {
          ++dimension_index;
        }
      }
    }

    return {vertex, is_forward ? next(vertex, dimension_index) : previous(vertex, dimension_index)};
  }

  // Returns the number of in-edges for [vertex]
  std::size_t in_degree(vertex_descriptor vertex) const { return out_degree(vertex); }

  // Returns an in-edge for [vertex] by index. Indices are in the
  // range [0, in_degree(vertex)).
  edge_descriptor in_edge_at(vertex_descriptor vertex, std::size_t in_edge_id) const {
    edge_descriptor out_edge = out_edge_at(vertex, in_edge_id);
    return {out_edge.target, out_edge.source};
  }

  // Pre-computes the number of vertices and edges
  void precalculate() {
    num_vertices_ = std::accumulate(dimension_lengths_.begin(), dimension_lengths_.end(), std::size_t{1}, std::multiplies<>());

    // Calculate number of edges in each dimension
    num_edges_ = 0;

    for (std::size_t dimension_index = 0; dimension_index < Dimensions; ++dimension_index) {
      if (wrapped(dimension_index)) {
        edge_count_[dimension_index] = num_vertices() * 2;
      } else {
        edge_count_[dimension_index] = (num_vertices() - (num_vertices() / length(dimension_index))) * 2;
      }

      num_edges_ += num_edges(dimension_index);
    }
  }

  const vertex_descriptor dimension_lengths_{};
  WrapDimensionArray wrap_dimension_{};
  std::size_t num_vertices_ = 0;

  std::array<std::size_t, Dimensions> edge_count_{};
  std::size_t num_edges_ = 0;

 public:
  // VertexListGraph

  friend auto vertices(const self& graph) {
    return std::ranges::iota_view{std::size_t{0}, graph.num_vertices()} |
           std::views::transform([&graph](std::size_t vertex_id) { return graph.vertex_at(vertex_id); });
  }

  friend auto num_vertices(const self& graph) { return graph.num_vertices(); }

  friend auto vertex(std::size_t vertex_id, const self& graph) { return graph.vertex_at(vertex_id); }

  // IncidenceGraph

  friend auto out_edges(vertex_descriptor vertex, const self& graph) {
    return std::ranges::iota_view{std::size_t{0}, graph.out_degree(vertex)} |
           std::views::transform(
               [vertex, &graph](std::size_t out_edge_id) { return graph.out_edge_at(vertex, out_edge_id); });
  }

  friend std::size_t out_degree(vertex_descriptor vertex, const self& graph) { return graph.out_degree(vertex); }

  friend auto out_edge_at(vertex_descriptor vertex, std::size_t out_edge_id, const self& graph) {
    return graph.out_edge_at(vertex, out_edge_id);
  }

  // AdjacencyGraph

  friend auto adjacent_vertices(vertex_descriptor vertex, const self& graph) {
    return std::ranges::iota_view{std::size_t{0}, graph.out_degree(vertex)} |
           std::views::transform(
               [vertex, &graph](std::size_t out_edge_id) { return graph.out_edge_at(vertex, out_edge_id).target; });
  }

  // EdgeListGraph

  friend std::size_t num_edges(const self& graph) { return (graph.num_edges()); }

  friend auto edge_at(std::size_t edge_id, const self& graph) { return graph.edge_at(edge_id); }

  friend auto edges(const self& graph) {
    return std::ranges::iota_view{std::size_t{0}, graph.num_edges()} |
           std::views::transform([&graph](std::size_t edge_id) { return graph.edge_at(edge_id); });
  }

  // BiDirectionalGraph

  friend auto in_edges(vertex_descriptor vertex, const self& graph) {
    return std::ranges::iota_view{std::size_t{0}, graph.in_degree(vertex)} |
           std::views::transform(
               [vertex, &graph](std::size_t out_edge_id) { return graph.in_edge_at(vertex, out_edge_id).target; });
  }

  friend std::size_t in_degree(vertex_descriptor vertex, const self& graph) { return graph.in_degree(vertex); }

  friend std::size_t degree(vertex_descriptor vertex, const self& graph) { return graph.out_degree(vertex) * 2; }

  friend auto in_edge_at(vertex_descriptor vertex, std::size_t in_edge_id, const self& graph) {
    return graph.in_edge_at(vertex, in_edge_id);
  }

  // Adjacency Matrix

  friend auto edge(vertex_descriptor source_vertex, vertex_descriptor target_vertex, const self& graph) {
    auto edge_exists = std::pair{edge_descriptor{source_vertex, target_vertex}, false};

    for (std::size_t dimension_index = 0; dimension_index < Dimensions; ++dimension_index) {
      std::size_t source_dim = source_vertex[dimension_index];
      std::size_t dest_dim = target_vertex[dimension_index];

      std::size_t dim_difference = (source_dim > dest_dim) ? (source_dim - dest_dim) : (dest_dim - source_dim);

      if (dim_difference > 0) {
        // If we've already found a valid edge, this would mean that
        // the vertices are really diagonal across dimensions and
        // therefore not connected.
        if (edge_exists.second) {
          edge_exists.second = false;
          break;
        }

        // If the difference is one, the vertices are right next to
        // each other and the edge is valid.  The edge is still
        // valid, though, if the dimension wraps and the vertices
        // are on opposite ends.
        if ((dim_difference == 1) || (graph.wrapped(dimension_index) &&
                                      (((source_dim == 0) && (dest_dim + 1 == graph.length(dimension_index))) ||
                                       ((dest_dim == 0) && (source_dim + 1 == graph.length(dimension_index)))))) {
          edge_exists.second = true;
          // Stay in the loop to check for diagonal vertices.
        } else {
          // Stop checking - the vertices are too far apart.
          edge_exists.second = false;
          break;
        }
      }

    }  // for dimension_index

    return edge_exists;
  }

  // Index Property Map Functions

  friend std::size_t get(vertex_index_t /*unused*/, const self& graph, vertex_descriptor vertex) {
    return graph.index_of(vertex);
  }

  friend std::size_t get(edge_index_t /*unused*/, const self& graph, edge_descriptor edge) {
    return graph.index_of(edge);
  }

  friend auto get(vertex_index_t /*unused*/, const self& graph) {
    return transformed_property_map<vertex_descriptor>([&graph](vertex_descriptor v) { return graph.index_of(v); });
  }

  friend auto get(edge_index_t /*unused*/, const self& graph) {
    return transformed_property_map<edge_descriptor>([&graph](edge_descriptor e) { return graph.index_of(e); });
  }

  friend auto get(edge_reverse_t /*unused*/, const self& /*graph*/) {
    return transformed_property_map<edge_descriptor>([](edge_descriptor e) { return edge_descriptor{e.target, e.source}; });
  }

};  // grid_graph


// Index Property Map

template <std::size_t Dimensions>
struct property_map<grid_graph<Dimensions>, vertex_index_t> {
  using type = decltype(get(vertex_index, std::declval<const grid_graph<Dimensions>&>()));
  using const_type = type;
};

template <std::size_t Dimensions>
struct property_map<grid_graph<Dimensions>, edge_index_t> {
  using type = decltype(get(edge_index, std::declval<const grid_graph<Dimensions>&>()));
  using const_type = type;
};

// Reverse Edge Property Map

template <std::size_t Dimensions>
struct property_map<grid_graph<Dimensions>, edge_reverse_t> {
  using type = decltype(get(edge_reverse, std::declval<const grid_graph<Dimensions>&>()));
  using const_type = type;
};

}  // namespace bagl

#endif  // BAGL_BAGL_GRID_GRAPH_H_
