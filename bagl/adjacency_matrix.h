// Copyright 2001 University of Notre Dame.
// Copyright 2006 Trustees of Indiana University
// Authors: Jeremy G. Siek and Douglas Gregor <dgregor@cs.indiana.edu>
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ADJACENCY_MATRIX_H_
#define BAGL_BAGL_ADJACENCY_MATRIX_H_

#include <cassert>
#include <iterator>
#include <limits>
#include <ranges>
#include <type_traits>
#include <vector>

#include "bagl/graph_mutability_traits.h"
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/more_property_maps.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

namespace adj_matrix_detail {

class undir_edges_range : public std::ranges::view_interface<undir_edges_range> {
 public:
  using self_range = undir_edges_range;

  explicit undir_edges_range(std::size_t num_v) : num_v_(num_v) {}
  undir_edges_range() = default;

  struct iterator {
    using value_type = std::pair<std::size_t, std::size_t>;
    using reference = const value_type&;
    using pointer = const value_type*;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::bidirectional_iterator_tag;

    iterator& operator++() {
      if (cur_v == cur_u) {
        // If at diagonal, move to start of next row.
        cur_v = 0;
        ++cur_u;
      } else {
        ++cur_v;
      }
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator& operator--() {
      if (cur_v == 0) {
        // If at start of row, move to diagonal of previous row.
        --cur_u;
        cur_v = cur_u;
      } else {
        --cur_v;
      }
      return *this;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      --(*this);
      return tmp;
    }
    [[nodiscard]] bool operator==(const iterator& rhs) const { return (cur_u == rhs.cur_u) && (cur_v == rhs.cur_v); }
    [[nodiscard]] bool operator!=(const iterator& rhs) const { return !(*this == rhs); }
    value_type operator*() const { return {cur_u, cur_v}; }
    const value_type* operator->() const {
      assert(false && "This function only exists to comply with legacy iterator concepts.");
      static auto referrable = *(*this);
      referrable = *(*this);
      return &referrable;
    }

    iterator(std::size_t u, std::size_t v) : cur_u(u), cur_v(v) {}

    iterator() = default;

    std::size_t cur_u = 0;
    std::size_t cur_v = 0;
  };

  [[nodiscard]] iterator begin() const { return iterator{0, 0}; }
  [[nodiscard]] iterator begin() { return iterator{0, 0}; }
  [[nodiscard]] iterator end() const { return iterator{num_v_, 0}; }
  [[nodiscard]] iterator end() { return iterator{num_v_, 0}; }

  std::size_t num_v_;
};

}  // namespace adj_matrix_detail

template <bool IsDirected>
struct adjacency_matrix_edge_descriptor {
  std::size_t source = 0;
  std::size_t target = 0;
  std::size_t edge_index = 0;

  auto operator<=>(const adjacency_matrix_edge_descriptor& rhs) const { return edge_index <=> rhs.edge_index; }

  adjacency_matrix_edge_descriptor() = default;
  adjacency_matrix_edge_descriptor(std::size_t u, std::size_t v, std::size_t num_v) : source(u), target(v) {
    if constexpr (IsDirected) {
      // Square matrix index.
      edge_index = u * num_v + v;
    } else {
      if (v > u) {
        std::swap(u, v);
      }
      // Triangular/symmetric matrix index.
      edge_index = u * (u + 1) / 2 + v;
    }
  }
};

//=========================================================================
// Adjacency Matrix Traits
template <typename Directed = directed_s>
class adjacency_matrix_traits {
 public:
  // The bidirectionalS tag is not allowed with the adjacency_matrix
  // graph type. Instead, use directed_s, which also provides the
  // functionality required for a Bidirectional Graph (in_edges,
  // in_degree, etc.).
  static_assert(!std::is_same_v<Directed, bidirectional_s>);

  static constexpr bool is_directed_v = std::is_same_v<Directed, directed_s>;

  using directed_category = std::conditional_t<is_directed_v, bidirectional_tag, undirected_tag>;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertex_descriptor = std::size_t;

  using edge_descriptor = adjacency_matrix_edge_descriptor<is_directed_v>;
};

struct adjacency_matrix_class_tag {};

struct adj_matrix_traversal_tag : public virtual adjacency_matrix_tag,
                                  public virtual vertex_list_graph_tag,
                                  public virtual incidence_graph_tag,
                                  public virtual adjacency_graph_tag,
                                  public virtual edge_list_graph_tag {};

//=========================================================================
// Adjacency Matrix Class
template <typename Directed = directed_s, typename VertexProperty = no_property, typename EdgeProperty = no_property,
          typename GraphProperty = no_property, typename Allocator = std::allocator<bool>>
class adjacency_matrix {
  using self = adjacency_matrix;
  using Traits = adjacency_matrix_traits<Directed>;

 public:
  static constexpr bool is_directed_v = Traits::is_directed_v;

  using graph_property_type = GraphProperty;
  using graph_bundled = lookup_one_property_t<GraphProperty, graph_bundle_t>;

  using vertex_property_type = VertexProperty;
  using vertex_bundled = lookup_one_property_t<VertexProperty, vertex_bundle_t>;

  using edge_property_type = EdgeProperty;
  using edge_bundled = lookup_one_property_t<EdgeProperty, edge_bundle_t>;

 private:
  using StoredEdge = std::optional<edge_property_type>;
  using Alloc = typename std::allocator_traits<Allocator>::template rebind_alloc<StoredEdge>;
  using Matrix = std::vector<StoredEdge, Alloc>;
  using size_type = std::size_t;

 public:
  // Graph concept required types
  using vertex_descriptor = typename Traits::vertex_descriptor;
  using edge_descriptor = typename Traits::edge_descriptor;
  using directed_category = typename Traits::directed_category;
  using edge_parallel_category = typename Traits::edge_parallel_category;
  using traversal_category = adj_matrix_traversal_tag;

  static vertex_descriptor null_vertex() { return (std::numeric_limits<vertex_descriptor>::max)(); }

  using degree_size_type = size_type;
  using vertices_size_type = size_type;
  using edges_size_type = size_type;

  // PropertyGraph required types
  using graph_tag = adjacency_matrix_class_tag;

  // Constructor required by MutableGraph
  explicit adjacency_matrix(vertices_size_type n_vertices, const GraphProperty& p = GraphProperty())
      : matrix_(edge_descriptor{n_vertices, 0, n_vertices}.edge_index, std::nullopt),
        vertex_properties_(n_vertices),
        property_(p) {}

  // Construct from a given number of vertices and an edge range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>>
  adjacency_matrix(vertices_size_type num_vertices, const EdgeRange& e_range, graph_property_type graph_prop = {})
      : matrix_(edge_descriptor{num_vertices, 0, num_vertices}.edge_index, std::nullopt),
        vertex_properties_(num_vertices),
        property_(std::move(graph_prop)) {
    for (auto [u, v] : e_range) {
      add_edge(u, v, edge_property_type{}, *this);
    }
  }

  // Construct from a given number of vertices and an edge and edge-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>> &&
      std::convertible_to<std::ranges::range_reference_t<EdgePropRange>, edge_property_type>
      adjacency_matrix(vertices_size_type num_vertices, const EdgeRange& e_range, const EdgePropRange& ep_range,
                       graph_property_type graph_prop = {})
      : matrix_(edge_descriptor{num_vertices, 0, num_vertices}.edge_index, std::nullopt),
        vertex_properties_(num_vertices),
        property_(std::move(graph_prop)) {
    for (const auto& [e, ep] : zip_range(e_range, ep_range)) {
      auto [u, v] = e;
      add_edge(u, v, ep, *this);
    }
  }

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>> &&
      std::convertible_to<std::ranges::range_reference_t<VertexPropRange>, vertex_property_type>
      adjacency_matrix(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                       graph_property_type graph_prop = {})
      : matrix_(edge_descriptor{num_vertices, 0, num_vertices}.edge_index, std::nullopt),
        vertex_properties_(num_vertices),
        property_(std::move(graph_prop)) {
    for (auto& [v_store, v_arg] : zip_range(vertex_properties_, vp_range)) {
      v_store = v_arg;
    }

    for (auto [u, v] : e_range) {
      add_edge(u, v, edge_property_type{}, *this);
    }
  }

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange,
            std::ranges::input_range EdgePropRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>> &&
      std::convertible_to<std::ranges::range_reference_t<VertexPropRange>, vertex_property_type> &&
      std::convertible_to<std::ranges::range_reference_t<EdgePropRange>, edge_property_type>
      adjacency_matrix(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                       const EdgePropRange& ep_range, graph_property_type graph_prop = {})
      : matrix_(edge_descriptor{num_vertices, 0, num_vertices}.edge_index, std::nullopt),
        vertex_properties_(num_vertices),
        property_(std::move(graph_prop)) {
    for (auto& [v_store, v_arg] : zip_range(vertex_properties_, vp_range)) {
      v_store = v_arg;
    }

    for (const auto& [e, ep] : zip_range(e_range, ep_range)) {
      auto [u, v] = e;
      add_edge(u, v, ep, *this);
    }
  }

  // Indexing operator. Returns a reference to the vertex-bundle associated to the given vertex descriptor.
  [[nodiscard]] auto& operator[](vertex_descriptor v) {
    return get_property_value(vertex_properties_[v], vertex_bundle);
  }

  // Indexing operator. Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
  [[nodiscard]] const auto& operator[](vertex_descriptor v) const {
    return get_property_value(vertex_properties_[v], vertex_bundle);
  }

  // Indexing operator. Returns a reference to the edge-bundle associated to the given edge descriptor.
  [[nodiscard]] auto& operator[](const edge_descriptor& e) {
    auto& ep = get_edge(e);
    assert(ep.has_value());
    return get_property_value(ep.value(), edge_bundle);
  }

  // Indexing operator. Returns a const-reference to the edge-bundle associated to the given edge descriptor.
  [[nodiscard]] const auto& operator[](const edge_descriptor& e) const {
    auto& ep = get_edge(e);
    assert(ep.has_value());
    return get_property_value(ep.value(), edge_bundle);
  }

  // Indexing operator. Returns a reference to the graph-bundle associated to the graph.
  [[nodiscard]] auto& operator[](graph_bundle_t /*unused*/) { return get_property_value(property_, graph_bundle); }

  // Indexing operator. Returns a const-reference to the graph-bundle associated to the graph.
  [[nodiscard]] const auto& operator[](graph_bundle_t /*unused*/) const {
    return get_property_value(property_, graph_bundle);
  }

  // Get a reference to the vertex-property associated to the given vertex descriptor.
  [[nodiscard]] auto& get_property(vertex_descriptor v) { return vertex_properties_[v]; }
  // Get a const-reference to the vertex-property associated to the given vertex descriptor.
  [[nodiscard]] const auto& get_property(vertex_descriptor v) const { return vertex_properties_[v]; }

  // Get a reference to the edge-property associated to the given edge descriptor.
  [[nodiscard]] auto& get_property(const edge_descriptor& e) { return get_edge(e).value(); }
  // Get a const-reference to the edge-property associated to the given edge descriptor.
  [[nodiscard]] const auto& get_property(const edge_descriptor& e) const { return get_edge(e).value(); }

  // Get a reference to the graph-property associated to the graph.
  [[nodiscard]] auto& get_property(graph_all_t /*unused*/) { return property_; }
  // Get a const-reference to the graph-property associated to the graph.
  [[nodiscard]] const auto& get_property(graph_all_t /*unused*/) const { return property_; }

  // private: if friends worked, these would be private

  [[nodiscard]] edge_descriptor make_edge_descriptor(vertex_descriptor u, vertex_descriptor v) const {
    return {u, v, vertex_properties_.size()};
  }
  [[nodiscard]] bool edge_exists(const edge_descriptor& e) const { return matrix_[e.edge_index].has_value(); }
  [[nodiscard]] const auto& get_edge(const edge_descriptor& e) const { return matrix_[e.edge_index]; }
  [[nodiscard]] auto& get_edge(const edge_descriptor& e) { return matrix_[e.edge_index]; }

  Matrix matrix_;
  std::vector<vertex_property_type> vertex_properties_;
  size_type num_edges_ = 0;
  graph_property_type property_;
};

#define BAGL_ADJACENCY_MATRIX_PARAMS typename D, typename VP, typename EP, typename GP, typename A
#define BAGL_ADJACENCY_MATRIX adjacency_matrix<D, VP, EP, GP, A>

//=========================================================================
// Functions required by the AdjacencyMatrix concept

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto edge(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v,
          const BAGL_ADJACENCY_MATRIX& g) {
  auto e = g.make_edge_descriptor(u, v);
  return std::pair{e, g.edge_exists(e)};
}

//=========================================================================
// Functions required by the IncidenceGraph concept

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto out_edges(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, const BAGL_ADJACENCY_MATRIX& g) {
  return std::views::iota(std::size_t{0}, num_vertices(g)) |
         std::views::transform([&g, u](std::size_t v) { return g.make_edge_descriptor(u, v); }) |
         std::views::filter([&g](const auto& e) { return g.edge_exists(e); });
}

// O(N)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto out_degree(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, const BAGL_ADJACENCY_MATRIX& g) {
  return std::ranges::distance(out_edges(u, g));
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto source(const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e, const BAGL_ADJACENCY_MATRIX& /*unused*/) {
  return e.source;
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto target(const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e, const BAGL_ADJACENCY_MATRIX& /*unused*/) {
  return e.target;
}

//=========================================================================
// Functions required by the BidirectionalGraph concept

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto in_edges(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v, const BAGL_ADJACENCY_MATRIX& g) {
  return std::views::iota(std::size_t{0}, num_vertices(g)) |
         std::views::transform([&g, v](std::size_t u) { return g.make_edge_descriptor(u, v); }) |
         std::views::filter([&g](const auto& e) { return g.edge_exists(e); });
}

// O(N)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto in_degree(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, const BAGL_ADJACENCY_MATRIX& g) {
  return std::ranges::distance(in_edges(u, g));
}

//=========================================================================
// Functions required by the AdjacencyGraph concept

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto adjacent_vertices(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, const BAGL_ADJACENCY_MATRIX& g) {
  return out_edges(u, g) | std::views::transform([](const auto& e) { return e.target; });
}

//=========================================================================
// Functions required by the VertexListGraph concept

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto vertices(const BAGL_ADJACENCY_MATRIX& g) {
  return std::views::iota(std::size_t{0}, num_vertices(g));
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
std::size_t num_vertices(const BAGL_ADJACENCY_MATRIX& g) {
  return g.vertex_properties_.size();
}

//=========================================================================
// Functions required by the EdgeListGraph concept

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto edges(const BAGL_ADJACENCY_MATRIX& g) {
  if constexpr (BAGL_ADJACENCY_MATRIX::is_directed_v) {
    // Iterate edge indices from [0, n^2), split into (u,v) for square matrix, then map to edges and filter.
    std::size_t num_v = num_vertices(g);
    return std::views::iota(std::size_t{0}, num_v * num_v) | std::views::transform([num_v](std::size_t id) {
             return std::pair{id / num_v, id % num_v};
           }) |
           std::views::transform([&g](auto uv) { return g.make_edge_descriptor(uv.first, uv.second); }) |
           std::views::filter([&g](const auto& e) { return g.edge_exists(e); });
  } else {
    // Iterate (u,v) indices with a triangular number iterator, then map to edges and filter.
    return adj_matrix_detail::undir_edges_range(num_vertices(g)) |
           std::views::transform([&g](auto uv) { return g.make_edge_descriptor(uv.first, uv.second); }) |
           std::views::filter([&g](const auto& e) { return g.edge_exists(e); });
  }
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto num_edges(const BAGL_ADJACENCY_MATRIX& g) {
  return g.num_edges_;
}

//=========================================================================
// Functions required by the MutableGraph concept

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS, typename EProp>
auto add_edge(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v,
              EProp&& ep, BAGL_ADJACENCY_MATRIX& g) {
  auto e = g.make_edge_descriptor(u, v);
  if (g.edge_exists(e)) {
    return std::pair{e, false};
  }
  g.get_edge(e) = std::forward<EProp>(ep);
  ++(g.num_edges_);
  return std::pair{e, true};
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto add_edge(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v,
              BAGL_ADJACENCY_MATRIX& g) {
  return add_edge(u, v, EP{}, g);
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
void remove_edge(typename BAGL_ADJACENCY_MATRIX::edge_descriptor e, BAGL_ADJACENCY_MATRIX& g) {
  if (g.edge_exists(e)) {
    --(g.num_edges_);
    g.get_edge(e) = std::nullopt;
  }
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
void remove_edge(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u,
                 typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v, BAGL_ADJACENCY_MATRIX& g) {
  remove_edge(g.make_edge_descriptor(u, v), g);
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS, typename EProp>
void remove_edge(typename BAGL_ADJACENCY_MATRIX::edge_descriptor e, EProp* ep, BAGL_ADJACENCY_MATRIX& g) {
  if (g.edge_exists(e)) {
    --(g.num_edges_);
    if (ep != nullptr) {
      *ep = std::move(g.get_edge(e).value());
    }
    g.get_edge(e) = std::nullopt;
  }
}

// O(1)
template <BAGL_ADJACENCY_MATRIX_PARAMS, typename EProp>
void remove_edge(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u,
                 typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v, EProp* ep, BAGL_ADJACENCY_MATRIX& g) {
  remove_edge(g.make_edge_descriptor(u, v), ep, g);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto add_vertex(BAGL_ADJACENCY_MATRIX& g) {
  // UNDER CONSTRUCTION
  assert(false && "Not implemented.");
  return *vertices(g).begin();
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename VP2>
auto add_vertex(const VP2&& /*vp*/, BAGL_ADJACENCY_MATRIX& g) {
  // UNDER CONSTRUCTION
  assert(false && "Not implemented.");
  return *vertices(g).begin();
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
void remove_vertex(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor /*u*/, BAGL_ADJACENCY_MATRIX& /*g*/) {
  // UNDER CONSTRUCTION
  assert(false && "Not implemented.");
}

// O(V)
template <BAGL_ADJACENCY_MATRIX_PARAMS>
void clear_vertex(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor u, BAGL_ADJACENCY_MATRIX& g) {
  for (auto e : out_edges(u, g)) {
    remove_edge(e, g);
  }
  if constexpr (BAGL_ADJACENCY_MATRIX::is_directed_v) {
    for (auto e : in_edges(u, g)) {
      remove_edge(e, g);
    }
  }
}

/***********************************************************************************************
 *                             Property Maps (from bundles)
 * ********************************************************************************************/

// Returns a const-reference to the vertex-bundle associated to the given vertex descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
const auto& get(const BAGL_ADJACENCY_MATRIX& g, typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v) {
  return g[v];
}
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto& get(BAGL_ADJACENCY_MATRIX& g, typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v) {
  return g[v];
}

// Returns a const-reference to the edge-bundle associated to the given edge descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
const auto& get(const BAGL_ADJACENCY_MATRIX& g, const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e) {
  return g[e];
}
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto& get(BAGL_ADJACENCY_MATRIX& g, const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e) {
  return g[e];
}

// Returns a const-reference to the graph-bundle associated to the graph.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
const auto& get(const BAGL_ADJACENCY_MATRIX& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto& get(BAGL_ADJACENCY_MATRIX& g, graph_bundle_t /*unused*/) {
  return g[graph_bundle];
}

// Sets the vertex-bundle associated to the given vertex descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS, typename VProp>
void put(BAGL_ADJACENCY_MATRIX& g, typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v, VProp&& value) {
  g[v] = std::forward<VProp>(value);
}

// Sets the edge-bundle associated to the given edge descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS, typename EProp>
void put(BAGL_ADJACENCY_MATRIX& g, const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e, EProp&& value) {
  g[e] = std::forward<EProp>(value);
}

// Sets the graph-bundle associated to the graph.
template <BAGL_ADJACENCY_MATRIX_PARAMS, typename GProp>
void put(BAGL_ADJACENCY_MATRIX& g, graph_bundle_t /*unused*/, GProp&& value) {
  g[graph_bundle] = std::forward<GProp>(value);
}

// Returns a reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto& get_property(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v, BAGL_ADJACENCY_MATRIX& g) {
  return g.get_property(v);
}

// Returns a const-reference to the vertex-property associated to the given vertex descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
const auto& get_property(typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v, const BAGL_ADJACENCY_MATRIX& g) {
  return g.get_property(v);
}

// Returns a reference to the edge-property associated to the given edge descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto& get_property(const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e, BAGL_ADJACENCY_MATRIX& g) {
  return g.get_property(e);
}

// Returns a const-reference to the edge-property associated to the given edge descriptor.
template <BAGL_ADJACENCY_MATRIX_PARAMS>
const auto& get_property(const typename BAGL_ADJACENCY_MATRIX::edge_descriptor& e, const BAGL_ADJACENCY_MATRIX& g) {
  return g.get_property(e);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto& get_property(BAGL_ADJACENCY_MATRIX& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
const auto& get_property(const BAGL_ADJACENCY_MATRIX& g, graph_all_t /*unused*/) {
  return g.get_property(graph_all);
}

// Handle graph property tags, also handles graph_bundle_t.

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>, lookup_one_property_t<GP, Tag>&>
get_property(BAGL_ADJACENCY_MATRIX& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Tag>
std::enable_if_t<std::is_same_v<property_kind_t<Tag>, graph_property_tag>, const lookup_one_property_t<GP, Tag>&>
get_property(const BAGL_ADJACENCY_MATRIX& g, Tag /*unused*/) {
  return get_property_value(g.get_property(graph_all), Tag{});
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename T, typename Bundle>
struct property_map<BAGL_ADJACENCY_MATRIX, T Bundle::*> {
  using non_const_Bundle = std::remove_cv_t<Bundle>;
  using non_const_T = std::remove_cv_t<T>;
  static constexpr bool is_vertex_bundle_v =
      std::is_convertible_v<typename BAGL_ADJACENCY_MATRIX::vertex_bundled*, non_const_Bundle*>;
  static constexpr bool is_edge_bundle_v =
      std::is_convertible_v<typename BAGL_ADJACENCY_MATRIX::edge_bundled*, non_const_Bundle*>;
  using tag_type = std::conditional_t<is_vertex_bundle_v, vertex_bundle_t,
                                      std::conditional_t<is_edge_bundle_v, edge_bundle_t, graph_bundle_t>>;
  using type = bundle_member_property_map<non_const_T, BAGL_ADJACENCY_MATRIX, tag_type>;
  using const_type = bundle_member_property_map<const non_const_T, const BAGL_ADJACENCY_MATRIX, tag_type>;
};

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename T, typename Bundle>
auto get(T Bundle::*p, BAGL_ADJACENCY_MATRIX& g) {
  return typename property_map<BAGL_ADJACENCY_MATRIX, T Bundle::*>::type{&g, p};
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename T, typename Bundle>
auto get(T Bundle::*p, const BAGL_ADJACENCY_MATRIX& g) {
  return typename property_map<BAGL_ADJACENCY_MATRIX, T Bundle::*>::const_type{&g, p};
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename T, typename Bundle, typename Key>
const std::remove_cv_t<T>& get(T Bundle::*p, const BAGL_ADJACENCY_MATRIX& g, const Key& k) {
  return (g[k]).*p;
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename T, typename Bundle, typename Key>
void put(T Bundle::*p, BAGL_ADJACENCY_MATRIX& g, const Key& k, T&& val) {
  (g[k]).*p = std::forward<T>(val);
}

/***********************************************************************************************
 *                             Property Maps (from tags)
 * ********************************************************************************************/

struct adj_matrix_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using value_type = typename property_value<Property, Tag>::type;

    using type = tagged_in_property_property_map<value_type, Graph, Tag>;
    using const_type = tagged_in_property_property_map<const value_type, const Graph, Tag>;
  };
};

/* specializations used by graph/properties.hpp */
template <>
struct vertex_property_selector<adjacency_matrix_tag> {
  using type = adj_matrix_property_selector;
};

template <>
struct edge_property_selector<adjacency_matrix_tag> {
  using type = adj_matrix_property_selector;
};

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Property>
auto get(Property p, BAGL_ADJACENCY_MATRIX& g) {
  using Map = typename property_map<BAGL_ADJACENCY_MATRIX, Property>::type;
  return Map(&g, p);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Property>
auto get(Property p, const BAGL_ADJACENCY_MATRIX& g) {
  using Map = typename property_map<BAGL_ADJACENCY_MATRIX, Property>::const_type;
  return Map(&g, p);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Property, typename Key>
decltype(auto) get(Property p, const BAGL_ADJACENCY_MATRIX& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_ADJACENCY_MATRIX& g, const Key& k) {
  return get_property_value(g.get_property(k), p);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS, typename Property, typename Key, typename Value>
void put(Property p, BAGL_ADJACENCY_MATRIX& g, const Key& k, Value&& val) {
  get_property_value(g.get_property(k), p) = std::forward<Value>(val);
}

//=========================================================================
// Vertex Index Property Map

template <BAGL_ADJACENCY_MATRIX_PARAMS>
struct property_map<BAGL_ADJACENCY_MATRIX, vertex_index_t> {
  using Vertex = typename BAGL_ADJACENCY_MATRIX::vertex_descriptor;
  using type = typed_identity_property_map<Vertex>;
  using const_type = type;
};

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(vertex_index_t /*unused*/, BAGL_ADJACENCY_MATRIX& /*unused*/) {
  return typed_identity_property_map<typename BAGL_ADJACENCY_MATRIX::vertex_descriptor>();
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(vertex_index_t /*unused*/, const BAGL_ADJACENCY_MATRIX& /*unused*/) {
  return typed_identity_property_map<typename BAGL_ADJACENCY_MATRIX::vertex_descriptor>();
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(vertex_index_t /*unused*/, BAGL_ADJACENCY_MATRIX& /*unused*/,
         typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v) {
  return v;
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(vertex_index_t /*unused*/, const BAGL_ADJACENCY_MATRIX& /*unused*/,
         typename BAGL_ADJACENCY_MATRIX::vertex_descriptor v) {
  return v;
}

//=========================================================================
// Edge Index Property Map

template <BAGL_ADJACENCY_MATRIX_PARAMS>
struct property_map<BAGL_ADJACENCY_MATRIX, edge_index_t> {
  using Edge = typename BAGL_ADJACENCY_MATRIX::edge_descriptor;
  using type = data_member_property_map<const std::size_t, const Edge>;
  using const_type = type;
};

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(edge_index_t /*unused*/, BAGL_ADJACENCY_MATRIX& /*unused*/) {
  using Edge = typename BAGL_ADJACENCY_MATRIX::edge_descriptor;
  return data_member_property_map<const std::size_t, const Edge>(&Edge::edge_index);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(edge_index_t /*unused*/, const BAGL_ADJACENCY_MATRIX& /*unused*/) {
  using Edge = typename BAGL_ADJACENCY_MATRIX::edge_descriptor;
  return data_member_property_map<const std::size_t, const Edge>(&Edge::edge_index);
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(edge_index_t /*unused*/, BAGL_ADJACENCY_MATRIX& /*unused*/,
         typename BAGL_ADJACENCY_MATRIX::edge_descriptor e) {
  return e.edge_index;
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto get(edge_index_t /*unused*/, const BAGL_ADJACENCY_MATRIX& /*unused*/,
         typename BAGL_ADJACENCY_MATRIX::edge_descriptor e) {
  return e.edge_index;
}

//=========================================================================
// Other Functions

template <BAGL_ADJACENCY_MATRIX_PARAMS>
auto vertex(std::size_t n, const BAGL_ADJACENCY_MATRIX& /*unused*/) {
  return n;
}

template <BAGL_ADJACENCY_MATRIX_PARAMS>
struct graph_mutability_traits<BAGL_ADJACENCY_MATRIX> {
  using category = mutable_edge_property_graph_tag;
};

#undef BAGL_ADJACENCY_MATRIX_PARAMS
#undef BAGL_ADJACENCY_MATRIX

}  // namespace bagl

#endif  // BAGL_BAGL_ADJACENCY_MATRIX_H_
