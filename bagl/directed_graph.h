// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DIRECTED_GRAPH_H_
#define BAGL_BAGL_DIRECTED_GRAPH_H_

#include <limits>
#include <ranges>
#include <type_traits>

#include "bagl/adjacency_list.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/transform_value_property_map.h"

namespace bagl {
struct directed_graph_tag {};

/**
 * The directed_graph class template is a simplified version of the BGL
 * adjacency list. This class is provided for ease of use, but may not
 * perform as well as custom-defined adjacency list classes. Instances of
 * this template model the BidirectionalGraph, VertexIndexGraph, and
 * EdgeIndexGraph concepts. The graph is also fully mutable, supporting
 * both insertions and removals of vertices and edges.
 *
 * @note Special care must be taken when removing vertices or edges since
 * those operations can invalidate the numbering of vertices.
 */
template <typename VertexProp = no_property, typename EdgeProp = no_property, typename GraphProp = no_property>
class directed_graph {
 public:
  using graph_property_type = GraphProp;
  using vertex_property_type = VertexProp;
  using edge_property_type = EdgeProp;
  using graph_bundled = lookup_one_property_t<GraphProp, graph_bundle_t>;
  using vertex_bundled = lookup_one_property_t<VertexProp, vertex_bundle_t>;
  using edge_bundled = lookup_one_property_t<EdgeProp, edge_bundle_t>;

  // Embed indices into the vertex type.
  using internal_vertex_property = property<vertex_index_t, std::size_t, vertex_property_type>;
  using internal_edge_property = property<edge_index_t, std::size_t, edge_property_type>;

  using graph_type =
      adjacency_list<list_s, list_s, bidirectional_s, internal_vertex_property, internal_edge_property, GraphProp>;

 private:
  // storage selectors
  using vertex_list_selector = typename graph_type::vertex_list_selector;
  using out_edge_list_selector = typename graph_type::out_edge_list_selector;
  using directed_selector = typename graph_type::directed_selector;

  // Hidden delegated-forwarding constructor.
  // So that we can transform some args in public constructors.
  struct non_empty_graph_ctor_t {};
  template <typename... Args>
  explicit directed_graph(non_empty_graph_ctor_t /*unused*/, Args&&... args)
      : graph_(std::forward<Args>(args)...),
        num_vertices_(num_vertices(graph_)),
        num_edges_(num_edges(graph_)),
        max_vertex_index_(num_vertices_),
        max_edge_index_(num_edges_) {
    renumber_indices();
  }

 public:
  // more commonly used graph types
  using vertices_size_type = typename graph_type::vertices_size_type;
  using edges_size_type = typename graph_type::edges_size_type;
  using degree_size_type = typename graph_type::degree_size_type;
  using vertex_descriptor = typename graph_type::vertex_descriptor;
  using edge_descriptor = typename graph_type::edge_descriptor;

  // miscellaneous types
  using graph_tag = directed_graph_tag;
  using directed_category = typename graph_type::directed_category;
  using edge_parallel_category = typename graph_type::edge_parallel_category;
  using traversal_category = typename graph_type::traversal_category;

  using vertex_index_type = std::size_t;
  using edge_index_type = std::size_t;

  explicit directed_graph(GraphProp const& p = GraphProp()) : graph_(p) {}
  ~directed_graph() = default;

  directed_graph(const directed_graph& x) = default;
  directed_graph(directed_graph&& x) noexcept = default;
  directed_graph& operator=(const directed_graph& x) = default;
  directed_graph& operator=(directed_graph&& x) noexcept = default;

  // Construct from a given number of vertices and an edge range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::size_t> directed_graph(
      vertices_size_type num_vertices, const EdgeRange& e_range, graph_property_type graph_prop = {})
      : directed_graph(non_empty_graph_ctor_t{}, num_vertices, e_range, std::move(graph_prop)) {}

  // Construct from a given number of vertices and an edge and edge-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::size_t> &&
      std::convertible_to<std::ranges::range_reference_t<EdgePropRange>, edge_property_type>
      directed_graph(vertices_size_type num_vertices, const EdgeRange& e_range, const EdgePropRange& ep_range,
                     graph_property_type graph_prop = {})
      : directed_graph(non_empty_graph_ctor_t{}, num_vertices, e_range,
                       ep_range | std::views::transform(top_property_value_adder<edge_index_t, std::size_t>{0}),
                       std::move(graph_prop)) {}

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::size_t> &&
      std::convertible_to<std::ranges::range_reference_t<VertexPropRange>, vertex_property_type>
      directed_graph(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                     graph_property_type graph_prop = {})
      : directed_graph(non_empty_graph_ctor_t{}, num_vertices,
                       vp_range | std::views::transform(top_property_value_adder<vertex_index_t, std::size_t>{0}),
                       e_range, std::move(graph_prop)) {}

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange,
            std::ranges::input_range EdgePropRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::size_t> &&
      std::convertible_to<std::ranges::range_reference_t<VertexPropRange>, vertex_property_type> &&
      std::convertible_to<std::ranges::range_reference_t<EdgePropRange>, edge_property_type>
      directed_graph(vertices_size_type num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                     const EdgePropRange& ep_range, graph_property_type graph_prop = {})
      : directed_graph(non_empty_graph_ctor_t{}, num_vertices,
                       vp_range | std::views::transform(top_property_value_adder<vertex_index_t, std::size_t>{0}),
                       e_range,
                       ep_range | std::views::transform(top_property_value_adder<edge_index_t, std::size_t>{0}),
                       std::move(graph_prop)) {}

  // The impl_() methods are not part of the public interface.
  graph_type& impl() { return graph_; }

  graph_type const& impl() const { return graph_; }

  // The following methods are not part of the public interface
  [[nodiscard]] vertices_size_type num_vertices() const { return num_vertices_; }

 private:
  // This helper function manages the attribution of vertex indices.
  vertex_descriptor make_index(vertex_descriptor v) {
    put(vertex_index, graph_, v, max_vertex_index_);
    num_vertices_++;
    max_vertex_index_++;
    return v;
  }

 public:
  vertex_descriptor add_vertex() { return make_index(add_vertex(graph_)); }

  template <typename VProp>
  vertex_descriptor add_vertex(VProp&& p) {
    return make_index(add_vertex(internal_vertex_property(0, std::forward<VProp>(p)), graph_));
  }

  void clear_vertex(vertex_descriptor v) {
    num_edges_ -= degree(v, graph_);
    clear_vertex(v, graph_);
  }

  void remove_vertex(vertex_descriptor v) {
    remove_vertex(v, graph_);
    --num_vertices_;
  }

  [[nodiscard]] edges_size_type num_edges() const { return num_edges_; }

 private:
  // A helper function for managing edge index attributes.
  std::pair<edge_descriptor, bool> const& make_index(std::pair<edge_descriptor, bool> const& x) {
    if (x.second) {
      put(edge_index, graph_, x.first, max_edge_index_);
      ++num_edges_;
      ++max_edge_index_;
    }
    return x;
  }

 public:
  std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u, vertex_descriptor v) {
    return make_index(add_edge(u, v, graph_));
  }

  template <typename EProp>
  std::pair<edge_descriptor, bool> add_edge(vertex_descriptor u, vertex_descriptor v, EProp&& p) {
    return make_index(add_edge(u, v, internal_edge_property(0, std::forward<EProp>(p)), graph_));
  }

  void remove_edge(vertex_descriptor u, vertex_descriptor v) {
    remove_out_edge_if(u, [this, v](edge_descriptor e) { return target(e, graph_) == v; });
  }

  void remove_edge(edge_descriptor e) {
    remove_edge(e, graph_);
    --num_edges_;
  }

  template <typename Predicate>
  void remove_edge_if(Predicate pred) {
    // find all edges matching predicate.
    auto e_rg = edges(graph_);
    for (auto e_i = e_rg.begin(); e_i != e_rg.end(); /*increment before remove*/) {
      auto e = *e_i++;
      if (pred(e)) {
        // This loop breaks of selections change.
        static_assert(std::is_same_v<vertex_list_selector, list_s> && std::is_same_v<out_edge_list_selector, list_s>);
        remove_edge(e);
      }
    }
  }

  template <typename Predicate>
  void remove_out_edge_if(vertex_descriptor u, Predicate pred) {
    // find all edges matching predicate.
    auto oe_rg = out_edges(u, graph_);
    for (auto oe_i = oe_rg.begin(); oe_i != oe_rg.end(); /*increment before remove*/) {
      auto e = *oe_i++;
      if (pred(e)) {
        // This loop breaks of selections change.
        static_assert(std::is_same_v<vertex_list_selector, list_s> && std::is_same_v<out_edge_list_selector, list_s>);
        remove_edge(e);
      }
    }
  }

  template <typename Predicate>
  void remove_in_edge_if(vertex_descriptor v, Predicate pred) {
    // find all edges matching predicate.
    auto ie_rg = in_edges(v, graph_);
    for (auto ie_i = ie_rg.begin(); ie_i != ie_rg.end(); /*increment before remove*/) {
      auto e = *ie_i++;
      if (pred(e)) {
        // This loop breaks of selections change.
        static_assert(std::is_same_v<vertex_list_selector, list_s> && std::is_same_v<out_edge_list_selector, list_s> &&
                      std::is_same_v<directed_selector, bidirectional_s>);
        remove_edge(e);
      }
    }
  }

  [[nodiscard]] vertex_index_type max_vertex_index() const { return max_vertex_index_; }

  void renumber_vertex_indices() {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    auto v_rg = vertices(graph_);
    max_vertex_index_ = renumber_vertex_indices(*v_rg.begin(), *v_rg.end(), 0);
  }

  void remove_vertex_and_renumber_indices(vertex_descriptor i) {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    vertex_descriptor j = std::next(i);
    auto end = vertices(graph_).end();
    vertex_index_type n = get(vertex_index, graph_, i);

    // remove the offending vertex and renumber everything after
    remove_vertex(i);
    max_vertex_index_ = renumber_vertex_indices(j, end, n);
  }

  [[nodiscard]] edge_index_type max_edge_index() const { return max_edge_index_; }

  void renumber_edge_indices() {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    auto v_rg = vertices(graph_);
    max_edge_index_ = renumber_edge_indices(*v_rg.begin(), *v_rg.end(), 0);
  }

  void remove_edge_and_renumber_indices(edge_descriptor e) {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    auto v_rg = vertices(graph_);
    auto i = source(e, graph_);
    auto j = i;
    vertex_index_type n = 0;
    if (e != *out_edges(i, graph_).begin()) {
      n = get(edge_index, graph_, *out_edges(i, graph_).begin());
    } else {
      // Find re-indexing start vertex before source of `e`.
      while (true) {
        if (j == *v_rg.begin()) {
          n = 0;
          break;
        }
        --j;
        auto oe_rg = out_edges(j, graph_);
        if (!oe_rg.empty()) {
          n = get(edge_index, graph_, *oe_rg.begin());
          break;
        }
      }
    }
    // remove the offending edge and renumber everything after
    remove_edge(e);
    max_edge_index_ = renumber_edge_indices(j, *v_rg.end(), n);
  }

  void renumber_indices() {
    renumber_vertex_indices();
    renumber_edge_indices();
  }

  // bundled property support
  vertex_bundled& operator[](vertex_descriptor v) { return graph_[v]; }
  vertex_bundled const& operator[](vertex_descriptor v) const { return graph_[v]; }

  edge_bundled& operator[](edge_descriptor e) { return graph_[e]; }
  edge_bundled const& operator[](edge_descriptor e) const { return graph_[e]; }

  graph_bundled& operator[](graph_bundle_t /*unused*/) { return graph_[graph_bundle]; }
  graph_bundled const& operator[](graph_bundle_t /*unused*/) const { return graph_[graph_bundle]; }

  // Getting full properties out (mostly to extract tagged properties).
  vertex_property_type& get_property(vertex_descriptor v) { return remove_top_property(graph_.get_property(v)); }
  vertex_property_type const& get_property(vertex_descriptor v) const {
    return remove_top_property(graph_.get_property(v));
  }

  edge_property_type& get_property(edge_descriptor e) { return remove_top_property(graph_.get_property(e)); }
  edge_property_type const& get_property(edge_descriptor e) const {
    return remove_top_property(graph_.get_property(e));
  }

  graph_property_type& get_property(graph_all_t /*unused*/) { return graph_.get_property(graph_all); }
  graph_property_type const& get_property(graph_all_t /*unused*/) const { return graph_.get_property(graph_all); }

  // Graph concepts
  static vertex_descriptor null_vertex() { return graph_type::null_vertex(); }

  void clear() {
    graph_.clear();
    num_vertices_ = max_vertex_index_ = 0;
    num_edges_ = max_edge_index_ = 0;
  }

  void swap(directed_graph& g) {
    graph_.swap(g.graph_);
    std::swap(num_vertices_, g.num_vertices_);
    std::swap(max_vertex_index_, g.max_vertex_index_);
    std::swap(num_edges_, g.num_edges_);
    std::swap(max_edge_index_, g.max_edge_index_);
  }

 private:
  vertices_size_type renumber_vertex_indices(vertex_descriptor i, vertex_descriptor end, vertices_size_type n) {
    auto indices = get(vertex_index, graph_);
    for (; i != end; ++i) {
      indices[i] = n++;
    }
    return n;
  }

  edges_size_type renumber_edge_indices(vertex_descriptor i, vertex_descriptor end, edges_size_type n) {
    auto indices = get(edge_index, graph_);
    for (; i != end; ++i) {
      for (edge_descriptor e : out_edges(i, graph_)) {
        indices[e] = n++;
      }
    }
    return n;
  }

  graph_type graph_;
  vertices_size_type num_vertices_ = 0;
  edges_size_type num_edges_ = 0;
  vertex_index_type max_vertex_index_ = 0;
  edge_index_type max_edge_index_ = 0;
};

#define DIRECTED_GRAPH_PARAMS typename VP, typename EP, typename GP
#define DIRECTED_GRAPH directed_graph<VP, EP, GP>

// IncidenceGraph concepts
template <DIRECTED_GRAPH_PARAMS>
auto source(typename DIRECTED_GRAPH::edge_descriptor e, DIRECTED_GRAPH const& g) {
  return source(e, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto target(typename DIRECTED_GRAPH::edge_descriptor e, DIRECTED_GRAPH const& g) {
  return target(e, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto out_degree(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return out_degree(v, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto out_edges(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return out_edges(v, g.impl());
}

// BidirectionalGraph concepts
template <DIRECTED_GRAPH_PARAMS>
auto in_degree(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return in_degree(v, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto in_edges(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return in_edges(v, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto degree(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return degree(v, g.impl());
}

// AdjacencyGraph concepts
template <DIRECTED_GRAPH_PARAMS>
auto adjacent_vertices(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return adjacent_vertices(v, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto vertex(typename DIRECTED_GRAPH::vertices_size_type n, DIRECTED_GRAPH const& g) {
  return vertex(n, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto edge(typename DIRECTED_GRAPH::vertex_descriptor u, typename DIRECTED_GRAPH::vertex_descriptor v,
          DIRECTED_GRAPH const& g) {
  return edge(u, v, g.impl());
}

// VertexListGraph concepts
template <DIRECTED_GRAPH_PARAMS>
auto num_vertices(DIRECTED_GRAPH const& g) {
  return g.num_vertices();
}

template <DIRECTED_GRAPH_PARAMS>
auto vertices(DIRECTED_GRAPH const& g) {
  return vertices(g.impl());
}

// EdgeListGraph concepts
template <DIRECTED_GRAPH_PARAMS>
auto num_edges(DIRECTED_GRAPH const& g) {
  return g.num_edges();
}

template <DIRECTED_GRAPH_PARAMS>
auto edges(DIRECTED_GRAPH const& g) {
  return edges(g.impl());
}

// MutableGraph concepts
template <DIRECTED_GRAPH_PARAMS>
auto add_vertex(DIRECTED_GRAPH& g) {
  return g.add_vertex();
}

template <DIRECTED_GRAPH_PARAMS, typename VProp>
auto add_vertex(VProp&& p, DIRECTED_GRAPH& g) {
  return g.add_vertex(std::forward<VProp>(p));
}

template <DIRECTED_GRAPH_PARAMS>
void clear_vertex(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH& g) {
  return g.clear_vertex(v);
}

template <DIRECTED_GRAPH_PARAMS>
void remove_vertex(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH& g) {
  return g.remove_vertex(v);
}

template <DIRECTED_GRAPH_PARAMS>
auto add_edge(typename DIRECTED_GRAPH::vertex_descriptor u, typename DIRECTED_GRAPH::vertex_descriptor v,
              DIRECTED_GRAPH& g) {
  return g.add_edge(u, v);
}

template <DIRECTED_GRAPH_PARAMS, typename EProp>
auto add_edge(typename DIRECTED_GRAPH::vertex_descriptor u, typename DIRECTED_GRAPH::vertex_descriptor v, EProp&& p,
              DIRECTED_GRAPH& g) {
  return g.add_edge(u, v, std::forward<EProp>(p));
}

template <DIRECTED_GRAPH_PARAMS>
void remove_edge(typename DIRECTED_GRAPH::vertex_descriptor u, typename DIRECTED_GRAPH::vertex_descriptor v,
                 DIRECTED_GRAPH& g) {
  return g.remove_edge(u, v);
}

template <DIRECTED_GRAPH_PARAMS, typename EProp>
void remove_edge(typename DIRECTED_GRAPH::vertex_descriptor u, typename DIRECTED_GRAPH::vertex_descriptor v, EProp& ep,
                 DIRECTED_GRAPH& g) {
  auto [e, e_found] = edge(u, v, g);
  if (e_found) {
    // We can only get property of first edge found.
    ep = std::move(g.get_property(e));
    // Remove all (u,v) edges.
    g.remove_edge(u, v);
  }
}

template <DIRECTED_GRAPH_PARAMS>
void remove_edge(typename DIRECTED_GRAPH::edge_descriptor e, DIRECTED_GRAPH& g) {
  return g.remove_edge(e);
}

template <DIRECTED_GRAPH_PARAMS, typename EProp>
void remove_edge(typename DIRECTED_GRAPH::edge_descriptor e, EProp& ep, DIRECTED_GRAPH& g) {
  ep = std::move(g.get_property(e));
  g.remove_edge(e);
}

template <DIRECTED_GRAPH_PARAMS, class Predicate>
void remove_edge_if(Predicate pred, DIRECTED_GRAPH& g) {
  return g.remove_edge_if(pred);
}

template <DIRECTED_GRAPH_PARAMS, class Predicate>
void remove_out_edge_if(typename DIRECTED_GRAPH::vertex_descriptor v, Predicate pred, DIRECTED_GRAPH& g) {
  return g.remove_out_edge_if(v, pred);
}

template <DIRECTED_GRAPH_PARAMS, class Predicate>
void remove_in_edge_if(typename DIRECTED_GRAPH::vertex_descriptor v, Predicate pred, DIRECTED_GRAPH& g) {
  return g.remove_in_edge_if(v, pred);
}

template <DIRECTED_GRAPH_PARAMS, typename Property>
struct property_map<DIRECTED_GRAPH, Property> : property_map<typename DIRECTED_GRAPH::graph_type, Property> {};

template <DIRECTED_GRAPH_PARAMS>
struct property_map<DIRECTED_GRAPH, vertex_all_t> {
  using const_type =
      transform_value_property_map<top_property_remover,
                                   property_map_const_t<typename DIRECTED_GRAPH::graph_type, vertex_all_t>>;
  using type = transform_value_property_map<top_property_remover,
                                            property_map_t<typename DIRECTED_GRAPH::graph_type, vertex_all_t>>;
};

template <DIRECTED_GRAPH_PARAMS>
struct property_map<DIRECTED_GRAPH, edge_all_t> {
  using const_type =
      transform_value_property_map<top_property_remover,
                                   property_map_const_t<typename DIRECTED_GRAPH::graph_type, edge_all_t>>;
  using type = transform_value_property_map<top_property_remover,
                                            property_map_t<typename DIRECTED_GRAPH::graph_type, edge_all_t>>;
};

// PropertyGraph concepts
template <DIRECTED_GRAPH_PARAMS, typename Property>
auto get(Property p, DIRECTED_GRAPH& g) {
  return get(p, g.impl());
}

template <DIRECTED_GRAPH_PARAMS, typename Property>
auto get(Property p, DIRECTED_GRAPH const& g) {
  return get(p, g.impl());
}

template <DIRECTED_GRAPH_PARAMS>
auto get(vertex_all_t /*unused*/, DIRECTED_GRAPH& g) {
  return property_map_t<DIRECTED_GRAPH, vertex_all_t>(top_property_remover{}, get(vertex_all, g.impl()));
}

template <DIRECTED_GRAPH_PARAMS>
auto get(vertex_all_t /*unused*/, DIRECTED_GRAPH const& g) {
  return property_map_const_t<DIRECTED_GRAPH, vertex_all_t>(top_property_remover{}, get(vertex_all, g.impl()));
}

template <DIRECTED_GRAPH_PARAMS>
auto get(edge_all_t /*unused*/, DIRECTED_GRAPH& g) {
  return property_map_t<DIRECTED_GRAPH, edge_all_t>(top_property_remover{}, get(edge_all, g.impl()));
}

template <DIRECTED_GRAPH_PARAMS>
auto get(edge_all_t /*unused*/, DIRECTED_GRAPH const& g) {
  return property_map_const_t<DIRECTED_GRAPH, edge_all_t>(top_property_remover{}, get(edge_all, g.impl()));
}

template <DIRECTED_GRAPH_PARAMS, typename Property, typename Key>
decltype(auto) get(Property p, DIRECTED_GRAPH const& g, Key const& k) {
  return get(p, g.impl(), k);
}
template <DIRECTED_GRAPH_PARAMS, typename Property, typename Key>
decltype(auto) get(Property p, DIRECTED_GRAPH& g, Key const& k) {
  return get(p, g.impl(), k);
}

template <DIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(vertex_all_t /*unused*/, DIRECTED_GRAPH const& g, Key const& k) {
  return remove_top_property(get(vertex_all, g.impl(), k));
}
template <DIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(vertex_all_t /*unused*/, DIRECTED_GRAPH& g, Key const& k) {
  return remove_top_property(get(vertex_all, g.impl(), k));
}

template <DIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(edge_all_t /*unused*/, DIRECTED_GRAPH const& g, Key const& k) {
  return remove_top_property(get(edge_all, g.impl(), k));
}
template <DIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(edge_all_t /*unused*/, DIRECTED_GRAPH& g, Key const& k) {
  return remove_top_property(get(edge_all, g.impl(), k));
}

template <DIRECTED_GRAPH_PARAMS, typename Property, typename Key, typename Value>
void put(Property p, DIRECTED_GRAPH& g, Key const& k, Value&& v) {
  put(p, g.impl(), k, std::forward<Value>(v));
}

template <DIRECTED_GRAPH_PARAMS, typename Key, typename Value>
void put(vertex_all_t /*unused*/, DIRECTED_GRAPH& g, Key const& k, Value&& v) {
  put(vertex_all, g.impl(), k, add_top_property(vertex_index, get(vertex_index, g.impl(), k), std::forward<Value>(v)));
}

template <DIRECTED_GRAPH_PARAMS, typename Key, typename Value>
void put(edge_all_t /*unused*/, DIRECTED_GRAPH& g, Key const& k, Value&& v) {
  put(edge_all, g.impl(), k, add_top_property(edge_index, get(edge_index, g.impl(), k), std::forward<Value>(v)));
}

template <DIRECTED_GRAPH_PARAMS, class Property>
auto& get_property(DIRECTED_GRAPH& g, Property p) {
  return get_property(g.impl(), p);
}

template <DIRECTED_GRAPH_PARAMS, class Property>
auto const& get_property(DIRECTED_GRAPH const& g, Property p) {
  return get_property(g.impl(), p);
}

template <DIRECTED_GRAPH_PARAMS, class Property, class Value>
void set_property(DIRECTED_GRAPH& g, Property p, Value&& v) {
  return set_property(g.impl(), p, std::forward<Value>(v));
}

// Vertex index management

template <DIRECTED_GRAPH_PARAMS>
auto get_vertex_index(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH const& g) {
  return get(vertex_index, g, v);
}

template <DIRECTED_GRAPH_PARAMS>
auto max_vertex_index(DIRECTED_GRAPH const& g) {
  return g.max_vertex_index();
}

template <DIRECTED_GRAPH_PARAMS>
void renumber_vertex_indices(DIRECTED_GRAPH& g) {
  g.renumber_vertex_indices();
}

template <DIRECTED_GRAPH_PARAMS>
void remove_vertex_and_renumber_indices(typename DIRECTED_GRAPH::vertex_descriptor v, DIRECTED_GRAPH& g) {
  g.remove_vertex_and_renumber_indices(v);
}

// Edge index management
template <DIRECTED_GRAPH_PARAMS>
auto get_edge_index(typename DIRECTED_GRAPH::edge_descriptor v, DIRECTED_GRAPH const& g) {
  return get(edge_index, g, v);
}

template <DIRECTED_GRAPH_PARAMS>
auto max_edge_index(DIRECTED_GRAPH const& g) {
  return g.max_edge_index();
}

template <DIRECTED_GRAPH_PARAMS>
void renumber_edge_indices(DIRECTED_GRAPH& g) {
  g.renumber_edge_indices();
}

template <DIRECTED_GRAPH_PARAMS>
void remove_edge_and_renumber_indices(typename DIRECTED_GRAPH::edge_descriptor e, DIRECTED_GRAPH& g) {
  g.remove_edge_and_renumber_indices(e);
}

// Index management
template <DIRECTED_GRAPH_PARAMS>
void renumber_indices(DIRECTED_GRAPH& g) {
  g.renumber_indices();
}

// Mutability Traits
template <DIRECTED_GRAPH_PARAMS>
struct graph_mutability_traits<DIRECTED_GRAPH> {
  using category = mutable_property_graph_tag;
};

#undef DIRECTED_GRAPH_PARAMS
#undef DIRECTED_GRAPH

}  // namespace bagl

#endif  // BAGL_BAGL_DIRECTED_GRAPH_H_
