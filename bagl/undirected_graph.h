// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_UNDIRECTED_GRAPH_H_
#define BAGL_BAGL_UNDIRECTED_GRAPH_H_

#include <limits>
#include <ranges>
#include <type_traits>

#include "bagl/adjacency_list.h"
#include "bagl/graph_selectors.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"

namespace bagl {
struct undirected_graph_tag {};

/**
 * The undirected_graph class template is a simplified version of the BGL
 * adjacency list. This class is provided for ease of use, but may not
 * perform as well as custom-defined adjacency list classes. Instances of
 * this template model the VertexIndexGraph, and EdgeIndexGraph concepts. The
 * graph is also fully mutable, supporting both insertions and removals of
 * vertices and edges.
 *
 * @note Special care must be taken when removing vertices or edges since
 * those operations can invalidate the numbering of vertices.
 */
template <typename VertexProp = no_property, typename EdgeProp = no_property, typename GraphProp = no_property>
class undirected_graph {
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

  // storage selectors
  using vertex_list_selector = list_s;
  using out_edge_list_selector = list_s;
  using directed_selector = undirected_s;

  using graph_type = adjacency_list<out_edge_list_selector, vertex_list_selector, directed_selector,
                                    internal_vertex_property, internal_edge_property, GraphProp>;

 private:
  // Hidden delegated-forwarding constructor.
  // So that we can transform some args in public constructors.
  struct non_empty_graph_ctor_t {};
  template <typename... Args>
  explicit undirected_graph(non_empty_graph_ctor_t /*unused*/, Args&&... args)
      : graph_(std::forward<Args>(args)...),
        max_vertex_index_(num_vertices(graph_)),
        max_edge_index_(num_edges(graph_)) {
    renumber_indices();
  }

 public:
  // more commonly used graph types
  using vertex_descriptor = typename graph_type::vertex_descriptor;
  using edge_descriptor = typename graph_type::edge_descriptor;

  // miscellaneous types
  using graph_tag = undirected_graph_tag;
  using directed_category = typename graph_type::directed_category;
  using edge_parallel_category = typename graph_type::edge_parallel_category;
  using traversal_category = typename graph_type::traversal_category;

  using vertex_index_type = std::size_t;
  using edge_index_type = std::size_t;

  undirected_graph() = default;
  explicit undirected_graph(const GraphProp& p) : graph_(0, p) {}
  explicit undirected_graph(std::size_t n, const GraphProp& p = {})
      : undirected_graph(non_empty_graph_ctor_t{}, n, p) {}
  ~undirected_graph() = default;

  undirected_graph(const undirected_graph& x) = default;
  undirected_graph(undirected_graph&& x) noexcept = default;
  undirected_graph& operator=(const undirected_graph& x) = default;
  undirected_graph& operator=(undirected_graph&& x) noexcept = default;

  // Construct from a given number of vertices and an edge range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>>
  undirected_graph(std::size_t num_vertices, const EdgeRange& e_range, graph_property_type graph_prop = {})
      : undirected_graph(non_empty_graph_ctor_t{}, num_vertices, e_range, std::move(graph_prop)) {}

  // Construct from a given number of vertices and an edge and edge-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range EdgeRange, std::ranges::input_range EdgePropRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>> &&
      std::convertible_to<std::ranges::range_reference_t<EdgePropRange>, edge_property_type>
      undirected_graph(std::size_t num_vertices, const EdgeRange& e_range, const EdgePropRange& ep_range,
                       graph_property_type graph_prop = {})
      : undirected_graph(non_empty_graph_ctor_t{}, num_vertices, e_range,
                         ep_range | std::views::transform(top_property_value_adder<edge_index_t, std::size_t>{0}),
                         std::move(graph_prop)) {}

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>> &&
      std::convertible_to<std::ranges::range_reference_t<VertexPropRange>, vertex_property_type>
      undirected_graph(std::size_t num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                       graph_property_type graph_prop = {})
      : undirected_graph(non_empty_graph_ctor_t{}, num_vertices,
                         vp_range | std::views::transform(top_property_value_adder<vertex_index_t, std::size_t>{0}),
                         e_range, std::move(graph_prop)) {}

  // Construct from a given number of vertices and an edge and vertex-property range.
  // Edges should be represented as pairs of vertex indices.
  template <std::ranges::input_range VertexPropRange, std::ranges::input_range EdgeRange,
            std::ranges::input_range EdgePropRange>
  requires std::convertible_to<std::ranges::range_value_t<EdgeRange>, std::pair<std::size_t, std::size_t>> &&
      std::convertible_to<std::ranges::range_reference_t<VertexPropRange>, vertex_property_type> &&
      std::convertible_to<std::ranges::range_reference_t<EdgePropRange>, edge_property_type>
      undirected_graph(std::size_t num_vertices, const VertexPropRange& vp_range, const EdgeRange& e_range,
                       const EdgePropRange& ep_range, graph_property_type graph_prop = {})
      : undirected_graph(non_empty_graph_ctor_t{}, num_vertices,
                         vp_range | std::views::transform(top_property_value_adder<vertex_index_t, std::size_t>{0}),
                         e_range,
                         ep_range | std::views::transform(top_property_value_adder<edge_index_t, std::size_t>{0}),
                         std::move(graph_prop)) {}

  // The impl_() methods are not part of the public interface.
  graph_type& impl() { return graph_; }

  graph_type const& impl() const { return graph_; }

 private:
  // This helper function manages the attribution of vertex indices.
  vertex_descriptor make_index(vertex_descriptor v) {
    put(vertex_index, graph_, v, max_vertex_index_);
    max_vertex_index_++;
    return v;
  }

 public:
  template <typename... VPArgs>
  vertex_descriptor add_vertex_and_update_index(VPArgs&&... vp_args) {
    return make_index(add_vertex(graph_, 0, std::forward<VPArgs>(vp_args)...));
  }

 private:
  // A helper function for managing edge index attributes.
  std::pair<edge_descriptor, bool> make_index(std::pair<edge_descriptor, bool> x) {
    if (x.second) {
      put(edge_index, graph_, x.first, max_edge_index_);
      ++max_edge_index_;
    }
    return x;
  }

 public:
  template <typename... EPArgs>
  std::pair<edge_descriptor, bool> add_edge_and_update_index(vertex_descriptor u, vertex_descriptor v,
                                                             EPArgs&&... ep_args) {
    return make_index(add_edge(u, v, graph_, 0, std::forward<EPArgs>(ep_args)...));
  }

  [[nodiscard]] vertex_index_type max_vertex_index() const { return max_vertex_index_; }

  void renumber_vertex_indices() {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    auto v_rg = bagl::vertices(graph_);
    max_vertex_index_ = renumber_vertex_indices(*v_rg.begin(), *v_rg.end(), 0);
  }

  void remove_vertex_and_renumber_indices(vertex_descriptor i) {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    vertex_descriptor j = std::next(i);
    auto end = bagl::vertices(graph_).end();
    vertex_index_type n = get(vertex_index, graph_, i);

    // remove the offending vertex and renumber everything after
    remove_vertex(i, graph_);
    max_vertex_index_ = renumber_vertex_indices(j, end, n);
  }

  [[nodiscard]] edge_index_type max_edge_index() const { return max_edge_index_; }

  void renumber_edge_indices() {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    auto v_rg = bagl::vertices(graph_);
    max_edge_index_ = renumber_edge_indices(*v_rg.begin(), *v_rg.end(), 0);
  }

  void remove_edge_and_renumber_indices(edge_descriptor e) {
    // With list_s we expect / know that vertex descriptors are actually list iterators.
    static_assert(std::input_iterator<vertex_descriptor>);
    auto v_rg = bagl::vertices(graph_);
    auto i = bagl::source(e, graph_);
    auto j = i;
    vertex_index_type n = 0;
    if (e != *bagl::out_edges(i, graph_).begin()) {
      n = get(edge_index, graph_, *bagl::out_edges(i, graph_).begin());
    } else {
      // Find re-indexing start vertex before source of `e`.
      while (true) {
        if (j == *v_rg.begin()) {
          n = 0;
          break;
        }
        --j;
        auto oe_rg = bagl::out_edges(j, graph_);
        if (!oe_rg.empty()) {
          n = get(edge_index, graph_, *oe_rg.begin());
          break;
        }
      }
    }
    // remove the offending edge and renumber everything after
    remove_edge(e, graph_);
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
    max_vertex_index_ = 0;
    max_edge_index_ = 0;
  }

  void swap(undirected_graph& g) {
    graph_.swap(g.graph_);
    std::swap(max_vertex_index_, g.max_vertex_index_);
    std::swap(max_edge_index_, g.max_edge_index_);
  }

 private:
  std::size_t renumber_vertex_indices(vertex_descriptor i, vertex_descriptor end, std::size_t n) {
    auto indices = get(vertex_index, graph_);
    for (; i != end; ++i) {
      indices[i] = n++;
    }
    return n;
  }

  std::size_t renumber_edge_indices(vertex_descriptor i, vertex_descriptor end, std::size_t n) {
    auto indices = get(edge_index, graph_);
    for (; i != end; ++i) {
      for (edge_descriptor e : bagl::out_edges(i, graph_)) {
        indices[e] = n++;
      }
    }
    return n;
  }

  graph_type graph_;
  vertex_index_type max_vertex_index_ = 0;
  edge_index_type max_edge_index_ = 0;
};

#define BAGL_UNDIRECTED_GRAPH_PARAMS typename VP, typename EP, typename GP
#define BAGL_UNDIRECTED_GRAPH undirected_graph<VP, EP, GP>

// IncidenceGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto source(typename BAGL_UNDIRECTED_GRAPH::edge_descriptor e, BAGL_UNDIRECTED_GRAPH const& g) {
  return source(e, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto target(typename BAGL_UNDIRECTED_GRAPH::edge_descriptor e, BAGL_UNDIRECTED_GRAPH const& g) {
  return target(e, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto out_degree(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return out_degree(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto out_edges(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return out_edges(v, g.impl());
}

// BidirectionalGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto in_degree(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return in_degree(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto in_edges(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return in_edges(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto incident_edges(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return out_edges(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto degree(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return degree(v, g.impl());
}

// AdjacencyGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto adjacent_vertices(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return adjacent_vertices(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto vertex(std::size_t n, BAGL_UNDIRECTED_GRAPH const& g) {
  return vertex(n, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto edge(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor u, typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v,
          BAGL_UNDIRECTED_GRAPH const& g) {
  return edge(u, v, g.impl());
}

// VertexListGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto num_vertices(BAGL_UNDIRECTED_GRAPH const& g) {
  return num_vertices(g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto vertices(BAGL_UNDIRECTED_GRAPH const& g) {
  return vertices(g.impl());
}

// EdgeListGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto num_edges(BAGL_UNDIRECTED_GRAPH const& g) {
  return num_edges(g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto edges(BAGL_UNDIRECTED_GRAPH const& g) {
  return edges(g.impl());
}

// MutableGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename... VPArgs>
auto add_vertex(BAGL_UNDIRECTED_GRAPH& g, VPArgs&&... vp_args) {
  return g.add_vertex_and_update_index(std::forward<VPArgs>(vp_args)...);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void clear_vertex(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH& g) {
  clear_vertex(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void remove_vertex(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH& g) {
  remove_vertex(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename VProp>
void remove_vertex(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH& g, VProp* vp) {
  if (vp != nullptr) {
    *vp = std::move(g.get_property(v));
  }
  remove_vertex(v, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename... EPArgs>
auto add_edge(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor u, typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v,
              BAGL_UNDIRECTED_GRAPH& g, EPArgs&&... ep_args) {
  return g.add_edge_and_update_index(u, v, std::forward<EPArgs>(ep_args)...);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void remove_edge(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor u,
                 typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH& g) {
  remove_out_edge_if(
      u, [&g, v](typename BAGL_UNDIRECTED_GRAPH::edge_descriptor e) { return target(e, g.impl()) == v; }, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename EProp>
void remove_edge(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor u,
                 typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH& g, EProp* ep) {
  auto [e, e_found] = edge(u, v, g);
  if (e_found) {
    // We can only get property of first edge found.
    if (ep != nullptr) {
      *ep = std::move(g.get_property(e));
    }
    // Remove all (u,v) edges.
    remove_edge(u, v, g);
  }
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void remove_edge(typename BAGL_UNDIRECTED_GRAPH::edge_descriptor e, BAGL_UNDIRECTED_GRAPH& g) {
  remove_edge(e, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename EProp>
void remove_edge(typename BAGL_UNDIRECTED_GRAPH::edge_descriptor e, BAGL_UNDIRECTED_GRAPH& g, EProp* ep) {
  if (ep != nullptr) {
    *ep = std::move(g.get_property(e));
  }
  remove_edge(e, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Predicate>
void remove_edge_if(Predicate pred, BAGL_UNDIRECTED_GRAPH& g) {
  remove_edge_if(pred, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Predicate>
void remove_incident_edge_if(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, Predicate pred,
                             BAGL_UNDIRECTED_GRAPH& g) {
  remove_out_edge_if(v, pred, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Predicate>
void remove_out_edge_if(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, Predicate pred, BAGL_UNDIRECTED_GRAPH& g) {
  remove_out_edge_if(v, pred, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Predicate>
void remove_in_edge_if(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, Predicate pred, BAGL_UNDIRECTED_GRAPH& g) {
  remove_in_edge_if(v, pred, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Property>
struct property_map<BAGL_UNDIRECTED_GRAPH, Property>
    : property_map<typename BAGL_UNDIRECTED_GRAPH::graph_type, Property> {};

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
struct property_map<BAGL_UNDIRECTED_GRAPH, vertex_all_t> {
  using const_type =
      transform_property_map<top_property_remover,
                             property_map_const_t<typename undirected_graph<VP, EP, GP>::graph_type, vertex_all_t>>;
  using type = transform_property_map<top_property_remover,
                                      property_map_t<typename undirected_graph<VP, EP, GP>::graph_type, vertex_all_t>>;
};

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
struct property_map<BAGL_UNDIRECTED_GRAPH, edge_all_t> {
  using const_type =
      transform_property_map<top_property_remover,
                             property_map_const_t<typename undirected_graph<VP, EP, GP>::graph_type, edge_all_t>>;
  using type = transform_property_map<top_property_remover,
                                      property_map_t<typename undirected_graph<VP, EP, GP>::graph_type, edge_all_t>>;
};

// PropertyGraph concepts
template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Property>
auto get(Property p, BAGL_UNDIRECTED_GRAPH& g) {
  return get(p, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Property>
auto get(Property p, BAGL_UNDIRECTED_GRAPH const& g) {
  return get(p, g.impl());
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto get(vertex_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH& g) {
  return property_map_t<BAGL_UNDIRECTED_GRAPH, vertex_all_t>(top_property_remover{}, get(vertex_all, g.impl()));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto get(vertex_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH const& g) {
  return property_map_const_t<BAGL_UNDIRECTED_GRAPH, vertex_all_t>(top_property_remover{}, get(vertex_all, g.impl()));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto get(edge_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH& g) {
  return property_map_t<BAGL_UNDIRECTED_GRAPH, edge_all_t>(top_property_remover{}, get(edge_all, g.impl()));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto get(edge_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH const& g) {
  return property_map_const_t<BAGL_UNDIRECTED_GRAPH, edge_all_t>(top_property_remover{}, get(edge_all, g.impl()));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_UNDIRECTED_GRAPH const& g, Key const& k) {
  return get(p, g.impl(), k);
}
template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Property, typename Key>
decltype(auto) get(Property p, BAGL_UNDIRECTED_GRAPH& g, Key const& k) {
  return get(p, g.impl(), k);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(vertex_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH const& g, Key const& k) {
  return remove_top_property(get(vertex_all, g.impl(), k));
}
template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(vertex_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH& g, Key const& k) {
  return remove_top_property(get(vertex_all, g.impl(), k));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(edge_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH const& g, Key const& k) {
  return remove_top_property(get(edge_all, g.impl(), k));
}
template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Key>
decltype(auto) get(edge_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH& g, Key const& k) {
  return remove_top_property(get(edge_all, g.impl(), k));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Property, typename Key, typename Value>
void put(Property p, BAGL_UNDIRECTED_GRAPH& g, Key const& k, Value&& v) {
  put(p, g.impl(), k, std::forward<Value>(v));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Key, typename Value>
void put(vertex_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH& g, Key const& k, Value&& v) {
  put(vertex_all, g.impl(), k, add_top_property(vertex_index, get(vertex_index, g.impl(), k), std::forward<Value>(v)));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, typename Key, typename Value>
void put(edge_all_t /*unused*/, BAGL_UNDIRECTED_GRAPH& g, Key const& k, Value&& v) {
  put(edge_all, g.impl(), k, add_top_property(edge_index, get(edge_index, g.impl(), k), std::forward<Value>(v)));
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, class Property>
auto& get_property(BAGL_UNDIRECTED_GRAPH& g, Property p) {
  return get_property(g.impl(), p);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, class Property>
auto const& get_property(BAGL_UNDIRECTED_GRAPH const& g, Property p) {
  return get_property(g.impl(), p);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS, class Property, class Value>
void set_property(BAGL_UNDIRECTED_GRAPH& g, Property p, Value&& v) {
  return set_property(g.impl(), p, std::forward<Value>(v));
}

// Indexed Vertex graph

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto get_vertex_index(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return get(vertex_index, g, v);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto max_vertex_index(BAGL_UNDIRECTED_GRAPH const& g) {
  return g.max_vertex_index();
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void renumber_vertex_indices(BAGL_UNDIRECTED_GRAPH& g) {
  g.renumber_vertex_indices();
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void remove_vertex_and_renumber_indices(typename BAGL_UNDIRECTED_GRAPH::vertex_descriptor v, BAGL_UNDIRECTED_GRAPH& g) {
  g.remove_vertex_and_renumber_indices(v);
}

// Edge index management

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto get_edge_index(typename BAGL_UNDIRECTED_GRAPH::edge_descriptor v, BAGL_UNDIRECTED_GRAPH const& g) {
  return get(edge_index, g, v);
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
auto max_edge_index(BAGL_UNDIRECTED_GRAPH const& g) {
  return g.max_edge_index();
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void renumber_edge_indices(BAGL_UNDIRECTED_GRAPH& g) {
  g.renumber_edge_indices();
}

template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void remove_edge_and_renumber_indices(typename BAGL_UNDIRECTED_GRAPH::edge_descriptor e, BAGL_UNDIRECTED_GRAPH& g) {
  g.remove_edge_and_renumber_indices(e);
}

// Index management
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
void renumber_indices(BAGL_UNDIRECTED_GRAPH& g) {
  g.renumber_indices();
}

// Mutability Traits
template <BAGL_UNDIRECTED_GRAPH_PARAMS>
struct graph_mutability_traits<BAGL_UNDIRECTED_GRAPH> {
  using category = mutable_property_graph_tag;
};

#undef BAGL_UNDIRECTED_GRAPH_PARAMS
#undef BAGL_UNDIRECTED_GRAPH

}  // namespace bagl

#endif  // BAGL_BAGL_UNDIRECTED_GRAPH_H_
