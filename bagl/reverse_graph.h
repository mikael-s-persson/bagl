// Copyright David Abrahams 2000.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_REVERSE_GRAPH_H_
#define BAGL_BAGL_REVERSE_GRAPH_H_

#include <ranges>
#include <tuple>
#include <type_traits>

#include "bagl/adjacency_range.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

struct reverse_graph_tag {};

namespace reverse_graph_detail {

template <typename EdgeDesc>
class reverse_graph_edge_descriptor {
 public:
  EdgeDesc underlying_descx_;  // Odd name is because this needs to be public but
                               // shouldn't be exposed to users anymore

  explicit reverse_graph_edge_descriptor(const EdgeDesc& underlying_descx = EdgeDesc())
      : underlying_descx_(underlying_descx) {}

  friend auto operator<=>(const reverse_graph_edge_descriptor& a, const reverse_graph_edge_descriptor& b) {
    return a.underlying_descx_ <=> b.underlying_descx_;
  }
};

template <typename EdgeRange>
auto reverse_edge_range(EdgeRange&& e_rg) {
  return std::forward<EdgeRange>(e_rg) | std::views::transform([](auto e) { return reverse_graph_edge_descriptor{e}; });
}

// Get the underlying descriptor from a vertex or edge descriptor
template <typename Desc>
Desc convert_to_underlying_descriptor(Desc d) {
  return d;
}
template <typename Desc>
Desc convert_to_underlying_descriptor(reverse_graph_edge_descriptor<Desc> d) {
  return d.underlying_descx_;
}

}  // namespace reverse_graph_detail

template <typename BidirectionalGraph, typename GraphRef = const BidirectionalGraph&>
class reverse_graph {
  using Self = reverse_graph<BidirectionalGraph, GraphRef>;
  using Traits = graph_traits<BidirectionalGraph>;

 public:
  using base_type = BidirectionalGraph;
  using base_ref_type = GraphRef;

  // Constructor
  explicit reverse_graph(GraphRef g) : m_g(g) {}
  // Conversion from reverse_graph on non-const reference to one on const
  // reference
  explicit reverse_graph(const reverse_graph<BidirectionalGraph, BidirectionalGraph&>& o) : m_g(o.m_g) {}

  // Graph requirements
  using vertex_descriptor = typename Traits::vertex_descriptor;
  using edge_descriptor = reverse_graph_detail::reverse_graph_edge_descriptor<typename Traits::edge_descriptor>;
  using directed_category = typename Traits::directed_category;
  using edge_parallel_category = typename Traits::edge_parallel_category;
  using traversal_category = typename Traits::traversal_category;

  using degree_size_type = typename Traits::degree_size_type;
  using vertices_size_type = typename Traits::vertices_size_type;
  using edges_size_type = typename Traits::edges_size_type;

  using vertex_property_type = vertex_property_type<base_type>;
  using edge_property_type = edge_property_type<base_type>;
  using graph_property_type = graph_property_type<base_type>;

  using vertex_bundle_type = vertex_bundle_type<base_type>;
  using edge_bundle_type = edge_bundle_type<base_type>;
  using graph_bundle_type = graph_bundle_type<base_type>;

  using graph_tag = reverse_graph_tag;

  // Bundled properties support
  template <typename Descriptor>
  auto& operator[](Descriptor x) {
    return m_g[reverse_graph_detail::convert_to_underlying_descriptor(x)];
  }

  template <typename Descriptor>
  const auto& operator[](Descriptor x) const {
    return m_g[reverse_graph_detail::convert_to_underlying_descriptor(x)];
  }

  static vertex_descriptor null_vertex() { return Traits::null_vertex(); }

  // would be private, but template friends aren't portable enough.
  // private:
  GraphRef m_g;
};

template <typename BidirectionalGraph>
auto make_reverse_graph(const BidirectionalGraph& g) {
  return reverse_graph<BidirectionalGraph>(g);
}

template <typename BidirectionalGraph>
auto make_reverse_graph(BidirectionalGraph& g) {
  return reverse_graph<BidirectionalGraph, BidirectionalGraph&>(g);
}

template <typename BidirectionalGraph, typename GRef>
auto vertices(const reverse_graph<BidirectionalGraph, GRef>& g) {
  return vertices(g.m_g);
}

template <typename BidirectionalGraph, typename GRef>
auto edges(const reverse_graph<BidirectionalGraph, GRef>& g) {
  return reverse_graph_detail::reverse_edge_range(edges(g.m_g));
}

template <typename BidirectionalGraph, typename GRef>
auto out_edges(const graph_vertex_descriptor_t<BidirectionalGraph> u,
               const reverse_graph<BidirectionalGraph, GRef>& g) {
  return reverse_graph_detail::reverse_edge_range(in_edges(u, g.m_g));
}

template <typename BidirectionalGraph, typename GRef>
auto num_vertices(const reverse_graph<BidirectionalGraph, GRef>& g) {
  return num_vertices(g.m_g);
}

template <typename BidirectionalGraph, typename GRef>
auto num_edges(const reverse_graph<BidirectionalGraph, GRef>& g) {
  return num_edges(g.m_g);
}

template <typename BidirectionalGraph, typename GRef>
auto out_degree(const graph_vertex_descriptor_t<BidirectionalGraph> u,
                const reverse_graph<BidirectionalGraph, GRef>& g) {
  return in_degree(u, g.m_g);
}

template <typename BidirectionalGraph, typename GRef>
graph_vertex_descriptor_t<BidirectionalGraph> vertex(std::size_t v, const reverse_graph<BidirectionalGraph, GRef>& g) {
  return vertex(v, g.m_g);
}

template <typename BidirectionalGraph, typename GRef>
auto edge(graph_vertex_descriptor_t<BidirectionalGraph> u, graph_vertex_descriptor_t<BidirectionalGraph> v,
          const reverse_graph<BidirectionalGraph, GRef>& g) {
  auto [e, found] = edge(v, u, g.m_g);
  return std::pair{reverse_graph_detail::reverse_graph_edge_descriptor(e), found};
}

template <typename BidirectionalGraph, typename GRef>
auto in_edges(graph_vertex_descriptor_t<BidirectionalGraph> u, const reverse_graph<BidirectionalGraph, GRef>& g) {
  return reverse_graph_detail::reverse_edge_range(out_edges(u, g.m_g));
}

template <typename BidirectionalGraph, typename GRef>
auto adjacent_vertices(graph_vertex_descriptor_t<BidirectionalGraph> u,
                       const reverse_graph<BidirectionalGraph, GRef>& g) {
  return adjacency_range(out_edges(u, g), g);
}

template <typename BidirectionalGraph, typename GRef>
auto in_degree(graph_vertex_descriptor_t<BidirectionalGraph> u, const reverse_graph<BidirectionalGraph, GRef>& g) {
  return out_degree(u, g.m_g);
}

template <typename Edge, typename BidirectionalGraph, typename GRef>
auto source(const reverse_graph_detail::reverse_graph_edge_descriptor<Edge>& e,
            const reverse_graph<BidirectionalGraph, GRef>& g) {
  return target(e.underlying_descx_, g.m_g);
}

template <typename Edge, typename BidirectionalGraph, typename GRef>
auto target(const reverse_graph_detail::reverse_graph_edge_descriptor<Edge>& e,
            const reverse_graph<BidirectionalGraph, GRef>& g) {
  return source(e.underlying_descx_, g.m_g);
}

template <typename BidirectionalGraph, typename GRef>
auto degree(graph_vertex_descriptor_t<BidirectionalGraph> u, const reverse_graph<BidirectionalGraph, GRef>& g) {
  return degree(u, g.m_g);
}

namespace reverse_graph_detail {

template <typename PM>
struct reverse_graph_edge_property_map {
 private:
  PM underlying_pm_;

 public:
  using key_type = reverse_graph_edge_descriptor<property_traits_key_t<PM>>;
  using value_type = property_traits_value_t<PM>;
  using reference = property_traits_reference_t<PM>;

  explicit reverse_graph_edge_property_map(const PM& pm) : underlying_pm_(pm) {}

  friend reference get(const reverse_graph_edge_property_map& m, const key_type& e) {
    return get(m.underlying_pm_, e.underlying_descx_);
  }

  friend void put(const reverse_graph_edge_property_map& m, const key_type& e, const value_type& v) {
    put(m.underlying_pm_, e.underlying_descx_, v);
  }

  reference operator[](const key_type& k) const { return (this->underlying_pm_)[k.underlying_descx_]; }
};

}  // namespace reverse_graph_detail

template <typename BidirGraph, typename GRef, typename Property>
struct property_map<reverse_graph<BidirGraph, GRef>, Property> {
  static constexpr bool is_edge_prop_v =
      std::is_same_v<typename properties_detail::property_kind_from_graph<BidirGraph, Property>::type,
                     edge_property_tag>;
  static constexpr bool is_ref_const_v = std::is_const_v<std::remove_reference_t<GRef>>;
  using orig_type = std::conditional_t<is_ref_const_v, property_map_const_t<BidirGraph, Property>,
                                       property_map_t<BidirGraph, Property>>;
  using orig_const_type = property_map_const_t<BidirGraph, Property>;
  using type =
      std::conditional_t<is_edge_prop_v, reverse_graph_detail::reverse_graph_edge_property_map<orig_type>, orig_type>;
  using const_type =
      std::conditional_t<is_edge_prop_v, reverse_graph_detail::reverse_graph_edge_property_map<orig_const_type>,
                         orig_const_type>;
};

template <typename BidirGraph, typename GRef, typename Property>
struct property_map<const reverse_graph<BidirGraph, GRef>, Property> {
  static constexpr bool is_edge_prop_v =
      std::is_same_v<typename properties_detail::property_kind_from_graph<BidirGraph, Property>::type,
                     edge_property_tag>;
  using orig_const_type = property_map_const_t<BidirGraph, Property>;
  using const_type =
      std::conditional_t<is_edge_prop_v, reverse_graph_detail::reverse_graph_edge_property_map<orig_const_type>,
                         orig_const_type>;
  using type = const_type;
};

template <typename BidirGraph, typename GRef, typename Property>
auto get(Property p, reverse_graph<BidirGraph, GRef>& g) {
  return get(p, g.m_g);
}

template <typename BidirGraph, typename GRef, typename Property>
auto get(Property p, const reverse_graph<BidirGraph, GRef>& g) {
  return get(p, std::as_const(g.m_g));  // in case GRef is non-const
}

template <typename BidirectionalGraph, typename GRef, typename Property, typename Key>
auto get(Property p, const reverse_graph<BidirectionalGraph, GRef>& g, const Key& k) {
  return get(get(p, g), k);
}

template <typename BidirGraph, typename GRef, typename Property, typename Key, typename Value>
void put(Property p, reverse_graph<BidirGraph, GRef>& g, const Key& k, const Value& val) {
  put(get(p, g), k, val);
}

// Get the underlying descriptor from a reverse_graph's wrapped edge descriptor

namespace reverse_graph_detail {
template <typename E>
struct underlying_edge_desc_map_type {
  E operator[](const reverse_graph_edge_descriptor<E>& k) const { return k.underlying_descx_; }
};

template <typename E>
E get(underlying_edge_desc_map_type<E> m, const reverse_graph_edge_descriptor<E>& k) {
  return m[k];
}
}  // namespace reverse_graph_detail

template <typename E>
struct property_traits<reverse_graph_detail::underlying_edge_desc_map_type<E>> {
  using key_type = reverse_graph_detail::reverse_graph_edge_descriptor<E>;
  using value_type = E;
  using reference = const E&;
};

template <typename Graph, typename GRef>
struct property_map<reverse_graph<Graph, GRef>, edge_underlying_t> {
 private:
  using ed = graph_edge_descriptor_t<Graph>;

 public:
  using type = reverse_graph_detail::underlying_edge_desc_map_type<ed>;
  using const_type = reverse_graph_detail::underlying_edge_desc_map_type<ed>;
};

template <typename BidirGraph, typename GRef>
auto get(edge_underlying_t /*unused*/, reverse_graph<BidirGraph, GRef>& /*g*/) {
  return reverse_graph_detail::underlying_edge_desc_map_type<graph_edge_descriptor_t<BidirGraph>>();
}

template <typename BidirGraph, typename GRef>
auto get(edge_underlying_t /*unused*/, reverse_graph<BidirGraph, GRef>& /*g*/,
         const graph_edge_descriptor_t<reverse_graph<BidirGraph, GRef>>& k) {
  return k.underlying_descx_;
}

template <typename BidirGraph, typename GRef>
auto get(edge_underlying_t /*unused*/, const reverse_graph<BidirGraph, GRef>& /*g*/) {
  return reverse_graph_detail::underlying_edge_desc_map_type<graph_edge_descriptor_t<BidirGraph>>();
}

template <typename BidirGraph, typename GRef>
auto get(edge_underlying_t /*unused*/, const reverse_graph<BidirGraph, GRef>& /*g*/,
         const graph_edge_descriptor_t<reverse_graph<BidirGraph, GRef>>& k) {
  return k.underlying_descx_;
}

// Access to wrapped graph's graph properties

template <typename BidirectionalGraph, typename GRef, typename Tag, typename Value>
void set_property(const reverse_graph<BidirectionalGraph, GRef>& g, Tag tag, const Value& value) {
  set_property(g.m_g, tag, value);
}

template <typename BidirectionalGraph, typename GRef, typename Tag>
auto get_property(const reverse_graph<BidirectionalGraph, GRef>& g, Tag tag) {
  return get_property(g.m_g, tag);
}

}  // namespace bagl

#endif  // BAGL_BAGL_REVERSE_GRAPH_H_
