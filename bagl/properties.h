// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PROPERTIES_H_
#define BAGL_BAGL_PROPERTIES_H_

#include <type_traits>

#include "bagl/constant_property_map.h"
#include "bagl/graph_traits.h"
#include "bagl/null_property_map.h"
#include "bagl/property.h"
#include "bagl/property_map.h"

namespace bagl {

namespace concepts {

template <typename PMap, typename G>
concept ReadableVertexPropertyMap = ReadablePropertyMap<PMap, graph_vertex_descriptor_t<G>>;
template <typename PMap, typename G>
concept ReadableEdgePropertyMap = ReadablePropertyMap<PMap, graph_edge_descriptor_t<G>>;

template <typename PMap, typename G>
concept WritableVertexPropertyMap = WritablePropertyMap<PMap, graph_vertex_descriptor_t<G>>;
template <typename PMap, typename G>
concept WritableEdgePropertyMap = WritablePropertyMap<PMap, graph_edge_descriptor_t<G>>;

template <typename PMap, typename G>
concept ReadWriteVertexPropertyMap = ReadWritePropertyMap<PMap, graph_vertex_descriptor_t<G>>;
template <typename PMap, typename G>
concept ReadWriteEdgePropertyMap = ReadWritePropertyMap<PMap, graph_edge_descriptor_t<G>>;

template <typename PMap, typename G>
concept LvalueVertexPropertyMap = LvaluePropertyMap<PMap, graph_vertex_descriptor_t<G>>;
template <typename PMap, typename G>
concept LvalueEdgePropertyMap = LvaluePropertyMap<PMap, graph_edge_descriptor_t<G>>;

template <typename PMap, typename G>
concept MutableLvalueVertexPropertyMap = MutableLvaluePropertyMap<PMap, graph_vertex_descriptor_t<G>>;
template <typename PMap, typename G>
concept MutableLvalueEdgePropertyMap = MutableLvaluePropertyMap<PMap, graph_edge_descriptor_t<G>>;

} // namespace concepts

enum class default_color_type { white_color, gray_color, green_color, red_color, black_color };

template <class ColorValue>
struct color_traits {
  static constexpr default_color_type white() { return default_color_type::white_color; }
  static constexpr default_color_type gray() { return default_color_type::gray_color; }
  static constexpr default_color_type green() { return default_color_type::green_color; }
  static constexpr default_color_type red() { return default_color_type::red_color; }
  static constexpr default_color_type black() { return default_color_type::black_color; }
};

// These functions are now obsolete, replaced by color_traits.
constexpr default_color_type white(default_color_type) { return default_color_type::white_color; }
constexpr default_color_type gray(default_color_type) { return default_color_type::gray_color; }
constexpr default_color_type green(default_color_type) { return default_color_type::green_color; }
constexpr default_color_type red(default_color_type) { return default_color_type::red_color; }
constexpr default_color_type black(default_color_type) { return default_color_type::black_color; }

// Define a bunch of pre-defined properties that come up a lot.

#define BAGL_DEF_PROPERTY(KIND, NAME) \
  struct KIND##_##NAME##_t {          \
    using kind = KIND##_property_tag; \
  };                                  \
  constexpr KIND##_##NAME##_t KIND##_##NAME = {};

BAGL_DEF_PROPERTY(vertex, index);
BAGL_DEF_PROPERTY(vertex, index1);
BAGL_DEF_PROPERTY(vertex, index2);
BAGL_DEF_PROPERTY(vertex, root);
BAGL_DEF_PROPERTY(edge, index);
BAGL_DEF_PROPERTY(edge, name);
BAGL_DEF_PROPERTY(edge, weight);
BAGL_DEF_PROPERTY(edge, weight2);
BAGL_DEF_PROPERTY(edge, color);
BAGL_DEF_PROPERTY(vertex, name);
BAGL_DEF_PROPERTY(graph, name);
BAGL_DEF_PROPERTY(vertex, distance);
BAGL_DEF_PROPERTY(vertex, distance2);
BAGL_DEF_PROPERTY(vertex, heuristic);
BAGL_DEF_PROPERTY(vertex, rhs);
BAGL_DEF_PROPERTY(vertex, key);
BAGL_DEF_PROPERTY(vertex, position);
BAGL_DEF_PROPERTY(vertex, color);
BAGL_DEF_PROPERTY(vertex, degree);
BAGL_DEF_PROPERTY(vertex, in_degree);
BAGL_DEF_PROPERTY(vertex, out_degree);
BAGL_DEF_PROPERTY(vertex, current_degree);
BAGL_DEF_PROPERTY(vertex, priority);
BAGL_DEF_PROPERTY(vertex, discover_time);
BAGL_DEF_PROPERTY(vertex, finish_time);
BAGL_DEF_PROPERTY(vertex, predecessor);
BAGL_DEF_PROPERTY(vertex, rank);
BAGL_DEF_PROPERTY(vertex, centrality);
BAGL_DEF_PROPERTY(vertex, lowpoint);
BAGL_DEF_PROPERTY(vertex, potential);
BAGL_DEF_PROPERTY(vertex, update);
BAGL_DEF_PROPERTY(vertex, underlying);
BAGL_DEF_PROPERTY(vertex, second_bundle);
BAGL_DEF_PROPERTY(edge, reverse);
BAGL_DEF_PROPERTY(edge, capacity);
BAGL_DEF_PROPERTY(edge, flow);
BAGL_DEF_PROPERTY(edge, residual_capacity);
BAGL_DEF_PROPERTY(edge, centrality);
BAGL_DEF_PROPERTY(edge, discover_time);
BAGL_DEF_PROPERTY(edge, update);
BAGL_DEF_PROPERTY(edge, finished);
BAGL_DEF_PROPERTY(edge, underlying);
BAGL_DEF_PROPERTY(graph, visitor);

// These tags are used to denote the owners and local descriptors
// for the vertices and edges of a distributed graph.
BAGL_DEF_PROPERTY(vertex, global);
BAGL_DEF_PROPERTY(vertex, owner);
BAGL_DEF_PROPERTY(vertex, local);
BAGL_DEF_PROPERTY(edge, global);
BAGL_DEF_PROPERTY(edge, owner);
BAGL_DEF_PROPERTY(edge, local);
BAGL_DEF_PROPERTY(vertex, local_index);
BAGL_DEF_PROPERTY(edge, local_index);

#undef BAGL_DEF_PROPERTY

namespace properties_detail {

template <typename G, typename Tag>
struct property_kind_from_graph : property_kind<Tag> {};

template <typename G, typename R, typename T>
struct property_kind_from_graph<G, R T::*> {
  using type =
      std::conditional_t<std::is_convertible_v<vertex_bundle_type<G>&, T&>, vertex_property_tag,
                         std::conditional_t<std::is_convertible_v<edge_bundle_type<G>&, T&>, edge_property_tag,
                                            std::conditional_t<std::is_convertible_v<graph_bundle_type<G>&, T&>,
                                                               graph_property_tag, void> > >;
};

struct dummy_edge_property_selector {
  template <typename Graph, typename Property, typename Tag>
  struct bind_ {
    using type = identity_property_map;
    using const_type = identity_property_map;
  };
};
struct dummy_vertex_property_selector {
  template <typename Graph, typename Property, typename Tag>
  struct bind_ {
    using type = identity_property_map;
    using const_type = identity_property_map;
  };
};

}  // namespace properties_detail

// Graph classes can either partially specialize property_map
// or they can specialize these two selector classes.
template <typename GraphTag>
struct edge_property_selector {
  using type = properties_detail::dummy_edge_property_selector;
};

template <typename GraphTag>
struct vertex_property_selector {
  using type = properties_detail::dummy_vertex_property_selector;
};

namespace properties_detail {

template <typename Graph, typename Enable = void>
struct graph_tag_or_void {
  using type = void;
};

template <typename Graph>
struct graph_tag_or_void<Graph, std::void_t<typename Graph::graph_tag> > {
  using type = typename Graph::graph_tag;
};

template <typename Graph, typename PropertyTag>
struct edge_property_map : edge_property_selector<typename graph_tag_or_void<Graph>::type>::type::template bind_<
                               Graph, typename edge_property_type<Graph>::type, PropertyTag> {};
template <typename Graph, typename PropertyTag>
struct vertex_property_map : vertex_property_selector<typename graph_tag_or_void<Graph>::type>::type::template bind_<
                                 Graph, typename vertex_property_type<Graph>::type, PropertyTag> {};
}  // namespace properties_detail

template <typename Graph, typename Property, typename Enable = void>
struct property_map
    : std::conditional_t<std::is_same_v<typename properties_detail::property_kind_from_graph<Graph, Property>::type,
                                        edge_property_tag>,
                         properties_detail::edge_property_map<Graph, Property>,
                         properties_detail::vertex_property_map<Graph, Property> > {};
template <typename Graph, typename Property>
using property_map_const_t = typename property_map<Graph, Property>::const_type;
template <typename Graph, typename Property>
using property_map_t = typename property_map<Graph, Property>::type;

// shortcut for accessing the value type of the property map
template <typename Graph, typename Property>
class property_map_value {
  using PMap = property_map_const_t<Graph, Property>;

 public:
  using type = property_traits_value_t<PMap>;
};

template <typename Graph, typename Property>
using graph_property = typename property_value<graph_property_type<Graph>, Property>::type;

template <typename Graph>
using vertex_property = vertex_property_type<Graph>;
template <typename Graph>
using edge_property = edge_property_type<Graph>;

template <typename Graph>
class degree_property_map : public put_get_helper<degree_property_map<Graph> > {
 public:
  using key_type = typename graph_traits<Graph>::vertex_descriptor;
  using value_type = typename graph_traits<Graph>::degree_size_type;
  using reference = value_type;
  using category = readable_property_map_tag;
  explicit degree_property_map(const Graph& g) : g_(g) {}
  value_type operator[](const key_type& v) const { return degree(v, g_); }

 private:
  const Graph& g_;
};
template <typename Graph>
auto make_degree_map(const Graph& g) {
  return degree_property_map<Graph>(g);
}

//========================================================================
// Iterator Property Map Generating Functions contributed by
// Kevin Vanhorn. (see also the property map generating functions
// in bagl/property_map.h)

// A helper function for creating a vertex property map out of a
// random access iterator and the internal vertex index map from a
// graph.
template <typename PropertyGraph, std::random_access_iterator RAIter>
auto make_iterator_vertex_map(RAIter iter, const PropertyGraph& g) {
  return make_iterator_property_map(iter, get(vertex_index, g));
}

// Use this next function when vertex_descriptor is known to be an
// integer type, with values ranging from 0 to num_vertices(g).
//
template <std::random_access_iterator RAIter>
auto make_iterator_vertex_map(RAIter iter) {
  return make_iterator_property_map(iter, identity_property_map());
}

template <typename PropertyGraph, typename RAContainer>
auto make_container_vertex_map(RAContainer& c, const PropertyGraph& g) {
  assert(c.size() >= num_vertices(g));
  return make_iterator_vertex_map(c.begin(), g);
}

template <typename RAContainer>
auto make_container_vertex_map(RAContainer& c) {
  return make_iterator_vertex_map(c.begin());
}

// NOTE: These functions are declared, but never defined since they need to
// be overloaded by graph implementations. However, we need them to be
// declared for the functions below.
template <typename Graph, typename Tag>
graph_property<Graph, graph_bundle_t>& get_property(Graph& g, Tag);

template <typename Graph, typename Tag>
const graph_property<Graph, graph_bundle_t>& get_property(Graph const& g, Tag);

// NOTE: This operation is a simple adaptor over the overloaded get_property
// operations.
template <typename Graph>
graph_property<Graph, graph_bundle_t>& get_property(Graph& g) {
  return get_property(g, graph_bundle);
}

template <typename Graph>
const graph_property<Graph, graph_bundle_t>& get_property(const Graph& g) {
  return get_property(g, graph_bundle);
}

}  // namespace bagl

#endif  // BAGL_BAGL_PROPERTIES_H_
