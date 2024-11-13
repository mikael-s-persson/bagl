// Copyright (C) 2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_LABELED_GRAPH_H_
#define BAGL_BAGL_LABELED_GRAPH_H_

#include <map>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/container_traits.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

// This file implements a utility for creating mappings from arbitrary
// identifiers to the vertices of a graph.

namespace bagl {

namespace graph_detail {

// Returns true if the selector is the default selector.
template <typename Selector>
constexpr bool is_default_v = std::is_same_v<Selector, default_s>;

// Choose the default map instance. If Label is an unsigned integral type
// the we can use a vector to store the information.
template <typename Label, typename Vertex>
using choose_default_map =
    std::conditional_t<std::is_unsigned_v<Label>, std::vector<Vertex>, std::unordered_map<Label, Vertex>>;

// Generate Label Map
// These type generators are responsible for instantiating an associative
// container for the the labeled graph. Note that the Selector must be
// select a pair associative container or a vecS, which is only valid if
// Label is an integral type.
template <typename Selector, typename Label, typename Vertex>
struct generate_label_map {
  using type = void;
};

template <typename Label, typename Vertex>
struct generate_label_map<vec_s, Label, Vertex> {
  using type = std::vector<Vertex>;
};

template <typename Label, typename Vertex>
struct generate_label_map<set_s, Label, Vertex> {
  using type = std::map<Label, Vertex>;
};

template <typename Label, typename Vertex>
struct generate_label_map<multiset_s, Label, Vertex> {
  using type = std::multimap<Label, Vertex>;
};

template <typename Label, typename Vertex>
struct generate_label_map<unordered_set_s, Label, Vertex> {
  using type = std::unordered_map<Label, Vertex>;
};

template <typename Label, typename Vertex>
struct generate_label_map<unordered_multiset_s, Label, Vertex> {
  using type = std::unordered_multimap<Label, Vertex>;
};

template <typename Selector, typename Label, typename Vertex>
using choose_custom_map = typename generate_label_map<Selector, Label, Vertex>::type;

/**
 * Choose and instantiate an "associative" container. Note that this can
 * also choose vector.
 */
template <typename Selector, typename Label, typename Vertex>
using choose_map = std::conditional_t<is_default_v<Selector>, choose_default_map<Label, Vertex>,
                                      choose_custom_map<Selector, Label, Vertex>>;

// Insert Labeled Vertex
// Tag dispatch on random access containers (i.e., vectors). This function
// basically requires a) that Container is vector<Label> and that Label
// is an unsigned integral value. Note that this will resize the vector
// to accommodate indices.
template <typename Container, typename Graph, typename Label, typename... Prop>
auto insert_labeled_vertex_impl(Container& c, Graph& g, Label const& l,
                                container_traits_detail::random_access_container_tag, Prop&&... p_args) {
  // If the label is out of bounds, resize the vector to accommodate.
  // Resize by 2x the index so we don't cause quadratic insertions over
  // time.
  if (l >= c.size()) {
    c.resize((l + 1) * 2);
  }
  c[l] = add_vertex(g, std::forward<Prop>(p_args)...);
  return std::pair{c[l], true};
}

// Tag dispatch on multi associative containers (i.e. multimaps).
template <typename Container, typename Graph, typename Label, typename... Prop>
auto insert_labeled_vertex_impl(Container& c, Graph& g, Label const& l,
                                container_traits_detail::multiple_associative_container_tag, Prop&&... p_args) {
  // Note that insertion always succeeds so we can add the vertex first
  // and then the mapping to the label.
  auto v = add_vertex(g, std::forward<Prop>(p_args)...);
  c.insert(std::make_pair(l, v));
  return std::pair{v, true};
}

// Tag dispatch on unique associative containers (i.e. maps).
template <typename Container, typename Graph, typename Label, typename... Prop>
auto insert_labeled_vertex_impl(Container& c, Graph& g, Label const& l,
                                container_traits_detail::unique_associative_container_tag, Prop&&... p_args) {
  // Here, we actually have to try the insertion first, and only add
  // the vertex if we get a new element.
  using Vertex = graph_vertex_descriptor_t<Graph>;
  auto [it, inserted] = c.emplace(l, Vertex());
  if (inserted) {
    it->second = add_vertex(g);
    using VProp = vertex_property_type<Graph>;
    put(vertex_all, g, it->second, VProp{std::forward<Prop>(p_args)...});
  }
  return std::pair{it->second, inserted};
}

// Dispatcher
template <typename Container, typename Graph, typename Label, typename... Prop>
auto insert_labeled_vertex(Container& c, Graph& g, Label const& l, Prop&&... p_args) {
  return insert_labeled_vertex_impl(c, g, l, container_traits_detail::container_category(c),
                                    std::forward<Prop>(p_args)...);
}

// Find Labeled Vertex
// Tag dispatch for sequential maps (i.e., vectors).
template <typename Container, typename Graph, typename Label>
auto find_labeled_vertex(Container const& c, Graph const&, Label const& l,
                         container_traits_detail::random_access_container_tag) {
  return l < c.size() ? c[l] : graph_traits<Graph>::null_vertex();
}

// Tag dispatch for pair associative maps (more or less).
template <typename Container, typename Graph, typename Label>
auto find_labeled_vertex(Container const& c, Graph const&, Label const& l,
                         container_traits_detail::associative_container_tag) {
  auto i = c.find(l);
  return i != c.end() ? i->second : graph_traits<Graph>::null_vertex();
}

// Dispatcher
template <typename Container, typename Graph, typename Label>
auto find_labeled_vertex(Container const& c, Graph const& g, Label const& l) {
  return find_labeled_vertex(c, g, l, container_traits_detail::container_category(c));
}

// Put Vertex Label
// Tag dispatch on vectors.
template <typename Container, typename Label, typename Graph, typename Vertex>
bool put_vertex_label(Container& c, Graph const&, Label const& l, Vertex v,
                      container_traits_detail::random_access_container_tag) {
  // If the element is already occupied, then we probably don't want to
  // overwrite it.
  if (c[l] == graph_traits<Graph>::null_vertex()) {
    return false;
  }
  c[l] = v;
  return true;
}

// Attempt the insertion and return its result.
template <typename Container, typename Label, typename Graph, typename Vertex>
bool put_vertex_label(Container& c, Graph const&, Label const& l, Vertex v,
                      container_traits_detail::unique_associative_container_tag) {
  return c.insert(std::make_pair(l, v)).second;
}

// Insert the pair and return true.
template <typename Container, typename Label, typename Graph, typename Vertex>
bool put_vertex_label(Container& c, Graph const&, Label const& l, Vertex v,
                      container_traits_detail::multiple_associative_container_tag) {
  c.insert(std::make_pair(l, v));
  return true;
}

// Dispatcher
template <typename Container, typename Label, typename Graph, typename Vertex>
bool put_vertex_label(Container& c, Graph const& g, Label const& l, Vertex v) {
  return put_vertex_label(c, g, l, v, container_traits_detail::container_category(c));
}

// Remove Labeled Vertex
// Tag dispatch on random access containers (i.e., vectors)
template <typename Container, typename Label, typename Graph>
void remove_labeled_vertex(Container& c, Graph& g, Label const& l,
                           container_traits_detail::random_access_container_tag) {
  if (l < c.size()) {
    remove_vertex(c[l], g);
    c.erase(c.begin() + l);
  }
}

// Tag dispatch on multi associative containers (i.e. multimaps).
template <typename Container, typename Label, typename Graph>
void remove_labeled_vertex(Container& c, Graph& g, Label const& l,
                           container_traits_detail::multiple_associative_container_tag) {
  auto c_it = c.find(l);
  if (c_it != c.end()) {
    remove_vertex(c_it->second, g);
    c.erase(c_it);
  }
}

// Tag dispatch on unique associative containers (i.e. maps).
template <typename Container, typename Label, typename Graph>
void remove_labeled_vertex(Container& c, Graph& g, Label const& l,
                           container_traits_detail::unique_associative_container_tag) {
  auto c_it = c.find(l);
  if (c_it != c.end()) {
    remove_vertex(c_it->second, g);
    c.erase(c_it);
  }
}

// Dispatcher
template <typename Container, typename Label, typename Graph>
void remove_labeled_vertex(Container& c, Graph& g, Label const& l) {
  remove_labeled_vertex(c, g, l, container_traits_detail::container_category(c));
}

}  // namespace graph_detail

struct labeled_graph_class_tag {};

// This class is responsible for the deduction and declaration of type names
// for the labeled_graph class template.
template <typename Graph, typename Label, typename Selector>
struct labeled_graph_types {
  using graph_type = Graph;

  // Label and maps
  using label_type = Label;
  using map_type = graph_detail::choose_map<Selector, Label, graph_vertex_descriptor_t<Graph>>;
};

// The labeled_graph class is a graph adaptor that maintains a mapping between
// vertex labels and vertex descriptors.
//
// TODO: This class is somewhat redundant for adjacency_list<*, vec_s>  if
// the intended label is an unsigned int (and perhaps some other cases), but
// it does avoid some weird ambiguities (i.e. adding a vertex with a label that
// does not match its target index).
template <typename Graph, typename Label, typename Selector = default_s>
class labeled_graph : protected labeled_graph_types<Graph, Label, Selector> {
  using Base = labeled_graph_types<Graph, Label, Selector>;

 public:
  using graph_tag = labeled_graph_class_tag;

  using graph_type = typename Base::graph_type;
  using traits_type = graph_traits<graph_type>;
  using vertex_descriptor = typename traits_type::vertex_descriptor;
  using edge_descriptor = typename traits_type::edge_descriptor;
  using directed_category = typename traits_type::directed_category;
  using edge_parallel_category = typename traits_type::edge_parallel_category;
  using traversal_category = typename traits_type::traversal_category;

  using degree_size_type = typename traits_type::degree_size_type;

  using vertices_size_type = typename traits_type::vertices_size_type;
  using edges_size_type = typename traits_type::edges_size_type;

  using graph_property_type = typename graph_type::graph_property_type;
  using graph_bundled = typename graph_type::graph_bundled;

  using vertex_property_type = typename graph_type::vertex_property_type;
  using vertex_bundled = typename graph_type::vertex_bundled;

  using edge_property_type = typename graph_type::edge_property_type;
  using edge_bundled = typename graph_type::edge_bundled;

  using label_type = typename Base::label_type;
  using map_type = typename Base::map_type;

 public:
  labeled_graph(graph_property_type const& gp = graph_property_type()) : graph_(0, gp), map_() {}

  labeled_graph(labeled_graph const& x) : graph_(x.graph_), map_(x.map_) {}

  // This constructor can only be used if map_type supports positional
  // range insertion (i.e. its a vector). This is the only case where we can
  // try to guess the intended labels for graph.
  labeled_graph(vertices_size_type n, graph_property_type const& gp = graph_property_type()) : graph_(n, gp), map_() {
    auto rng = vertices(graph_);
    map_.insert(map_.end(), rng.begin(), rng.end());
  }

  // Construct a graph over n vertices, each of which receives a label from
  // the range [l, l + n). Note that the graph is not directly constructed
  // over the n vertices, but added sequentially. This constructor is
  // necessarily slower than the underlying counterpart.
  template <typename LabelIter>
  labeled_graph(vertices_size_type n, LabelIter l, graph_property_type const& gp = graph_property_type()) : graph_(gp) {
    while (n-- > 0) {
      add_vertex(*l++);
    }
  }

  // Construct the graph over n vertices each of which has a label in the
  // range [l, l + n) and a property in the range [p, p + n).
  template <typename LabelIter, typename PropIter>
  labeled_graph(vertices_size_type n, LabelIter l, PropIter p, graph_property_type const& gp = graph_property_type())
      : graph_(gp) {
    while (n-- > 0) {
      add_vertex(*l++, *p++);
    }
  }

  labeled_graph& operator=(labeled_graph const& x) {
    graph_ = x.graph_;
    map_ = x.map_;
    return *this;
  }

  // Graph Accessors
  graph_type& graph() { return graph_; }
  graph_type const& graph() const { return graph_; }

  // Create a new label for the given vertex, returning false, if the label
  // cannot be created.
  bool label_vertex(vertex_descriptor v, Label const& l) { return graph_detail::put_vertex_label(map_, graph_, l, v); }

  // Add a vertex to the graph, returning the descriptor. If the vertices
  // are uniquely labeled and the label already exists within the graph,
  // then no vertex is added, and the returned descriptor refers to the
  // existing vertex. A vertex property can be given as a parameter, if
  // needed.
  vertex_descriptor add_vertex(Label const& l) {
    return graph_detail::insert_labeled_vertex(map_, graph_, l, vertex_property_type()).first;
  }

  vertex_descriptor add_vertex(Label const& l, vertex_property_type const& p) {
    return graph_detail::insert_labeled_vertex(map_, graph_, l, p).first;
  }

  // Insert a vertex into the graph, returning a pair containing the
  // descriptor of a vertex and a boolean value that describes whether or not
  // a new vertex was inserted. If vertices are not uniquely labeled, then
  // insertion will always succeed.
  std::pair<vertex_descriptor, bool> insert_vertex(Label const& l) {
    return graph_detail::insert_labeled_vertex(map_, graph_, l, vertex_property_type());
  }

  std::pair<vertex_descriptor, bool> insert_vertex(Label const& l, vertex_property_type const& p) {
    return graph_detail::insert_labeled_vertex(map_, graph_, l, p);
  }

  // Remove the vertex with the given label.
  void remove_vertex(Label const& l) { return graph_detail::remove_labeled_vertex(map_, graph_, l); }

  // Return a descriptor for the given label.
  vertex_descriptor vertex(Label const& l) const { return graph_detail::find_labeled_vertex(map_, graph_, l); }

  // Bundled Properties
  // Lookup the requested vertex and return the bundle.
  vertex_bundled& operator[](Label const& l) { return graph_[vertex(l)]; }

  vertex_bundled const& operator[](Label const& l) const { return graph_[vertex(l)]; }

  // Delegate edge lookup to the underlying graph.
  edge_bundled& operator[](edge_descriptor e) { return graph_[e]; }

  edge_bundled const& operator[](edge_descriptor e) const { return graph_[e]; }

  /** Return a null descriptor */
  static vertex_descriptor null_vertex() { return graph_traits<graph_type>::null_vertex(); }

 private:
  graph_type graph_;
  map_type map_;
};

// The partial specialization over graph pointers allows the construction
// of temporary labeled graph objects. In this case, the labels are destructed
// when the wrapper goes out of scope.
template <typename Graph, typename Label, typename Selector>
class labeled_graph<Graph*, Label, Selector> : protected labeled_graph_types<Graph, Label, Selector> {
  using Base = labeled_graph_types<Graph, Label, Selector>;

 public:
  using graph_tag = labeled_graph_class_tag;

  using graph_type = typename Base::graph_type;
  using traits_type = graph_traits<graph_type>;
  using vertex_descriptor = typename traits_type::vertex_descriptor;
  using edge_descriptor = typename traits_type::edge_descriptor;
  using directed_category = typename traits_type::directed_category;
  using edge_parallel_category = typename traits_type::edge_parallel_category;
  using traversal_category = typename traits_type::traversal_category;

  using degree_size_type = typename traits_type::degree_size_type;
  using vertices_size_type = typename traits_type::vertices_size_type;
  using edges_size_type = typename traits_type::edges_size_type;

  using vertex_property_type = typename graph_type::vertex_property_type;
  using edge_property_type = typename graph_type::edge_property_type;
  using graph_property_type = typename graph_type::graph_property_type;
  using vertex_bundled = typename graph_type::vertex_bundled;
  using edge_bundled = typename graph_type::edge_bundled;

  using label_type = typename Base::label_type;
  using map_type = typename Base::map_type;

  labeled_graph(graph_type* g) : graph_(g) {}

  // Graph Access
  graph_type& graph() { return *graph_; }
  graph_type const& graph() const { return *graph_; }

  // Create a new label for the given vertex, returning false, if the label
  // cannot be created.
  bool label_vertex(vertex_descriptor v, Label const& l) { return graph_detail::put_vertex_label(map_, *graph_, l, v); }

  // @name Add Vertex
  vertex_descriptor add_vertex(Label const& l) {
    return graph_detail::insert_labeled_vertex(map_, *graph_, l, vertex_property_type()).first;
  }

  vertex_descriptor add_vertex(Label const& l, vertex_property_type const& p) {
    return graph_detail::insert_labeled_vertex(map_, *graph_, l, p).first;
  }

  std::pair<vertex_descriptor, bool> insert_vertex(Label const& l) {
    return graph_detail::insert_labeled_vertex(map_, *graph_, l, vertex_property_type());
  }

  // Try to insert a vertex with the given label.
  std::pair<vertex_descriptor, bool> insert_vertex(Label const& l, vertex_property_type const& p) {
    return graph_detail::insert_labeled_vertex(map_, *graph_, l, p);
  }

  // Remove the vertex with the given label.
  void remove_vertex(Label const& l) { return remove_vertex(vertex(l), *graph_); }

  // Return a descriptor for the given label.
  vertex_descriptor vertex(Label const& l) const { return graph_detail::find_labeled_vertex(map_, *graph_, l); }

  // Bundled Properties
  // Lookup the requested vertex and return the bundle.
  vertex_bundled& operator[](Label const& l) { return (*graph_)[vertex(l)]; }

  vertex_bundled const& operator[](Label const& l) const { return (*graph_)[vertex(l)]; }

  // Delegate edge lookup to the underlying graph.
  edge_bundled& operator[](edge_descriptor e) { return (*graph_)[e]; }

  edge_bundled const& operator[](edge_descriptor e) const { return (*graph_)[e]; }

  static vertex_descriptor null_vertex() { return graph_traits<graph_type>::null_vertex(); }

 private:
  graph_type* graph_;
  map_type map_;
};

#define LABELED_GRAPH_PARAMS typename G, typename L, typename S
#define LABELED_GRAPH labeled_graph<G, L, S>

// Labeled Graph
template <LABELED_GRAPH_PARAMS>
bool label_vertex(typename LABELED_GRAPH::vertex_descriptor v, typename LABELED_GRAPH::label_type const l,
                  LABELED_GRAPH& g) {
  return g.label_vertex(v, l);
}

template <LABELED_GRAPH_PARAMS>
auto vertex_by_label(typename LABELED_GRAPH::label_type const l, LABELED_GRAPH& g) {
  return g.vertex(l);
}

// Graph
template <LABELED_GRAPH_PARAMS>
auto edge(typename LABELED_GRAPH::vertex_descriptor const& u, typename LABELED_GRAPH::vertex_descriptor const& v,
          LABELED_GRAPH const& g) {
  return edge(u, v, g.graph());
}

// Labeled Extensions
template <LABELED_GRAPH_PARAMS>
auto edge_by_label(typename LABELED_GRAPH::label_type const& u, typename LABELED_GRAPH::label_type const& v,
                   LABELED_GRAPH const& g) {
  return edge(g.vertex(u), g.vertex(v), g);
}

// Incidence Graph
template <LABELED_GRAPH_PARAMS>
auto out_edges(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH const& g) {
  return out_edges(v, g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto out_degree(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH const& g) {
  return out_degree(v, g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto source(typename LABELED_GRAPH::edge_descriptor e, LABELED_GRAPH const& g) {
  return source(e, g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto target(typename LABELED_GRAPH::edge_descriptor e, LABELED_GRAPH const& g) {
  return target(e, g.graph());
}

// Bidirectional Graph
template <LABELED_GRAPH_PARAMS>
auto in_edges(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH const& g) {
  return in_edges(v, g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto in_degree(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH const& g) {
  return in_degree(v, g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto degree(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH const& g) {
  return degree(v, g.graph());
}

// Adjacency Graph
template <LABELED_GRAPH_PARAMS>
auto adjacent_vertices(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH const& g) {
  return adjacent_vertices(v, g.graph());
}

// VertexListGraph
template <LABELED_GRAPH_PARAMS>
auto vertices(LABELED_GRAPH const& g) {
  return vertices(g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto num_vertices(LABELED_GRAPH const& g) {
  return num_vertices(g.graph());
}

// EdgeListGraph
template <LABELED_GRAPH_PARAMS>
auto edges(LABELED_GRAPH const& g) {
  return edges(g.graph());
}

template <LABELED_GRAPH_PARAMS>
auto num_edges(LABELED_GRAPH const& g) {
  return num_edges(g.graph());
}

// Property Lookup
namespace graph_detail {
struct labeled_graph_vertex_property_selector {
  template <class LabeledGraph, class Property, class Tag>
  struct bind_ {
    using Graph = typename LabeledGraph::graph_type;
    using PropertyMap = property_map<Graph, Tag>;
    using type = typename PropertyMap::type;
    using const_type = typename PropertyMap::const_type;
  };
};

struct labeled_graph_edge_property_selector {
  template <class LabeledGraph, class Property, class Tag>
  struct bind_ {
    using Graph = typename LabeledGraph::graph_type;
    using PropertyMap = property_map<Graph, Tag>;
    using type = typename PropertyMap::type;
    using const_type = typename PropertyMap::const_type;
  };
};
}  // namespace graph_detail

template <>
struct vertex_property_selector<labeled_graph_class_tag> {
  using type = graph_detail::labeled_graph_vertex_property_selector;
};

template <>
struct edge_property_selector<labeled_graph_class_tag> {
  using type = graph_detail::labeled_graph_edge_property_selector;
};

// Property Graph
template <LABELED_GRAPH_PARAMS, typename Prop>
auto get(Prop p, LABELED_GRAPH& g) {
  return get(p, g.graph());
}

template <LABELED_GRAPH_PARAMS, typename Prop>
auto get(Prop p, LABELED_GRAPH const& g) {
  return get(p, g.graph());
}

template <LABELED_GRAPH_PARAMS, typename Prop, typename Key>
auto get(Prop p, LABELED_GRAPH const& g, const Key& k) {
  return get(p, g.graph(), k);
}

template <LABELED_GRAPH_PARAMS, typename Prop, typename Key, typename Value>
void put(Prop p, LABELED_GRAPH& g, Key const& k, Value const& v) {
  put(p, g.graph(), k, v);
}

// Mutable Graph
template <LABELED_GRAPH_PARAMS>
auto add_edge(typename LABELED_GRAPH::vertex_descriptor const& u, typename LABELED_GRAPH::vertex_descriptor const& v,
              LABELED_GRAPH& g) {
  return add_edge(u, v, g.graph());
}

template <LABELED_GRAPH_PARAMS, typename... EPArgs>
auto add_edge(typename LABELED_GRAPH::vertex_descriptor const& u, typename LABELED_GRAPH::vertex_descriptor const& v,
              LABELED_GRAPH& g, EPArgs&&... ep_args) {
  return add_edge(u, v, g.graph(), std::forward<EPArgs>(ep_args)...);
}

template <LABELED_GRAPH_PARAMS>
void clear_vertex(typename LABELED_GRAPH::vertex_descriptor v, LABELED_GRAPH& g) {
  return clear_vertex(v, g.graph());
}

template <LABELED_GRAPH_PARAMS>
void remove_edge(typename LABELED_GRAPH::edge_descriptor e, LABELED_GRAPH& g) {
  return remove_edge(e, g.graph());
}

template <LABELED_GRAPH_PARAMS>
void remove_edge(typename LABELED_GRAPH::vertex_descriptor u, typename LABELED_GRAPH::vertex_descriptor v,
                 LABELED_GRAPH& g) {
  return remove_edge(u, v, g.graph());
}

// Labeled extensions
template <LABELED_GRAPH_PARAMS, typename... EPArgs>
auto add_edge_by_label(typename LABELED_GRAPH::label_type const& u, typename LABELED_GRAPH::label_type const& v,
                       LABELED_GRAPH& g, EPArgs&&... ep_args) {
  return add_edge(g.vertex(u), g.vertex(v), g, std::forward<EPArgs>(ep_args)...);
}

template <LABELED_GRAPH_PARAMS>
void clear_vertex_by_label(typename LABELED_GRAPH::label_type const& l, LABELED_GRAPH& g) {
  clear_vertex(g.vertex(l), g.graph());
}

template <LABELED_GRAPH_PARAMS>
void remove_edge_by_label(typename LABELED_GRAPH::label_type const& u, typename LABELED_GRAPH::label_type const& v,
                          LABELED_GRAPH& g) {
  remove_edge(g.vertex(u), g.vertex(v), g.graph());
}

// Labeled Mutable Graph
// The labeled mutable graph hides the add_ and remove_ vertex functions from
// the mutable graph concept. Note that the remove_vertex is hidden because
// removing the vertex without its key could leave a dangling reference in
// the map.
template <LABELED_GRAPH_PARAMS, typename... VPArgs>
auto add_vertex(LABELED_GRAPH& g, typename LABELED_GRAPH::label_type const& l, VPArgs&&... vp_args) {
  return g.add_vertex(l, std::forward<VPArgs>(vp_args)...);
}

template <LABELED_GRAPH_PARAMS>
void remove_vertex(typename LABELED_GRAPH::label_type const& l, LABELED_GRAPH& g) {
  g.remove_vertex(l);
}

#undef LABELED_GRAPH_PARAMS
#undef LABELED_GRAPH

}  // namespace bagl

// Pull the labeled graph traits
#include "bagl/detail/labeled_graph_traits.h"

#endif  // BAGL_BAGL_LABELED_GRAPH_H_
