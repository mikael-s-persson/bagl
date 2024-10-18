// Copyright 2001 University of Notre Dame.
// Authors: Jeremy G. Siek and Lie-Quan Lee
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SUBGRAPH_H_
#define BAGL_BAGL_SUBGRAPH_H_

// UNDER CONSTRUCTION

#include <cassert>
#include <list>
#include <map>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "bagl/graph_mutability_traits.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

struct subgraph_tag {};

// The local_property and global_property functions are used to create
// structures that determine the lookup strategy for properties in subgraphs.
// Note that the nested kind member is used to help interoperate with actual
// Property types.
template <typename T>
struct local_property {
  using kind = T;
  explicit local_property(T x) : value(x) {}
  T value;
};

template <typename T>
auto local(T x) {
  return local_property<T>(x);
}

template <typename T>
struct global_property {
  using kind = T;
  explicit global_property(T x) : value(x) {}
  T value;
};

template <typename T>
auto global(T x) {
  return global_property<T>(x);
}

// Invariants of an induced subgraph:
//   - If vertex u is in subgraph g, then u must be in g.parent().
//   - If edge e is in subgraph g, then e must be in g.parent().
//   - If edge e=(u,v) is in the root graph, then edge e
//     is also in any subgraph that contains both vertex u and v.

// The Graph template parameter must have a vertex_index and edge_index
// internal property. It is assumed that the vertex indices are assigned
// automatically by the graph during a call to add_vertex(). It is not
// assumed that the edge vertices are assigned automatically, they are
// explicitly assigned here.

template <typename Graph>
class subgraph {
  using Traits = graph_traits<Graph>;
  using ChildrenList = std::vector<std::unique_ptr<subgraph<Graph>>>;

 public:
  // Graph requirements
  using vertex_descriptor = typename Traits::vertex_descriptor;
  using edge_descriptor = typename Traits::edge_descriptor;
  using directed_category = typename Traits::directed_category;
  using edge_parallel_category = typename Traits::edge_parallel_category;
  using traversal_category = typename Traits::traversal_category;

  using degree_size_type = typename Traits::degree_size_type;
  using vertices_size_type = typename Traits::vertices_size_type;
  using edges_size_type = typename Traits::edges_size_type;

  using graph_tag = subgraph_tag;
  using graph_type = Graph;

  using edge_property_type = edge_property_type<Graph>;
  using vertex_property_type = vertex_property_type<Graph>;
  using graph_property_type = graph_property_type<Graph>;

  using vertex_bundle_type = vertex_bundle_type<Graph>;
  using edge_bundle_type = edge_bundle_type<Graph>;
  using graph_bundle_type = graph_bundle_type<Graph>;

  // Create the main graph, the root of the subgraph tree
  subgraph() : m_parent(0), m_edge_counter(0) {}

  explicit subgraph(const graph_property_type& p) : m_graph(p), m_parent(0), m_edge_counter(0) {}

  explicit subgraph(std::size_t n, const graph_property_type& p = graph_property_type())
      : m_graph(n, p), m_parent(0), m_edge_counter(0), m_global_vertex(n) {
    std::size_t i = 0;
    for (auto v : vertices(m_graph)) {
      m_global_vertex[i++] = v;
    }
  }

  // copy constructor
  subgraph(const subgraph& x) : m_parent(x.m_parent), m_edge_counter(0) {
    if (x.is_root()) {
      m_graph = x.m_graph;
      m_edge_counter = x.m_edge_counter;
      m_global_vertex = x.m_global_vertex;
      m_global_edge = x.m_global_edge;
    } else {
      get_property(*this) = get_property(x);
      for (auto v : vertices(x)) {
        add_vertex(x.local_to_global(v), *this);
      }
    }
    // Do a deep copy (recursive).
    // Only the root graph is copied, the subgraphs contain
    // only references to the global vertices they own.
    for (auto i : x.children()) {
      m_children.emplace_back(std::make_unique<subgraph<Graph>>(i));
      m_children.back()->m_parent = this;
    }
  }

  // Return a null vertex descriptor for the graph.
  static vertex_descriptor null_vertex() { return Traits::null_vertex(); }

  // Create a subgraph
  subgraph<Graph>& create_subgraph() {
    m_children.emplace_back(std::make_unique<subgraph<Graph>>());
    m_children.back()->m_parent = this;
    return *m_children.back();
  }

  // Create a subgraph with the specified vertex set.
  template <std::ranges::input_range VertexRange>
  subgraph<Graph>& create_subgraph(const VertexRange& v_range) {
    m_children.emplace_back(std::make_unique<subgraph<Graph>>());
    m_children.back()->m_parent = this;
    for (auto v : v_range) {
      add_vertex(v, *m_children.back());
    }
    return *m_children.back();
  }

  // local <-> global descriptor conversion functions
  [[nodiscard]] vertex_descriptor local_to_global(vertex_descriptor u_local) const {
    return is_root() ? u_local : m_global_vertex[u_local];
  }

  [[nodiscard]] vertex_descriptor global_to_local(vertex_descriptor u_global) const {
    if (is_root()) {
      return u_global;
    }
    auto [u_local, in_subgraph] = this->find_vertex(u_global);
    assert(in_subgraph == true);
    return u_local;
  }

  [[nodiscard]] edge_descriptor local_to_global(edge_descriptor e_local) const {
    return is_root() ? e_local : m_global_edge[get(get(edge_index, m_graph), e_local)];
  }

  [[nodiscard]] edge_descriptor global_to_local(edge_descriptor e_global) const {
    return is_root() ? e_global : (*m_local_edge.find(get(get(edge_index, root().m_graph), e_global))).second;
  }

  // Is vertex u (of the root graph) contained in this subgraph?
  // If so, return the matching local vertex.
  [[nodiscard]] std::pair<vertex_descriptor, bool> find_vertex(vertex_descriptor u_global) const {
    if (is_root()) {
      return {u_global, true};
    }
    auto i = m_local_vertex.find(u_global);
    bool valid = i != m_local_vertex.end();
    return {(valid ? (*i).second : null_vertex()), valid};
  }

  // Is edge e (of the root graph) contained in this subgraph?
  // If so, return the matching local edge.
  [[nodiscard]] std::pair<edge_descriptor, bool> find_edge(edge_descriptor e_global) const {
    if (is_root()) {
      return {e_global, true};
    }
    auto i = m_local_edge.find(get(get(edge_index, root().m_graph), e_global));
    bool valid = i != m_local_edge.end();
    return {(valid ? (*i).second : edge_descriptor()), valid};
  }

  // Return the parent graph.
  [[nodiscard]] subgraph& parent() { return *m_parent; }
  [[nodiscard]] const subgraph& parent() const { return *m_parent; }

  // Return true if this is the root subgraph
  [[nodiscard]] bool is_root() const { return m_parent == 0; }

  // Return the root graph of the subgraph tree.
  [[nodiscard]] subgraph& root() { return is_root() ? *this : m_parent->root(); }
  [[nodiscard]] Graph& underlying() { return m_graph; }

  [[nodiscard]] const subgraph& root() const { return is_root() ? *this : m_parent->root(); }
  [[nodiscard]] const Graph& underlying() const { return m_graph; }

  // Return the children subgraphs of this graph/subgraph.
  [[nodiscard]] auto children() const {
    return std::ranges::ref_view(m_children) |
           std::views::transform([](const auto& sg_ptr) { return std::as_const(*sg_ptr); });
  }

  [[nodiscard]] auto children() {
    return std::ranges::ref_view(m_children) | std::views::transform([](const auto& sg_ptr) { return *sg_ptr; });
  }

  [[nodiscard]] std::size_t num_children() const { return m_children.size(); }

  // Default property access delegates the lookup to global properties.
  template <typename Descriptor>
  [[nodiscard]] auto& operator[](Descriptor x) {
    return is_root() ? m_graph[x] : root().m_graph[local_to_global(x)];
  }

  template <typename Descriptor>
  [[nodiscard]] const auto& operator[](Descriptor x) const {
    return is_root() ? m_graph[x] : root().m_graph[local_to_global(x)];
  }

  // Local property access returns the local property of the given descripor.
  template <typename Descriptor>
  [[nodiscard]] auto& operator[](local_property<Descriptor> x) {
    return m_graph[x.value];
  }

  template <typename Descriptor>
  [[nodiscard]] const auto& operator[](local_property<Descriptor> x) const {
    return m_graph[x.value];
  }

  // Global property access returns the global property associated with the
  // given descriptor. This is an alias for the default bundled property
  // access operations.
  template <typename Descriptor>
  [[nodiscard]] auto& operator[](global_property<Descriptor> x) {
    return (*this)[x.value];
  }

  template <typename Descriptor>
  [[nodiscard]] const auto& operator[](global_property<Descriptor> x) const {
    return (*this)[x.value];
  }

  //  private:
  using EdgeIndexMap = property_map_t<Graph, edge_index_t>;
  using edge_index_type = property_traits_value_t<EdgeIndexMap>;
  static_assert(!std::is_same_v<edge_index_type, void>);

 private:
  using GlobalVertexList = std::vector<vertex_descriptor>;
  using GlobalEdgeList = std::vector<edge_descriptor>;
  using LocalVertexMap = std::unordered_map<vertex_descriptor, vertex_descriptor>;
  using LocalEdgeMap = std::unordered_map<edge_index_type, edge_descriptor>;

 public:  // Probably shouldn't be public....
  Graph m_graph;
  subgraph<Graph>* m_parent;
  edge_index_type m_edge_counter;  // for generating unique edge indices
  ChildrenList m_children;
  GlobalVertexList m_global_vertex;  // local -> global
  LocalVertexMap m_local_vertex;     // global -> local
  GlobalEdgeList m_global_edge;      // local -> global
  LocalEdgeMap m_local_edge;         // global -> local

  edge_descriptor local_add_edge(vertex_descriptor u_local, vertex_descriptor v_local, edge_descriptor e_global) {
    auto [e_local, inserted] = add_edge(u_local, v_local, m_graph);
    put(edge_index, m_graph, e_local, m_edge_counter++);
    m_global_edge.push_back(e_global);
    m_local_edge[get(get(edge_index, this->root()), e_global)] = e_local;
    return e_local;
  }
};

//===========================================================================
// Functions special to the Subgraph Class

template <typename G>
graph_vertex_descriptor_t<subgraph<G>> add_vertex(graph_vertex_descriptor_t<subgraph<G>> u_global, subgraph<G>& g) {
  assert(!g.is_root());
  auto [u_local, exists_local] = g.find_vertex(u_global);

  if (!exists_local) {
    // call recursion for parent subgraph
    if (!g.parent().is_root()) {
      add_vertex(u_global, g.parent());
    }

    u_local = add_vertex(g.m_graph);
    g.m_global_vertex.push_back(u_global);
    g.m_local_vertex[u_global] = u_local;

    subgraph<G>& r = g.root();

    // remember edge global and local maps
    for (auto e_global : out_edges(u_global, r)) {
      auto v_global = target(e_global, r);
      if (g.find_vertex(v_global).second) {
        g.local_add_edge(u_local, g.global_to_local(v_global), e_global);
      }
    }
    if (is_directed(g)) {  // not necessary for undirected graph
      for (auto v_global : vertices(r)) {
        if (v_global == u_global) {
          continue;  // don't insert self loops twice!
        }
        if (!g.find_vertex(v_global).second) {
          continue;  // not a subgraph vertex => try next one
        }
        for (auto e_global : out_edges(v_global, r)) {
          if (target(e_global, r) == u_global) {
            g.local_add_edge(g.global_to_local(v_global), u_local, e_global);
          }
        }
      }
    }
  }
  return u_local;
}

// NOTE: Descriptors are local unless otherwise noted.

//===========================================================================
// Functions required by the IncidenceGraph concept

template <typename G>
auto out_edges(graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return out_edges(v, g.m_graph);
}

template <typename G>
auto out_degree(graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return out_degree(v, g.m_graph);
}

template <typename G>
auto source(graph_edge_descriptor_t<G> e, const subgraph<G>& g) {
  return source(e, g.m_graph);
}

template <typename G>
auto target(graph_edge_descriptor_t<G> e, const subgraph<G>& g) {
  return target(e, g.m_graph);
}

//===========================================================================
// Functions required by the BidirectionalGraph concept

template <typename G>
auto in_edges(graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return in_edges(v, g.m_graph);
}

template <typename G>
auto in_degree(graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return in_degree(v, g.m_graph);
}

template <typename G>
auto degree(graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return degree(v, g.m_graph);
}

//===========================================================================
// Functions required by the AdjacencyGraph concept

template <typename G>
auto adjacent_vertices(graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return adjacent_vertices(v, g.m_graph);
}

//===========================================================================
// Functions required by the VertexListGraph concept

template <typename G>
auto vertices(const subgraph<G>& g) {
  return vertices(g.m_graph);
}

template <typename G>
auto num_vertices(const subgraph<G>& g) {
  return num_vertices(g.m_graph);
}

//===========================================================================
// Functions required by the EdgeListGraph concept

template <typename G>
auto edges(const subgraph<G>& g) {
  return edges(g.m_graph);
}

template <typename G>
auto num_edges(const subgraph<G>& g) {
  return num_edges(g.m_graph);
}

//===========================================================================
// Functions required by the AdjacencyMatrix concept

template <typename G>
auto edge(graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v, const subgraph<G>& g) {
  return edge(u, v, g.m_graph);
}

//===========================================================================
// Functions required by the MutableGraph concept

namespace detail {

template <typename Vertex, typename Edge, typename Graph>
void add_edge_recur_down(Vertex u_global, Vertex v_global, Edge e_global, subgraph<Graph>& g);

template <typename Vertex, typename Edge, typename Children, typename G>
void children_add_edge(Vertex u_global, Vertex v_global, Edge e_global, Children& c, subgraph<G>* orig) {
  for (auto i : c) {
    if (i->find_vertex(u_global).second && i->find_vertex(v_global).second) {
      add_edge_recur_down(u_global, v_global, e_global, *i, orig);
    }
  }
}

template <typename Vertex, typename Edge, typename Graph>
void add_edge_recur_down(Vertex u_global, Vertex v_global, Edge e_global, subgraph<Graph>& g, subgraph<Graph>* orig) {
  if (&g != orig) {
    // add local edge only if u_global and v_global are in subgraph g
    auto [u_local, u_in_subgraph] = g.find_vertex(u_global);
    auto [v_local, v_in_subgraph] = g.find_vertex(v_global);
    if (u_in_subgraph && v_in_subgraph) {
      g.local_add_edge(u_local, v_local, e_global);
    }
  }
  children_add_edge(u_global, v_global, e_global, g.m_children, orig);
}

template <typename Vertex, typename Graph>
std::pair<graph_edge_descriptor_t<Graph>, bool> add_edge_recur_up(Vertex u_global, Vertex v_global,
                                                                  const edge_property_type<Graph>& ep,
                                                                  subgraph<Graph>& g, subgraph<Graph>* orig) {
  if (g.is_root()) {
    auto [e_global, inserted] = add_edge(u_global, v_global, ep, g.m_graph);
    put(edge_index, g.m_graph, e_global, g.m_edge_counter++);
    g.m_global_edge.push_back(e_global);
    children_add_edge(u_global, v_global, e_global, g.m_children, orig);
    return {e_global, inserted};
  }
  return add_edge_recur_up(u_global, v_global, ep, *g.m_parent, orig);
}

}  // namespace detail

// Add an edge to the subgraph g, specified by the local vertex descriptors u
// and v. In addition, the edge will be added to any (all) other subgraphs that
// contain vertex descriptors u and v.

template <typename G>
std::pair<graph_edge_descriptor_t<G>, bool> add_edge(graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v,
                                                     const edge_property_type<G>& ep, subgraph<G>& g) {
  if (g.is_root()) {
    // u and v are really global
    return detail::add_edge_recur_up(u, v, ep, g, &g);
  }

  auto [e_global, inserted] = detail::add_edge_recur_up(g.local_to_global(u), g.local_to_global(v), ep, g, &g);
  auto e_local = g.local_add_edge(u, v, e_global);
  return {e_local, inserted};
}

template <typename G>
std::pair<graph_edge_descriptor_t<G>, bool> add_edge(graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v,
                                                     subgraph<G>& g) {
  return add_edge(u, v, edge_property_type<G>{}, g);
}

namespace subgraph_detail {
//-------------------------------------------------------------------------
// implementation of remove_edge(u,v,g)
template <typename Vertex, typename Graph>
void remove_edge_recur_down(Vertex u_global, Vertex v_global, subgraph<Graph>& g);

template <typename Vertex, typename Children>
void children_remove_edge(Vertex u_global, Vertex v_global, Children& c) {
  for (auto i : c) {
    if (i->find_vertex(u_global).second && i->find_vertex(v_global).second) {
      remove_edge_recur_down(u_global, v_global, *i);
    }
  }
}

template <typename Vertex, typename Graph>
void remove_edge_recur_down(Vertex u_global, Vertex v_global, subgraph<Graph>& g) {
  Vertex u_local = g.m_local_vertex[u_global];
  Vertex v_local = g.m_local_vertex[v_global];
  remove_edge(u_local, v_local, g.m_graph);
  children_remove_edge(u_global, v_global, g.m_children);
}

template <typename Vertex, typename Graph>
void remove_edge_recur_up(Vertex u_global, Vertex v_global, subgraph<Graph>& g) {
  if (g.is_root()) {
    remove_edge(u_global, v_global, g.m_graph);
    children_remove_edge(u_global, v_global, g.m_children);
  } else {
    remove_edge_recur_up(u_global, v_global, *g.m_parent);
  }
}

//-------------------------------------------------------------------------
// implementation of remove_edge(e,g)

template <typename G, typename Edge, typename Children>
void children_remove_edge(Edge e_global, Children& c) {
  for (auto i : c) {
    auto [e_local, found] = i->find_edge(e_global);
    if (!found) {
      continue;
    }
    children_remove_edge<G>(e_global, i->m_children);
    remove_edge(e_local, i->m_graph);
  }
}

}  // namespace subgraph_detail

template <typename G>
void remove_edge(graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v, subgraph<G>& g) {
  if (g.is_root()) {
    subgraph_detail::remove_edge_recur_up(u, v, g);
  } else {
    subgraph_detail::remove_edge_recur_up(g.local_to_global(u), g.local_to_global(v), g);
  }
}

template <typename G>
void remove_edge(typename subgraph<G>::edge_descriptor e, subgraph<G>& g) {
  auto e_global = g.local_to_global(e);
#ifndef NDEBUG
  auto [fe_e, fe_found] = g.find_edge(e_global);
  assert(fe_found && fe_e == e);
#endif                           // NDEBUG
  subgraph<G>& root = g.root();  // chase to root
  subgraph_detail::children_remove_edge<G>(e_global, root.m_children);
  remove_edge(e_global, root.m_graph);  // kick edge from root
}

// This is slow, but there may not be a good way to do it safely otherwise
template <typename Predicate, typename G>
void remove_edge_if(Predicate p, subgraph<G>& g) {
  while (true) {
    bool any_removed = false;
    for (auto e : edges(g)) {
      if (p(e)) {
        any_removed = true;
        remove_edge(e, g);
        break; /* Since iterators may be invalidated */
      }
    }
    if (!any_removed) {
      break;
    }
  }
}

template <typename G>
void clear_vertex(graph_vertex_descriptor_t<G> v, subgraph<G>& g) {
  while (true) {
    auto p_rg = out_edges(v, g);
    if (p_rg.empty()) {
      break;
    }
    remove_edge(*p_rg.begin(), g);
  }
}

namespace subgraph_detail {
template <typename G, typename... VPArgs>
auto add_vertex_recur_up(subgraph<G>& g, VPArgs&&... vp_args) {
  if (g.is_root()) {
    auto u_global = add_vertex(g.m_graph, std::forward<VPArgs>(vp_args)...);
    g.m_global_vertex.push_back(u_global);
    return u_global;
  }
  auto u_global = add_vertex_recur_up(*g.m_parent);
  auto u_local = add_vertex(g.m_graph);
  g.m_global_vertex.push_back(u_global);
  g.m_local_vertex[u_global] = u_local;
  return u_global;
}
}  // namespace subgraph_detail

template <typename G, typename... VPArgs>
auto add_vertex(subgraph<G>& g, VPArgs&&... vp_args) {
  if (g.is_root()) {
    auto u_global = add_vertex(g.m_graph, std::forward<VPArgs>(vp_args)...);
    g.m_global_vertex.push_back(u_global);
    return u_global;
  }
  auto u_global = subgraph_detail::add_vertex_recur_up(g.parent(), std::forward<VPArgs>(vp_args)...);
  auto u_local = add_vertex(g.m_graph);
  g.m_global_vertex.push_back(u_global);
  g.m_local_vertex[u_global] = u_local;
  return u_local;
}

// TODO: Under Construction (delete for overload resolution)
template <typename G>
void remove_vertex(graph_vertex_descriptor_t<G> u, subgraph<G>& g) = delete;

//===========================================================================
// Functions required by the PropertyGraph concept

// The global property map returns the global properties associated with local descriptors.
template <typename GraphPtr, typename PropertyMap, typename Tag>
class subgraph_global_property_map : public put_get_helper<subgraph_global_property_map<GraphPtr, PropertyMap, Tag>> {
  using Traits = property_traits<PropertyMap>;

 public:
  using value_type = typename Traits::value_type;
  using key_type = typename Traits::key_type;
  using reference = typename Traits::reference;

  subgraph_global_property_map() = default;

  subgraph_global_property_map(GraphPtr g, Tag tag) : m_g(g), m_tag(tag) {}

  reference operator[](key_type e) const {
    PropertyMap pmap = get(m_tag, m_g->root().m_graph);
    return m_g->is_root() ? pmap[e] : pmap[m_g->local_to_global(e)];
  }

  GraphPtr m_g;
  Tag m_tag;
};

// The local property map returns the local property associated with the local descriptors.
template <typename GraphPtr, typename PropertyMap, typename Tag>
class subgraph_local_property_map : public put_get_helper<subgraph_local_property_map<GraphPtr, PropertyMap, Tag>> {
  using Traits = property_traits<PropertyMap>;

 public:
  using value_type = typename Traits::value_type;
  using key_type = typename Traits::key_type;
  using reference = typename Traits::reference;

  using tag = Tag;
  using pmap = PropertyMap;

  subgraph_local_property_map() = default;

  subgraph_local_property_map(GraphPtr g, Tag tag) : m_g(g), m_tag(tag) {}

  reference operator[](key_type e) const {
    // Get property map on the underlying graph.
    PropertyMap pmap = get(m_tag, m_g->m_graph);
    return pmap[e];
  }

  GraphPtr m_g;
  Tag m_tag;
};

namespace subgraph_detail {
// Extract the actual tags from local or global property maps so we don't
// try to find non-properties.
template <typename P>
struct extract_lg_tag {
  using type = P;
};
template <typename P>
struct extract_lg_tag<local_property<P>> {
  using type = P;
};
template <typename P>
struct extract_lg_tag<global_property<P>> {
  using type = P;
};

// NOTE: Mysterious Property template parameter unused in both metafunction
// classes.
struct subgraph_global_pmap {
  template <class Tag, class SubGraph, class Property>
  struct bind_ {
    using Graph = typename SubGraph::graph_type;
    using SubGraphPtr = SubGraph*;
    using const_SubGraphPtr = const SubGraph*;
    using TagType = typename extract_lg_tag<Tag>::type;
    using PMap = property_map_t<Graph, TagType>;
    using const_PMap = property_map_const_t<Graph, TagType>;

   public:
    using type = subgraph_global_property_map<SubGraphPtr, PMap, TagType>;
    using const_type = subgraph_global_property_map<const_SubGraphPtr, const_PMap, TagType>;
  };
};

struct subgraph_local_pmap {
  template <class Tag, class SubGraph, class Property>
  struct bind_ {
    using Graph = typename SubGraph::graph_type;
    using SubGraphPtr = SubGraph*;
    using const_SubGraphPtr = const SubGraph*;
    using TagType = typename extract_lg_tag<Tag>::type;
    using PMap = property_map_t<Graph, TagType>;
    using const_PMap = property_map_const_t<Graph, TagType>;

   public:
    using type = subgraph_local_property_map<SubGraphPtr, PMap, TagType>;
    using const_type = subgraph_local_property_map<const_SubGraphPtr, const_PMap, TagType>;
  };
};

// These metafunctions select the corresponding metafunctions above, and
// are used by the choose_pmap metafunction below to specialize the choice
// of local/global property map. By default, we defer to the global
// property.
template <class Tag>
struct subgraph_choose_pmap_helper {
  using type = subgraph_global_pmap;
};
template <class Tag>
struct subgraph_choose_pmap_helper<local_property<Tag>> {
  using type = subgraph_local_pmap;
};
template <class Tag>
struct subgraph_choose_pmap_helper<global_property<Tag>> {
  using type = subgraph_global_pmap;
};

// As above, unless we're requesting vertex_index_t. Then it's always a
// local property map. This enables the correct translation of descriptors
// between local and global layers.
template <>
struct subgraph_choose_pmap_helper<vertex_index_t> {
  using type = subgraph_local_pmap;
};
template <>
struct subgraph_choose_pmap_helper<local_property<vertex_index_t>> {
  using type = subgraph_local_pmap;
};
template <>
struct subgraph_choose_pmap_helper<global_property<vertex_index_t>> {
  using type = subgraph_local_pmap;
};

// Determine the kind of property. If SameType<Tag, vertex_index_t>, then
// the property lookup is always local. Otherwise, the lookup is global.
// NOTE: Property parameter is basically unused.
template <class Tag, class Graph, class Property>
struct subgraph_choose_pmap {
  using Helper = typename subgraph_choose_pmap_helper<Tag>::type;
  using Bind = typename Helper::template bind_<Tag, Graph, Property>;
  using type = typename Bind::type;
  using const_type = typename Bind::const_type;
};

// Used by the vertex/edge property selectors to determine the kind(s) of
// property maps used by the property_map type generator.
struct subgraph_property_generator {
  template <class SubGraph, class Property, class Tag>
  struct bind_ {
    using Choice = subgraph_choose_pmap<Tag, SubGraph, Property>;
    using type = typename Choice::type;
    using const_type = typename Choice::const_type;
  };
};

}  // namespace subgraph_detail

template <>
struct vertex_property_selector<subgraph_tag> {
  using type = subgraph_detail::subgraph_property_generator;
};

template <>
struct edge_property_selector<subgraph_tag> {
  using type = subgraph_detail::subgraph_property_generator;
};

// ==================================================
// get(p, g), get(p, g, k), and put(p, g, k, v)
// ==================================================
template <typename G, typename Property>
auto get(Property p, subgraph<G>& g) {
  using PMap = property_map_t<subgraph<G>, Property>;
  return PMap(&g, p);
}

template <typename G, typename Property>
auto get(Property p, const subgraph<G>& g) {
  using PMap = property_map_const_t<subgraph<G>, Property>;
  return PMap(&g, p);
}

template <typename G, typename Property, typename Key>
auto get(Property p, const subgraph<G>& g, const Key& k) {
  using PMap = property_map_const_t<subgraph<G>, Property>;
  PMap pmap(&g, p);
  return pmap[k];
}

template <typename G, typename Property, typename Key, typename Value>
void put(Property p, subgraph<G>& g, const Key& k, const Value& val) {
  using PMap = property_map_t<subgraph<G>, Property>;
  PMap pmap(&g, p);
  pmap[k] = val;
}

// ==================================================
// get(global(p), g)
// NOTE: get(global(p), g, k) and put(global(p), g, k, v) not supported
// ==================================================
template <typename G, typename Property>
auto get(global_property<Property> p, subgraph<G>& g) {
  using Map = property_map_t<subgraph<G>, global_property<Property>>;
  return Map(&g, p.value);
}

template <typename G, typename Property>
auto get(global_property<Property> p, const subgraph<G>& g) {
  using Map = property_map_const_t<subgraph<G>, global_property<Property>>;
  return Map(&g, p.value);
}

// ==================================================
// get(local(p), g)
// NOTE: get(local(p), g, k) and put(local(p), g, k, v) not supported
// ==================================================
template <typename G, typename Property>
auto get(local_property<Property> p, subgraph<G>& g) {
  using Map = property_map_t<subgraph<G>, local_property<Property>>;
  return Map(&g, p.value);
}

template <typename G, typename Property>
auto get(local_property<Property> p, const subgraph<G>& g) {
  using Map = property_map_const_t<subgraph<G>, local_property<Property>>;
  return Map(&g, p.value);
}

template <typename G, typename Tag>
auto& get_property(subgraph<G>& g, Tag tag) {
  return get_property(g.m_graph, tag);
}

template <typename G, typename Tag>
const auto& get_property(const subgraph<G>& g, Tag tag) {
  return get_property(g.m_graph, tag);
}

//===========================================================================
// Miscellaneous Functions
template <typename G>
auto vertex(std::size_t n, const subgraph<G>& g) {
  return vertex(n, g.m_graph);
}

//===========================================================================
// Mutability Traits
// Just pull the mutability traits form the underlying graph. Note that this
// will probably fail (badly) for labeled graphs.
template <typename G>
struct graph_mutability_traits<subgraph<G>> {
  using category = graph_mutability_category_t<G>;
};

}  // namespace bagl

#endif  // BAGL_BAGL_SUBGRAPH_H_
