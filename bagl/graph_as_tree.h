// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_GRAPH_AS_TREE_H_
#define BAGL_BAGL_GRAPH_AS_TREE_H_

#include "bagl/breadth_first_search.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/tree_traits.h"
#include "bagl/visitors.h"

namespace bagl {

template <typename Graph, typename Derived>
class graph_as_tree_base {
  using Tree = Derived;

 public:
  using node_descriptor = graph_vertex_descriptor_t<Graph>;

  graph_as_tree_base(const Graph& g, node_descriptor root) : g_(g), root_(root) {}

  node_descriptor get_root() const { return root_; }
  const Graph& get_graph() const { return g_; }

 private:
  const Graph& g_;
  node_descriptor root_;
};

template <typename Graph, typename Derived>
auto root(const graph_as_tree_base<Graph, Derived>& t) {
  return t.get_root();
}

template <typename Graph, typename Derived>
auto children(graph_vertex_descriptor_t<Graph> n, const graph_as_tree_base<Graph, Derived>& t) {
  return adjacent_vertices(n, t.get_graph());
}

template <typename Graph, typename Derived>
auto parent(graph_vertex_descriptor_t<Graph> n, const graph_as_tree_base<Graph, Derived>& t) {
  return get(static_cast<const Derived&>(t).parent_pa(), n);
}

struct graph_as_tree_tag {};

template <typename Graph, typename ParentMap>
class graph_as_tree : public graph_as_tree_base<Graph, graph_as_tree<Graph, ParentMap>> {
  using self = graph_as_tree;
  using super = graph_as_tree_base<Graph, self>;

 public:
  using node_descriptor = graph_vertex_descriptor_t<Graph>;

  graph_as_tree(const Graph& g, node_descriptor root) : super(g, root) {}

  graph_as_tree(const Graph& g, node_descriptor root, ParentMap p) : super(g, root), p_(p) {
    breadth_first_search(g, std::ranges::single_view{root}, make_bfs_visitor(predecessor_recorder_on_tree_edge(p)));
  }
  ParentMap parent_pa() const { return p_; }
  using graph_tag = graph_as_tree_tag;  // for property_map
 protected:
  ParentMap p_;
};

template <typename Graph, typename ParentMap>
auto make_graph_as_tree(const Graph& g, graph_vertex_descriptor_t<Graph> root, ParentMap p) {
  return graph_as_tree<Graph, ParentMap>(g, root, p);
}

template <typename Graph, typename P>
struct graph_traits<graph_as_tree<Graph, P>> : graph_traits<Graph> {};

template <typename Graph, typename P, typename Tag>
struct property_map<graph_as_tree<Graph, P>, Tag> : property_map<Graph, Tag> {};

template <typename Graph, typename P, typename Property>
auto get(Property p, const graph_as_tree<Graph, P>& g) {
  return get(p, g.get_graph());
}

template <typename Graph, typename P, typename Property, typename Key>
auto get(Property p, const graph_as_tree<Graph, P>& g, Key&& k) {
  return get(p, g.get_graph(), std::forward<Key>(k));
}

template <typename Graph, typename P, typename Property, typename Key, typename Value>
void put(Property p, const graph_as_tree<Graph, P>& g, Key&& k, Value&& val) {
  put(p, g.get_graph(), std::forward<Key>(k), std::forward<Value>(val));
}

}  // namespace bagl

#endif  //  BAGL_BAGL_GRAPH_AS_TREE_H_
