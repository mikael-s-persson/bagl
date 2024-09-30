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

template <typename Graph, typename Node, typename Derived>
class graph_as_tree_base {
  using Tree = Derived;

 public:
  using node_descriptor = Node;

  graph_as_tree_base(Graph& g, Node root) : g_(g), root_(root) {}

  friend Node root(const Tree& t) { return t._root; }

  template <typename N>
  friend auto children(N n, const Tree& t) {
    return adjacent_vertices(n, t.g_);
  }

  template <typename N>
  friend Node parent(N n, const Tree& t) {
    return get(t.parent_pa(), n);
  }

 private:
  Graph& g_;
  Node root_;
};

struct graph_as_tree_tag {};

template <typename Graph, typename ParentMap, typename Node = graph_vertex_descriptor_t<Graph>>
class graph_as_tree : public graph_as_tree_base<Graph, Node, graph_as_tree<Graph, ParentMap, Node> > {
  using self = graph_as_tree;
  using super = graph_as_tree_base<Graph, Node, self>;

 public:
  graph_as_tree(Graph& g, Node root) : super(g, root) {}

  graph_as_tree(Graph& g, Node root, ParentMap p) : super(g, root), p_(p) {
    breadth_first_search(g, root, visitor(make_bfs_visitor(predecessor_recorder_on_tree_edge(p))));
  }
  ParentMap parent_pa() const { return p_; }
  using graph_tag = graph_as_tree_tag;  // for property_map
 protected:
  ParentMap p_;
};

template <typename Graph, typename Node, typename ParentMap>
auto make_graph_as_tree(Graph& g, Node root, ParentMap p) {
  return graph_as_tree<Graph, ParentMap, Node>(g, root, p);
}

namespace detail {

struct graph_as_tree_vertex_property_selector {
  template <typename GraphAsTree, typename Property, typename Tag>
  struct bind_ {
    using Graph = typename GraphAsTree::base_type;
    using PMap = property_map<Graph, Tag>;
    using type = typename PMap::type;
    using const_type = typename PMap::const_type;
  };
};

struct graph_as_tree_edge_property_selector {
  template <typename GraphAsTree, typename Property, typename Tag>
  struct bind_ {
    using Graph = typename GraphAsTree::base_type;
    using PMap = property_map<Graph, Tag>;
    using type = typename PMap::type;
    using const_type = typename PMap::const_type;
  };
};

}  // namespace detail

template <>
struct vertex_property_selector<graph_as_tree_tag> {
  using type = detail::graph_as_tree_vertex_property_selector;
};

template <>
struct edge_property_selector<graph_as_tree_tag> {
  using type = detail::graph_as_tree_edge_property_selector;
};

template <typename Graph, typename P, typename N, typename Property>
auto get(Property p, graph_as_tree<Graph, P, N>& g) {
  return get(p, g.g_);
}

template <typename Graph, typename P, typename N, typename Property>
auto get(Property p, const graph_as_tree<Graph, P, N>& g) {
  return get(p, g.g_);
}

template <typename Graph, typename P, typename N, typename Property, typename Key>
auto get(Property p, const graph_as_tree<Graph, P, N>& g, const Key& k) {
  return get(p, g.g_, k);
}

template <typename Graph, typename P, typename N, typename Property, typename Key, typename Value>
void put(Property p, const graph_as_tree<Graph, P, N>& g, const Key& k, const Value& val) {
  put(p, g.g_, k, val);
}

}  // namespace bagl

#endif  //  BAGL_BAGL_GRAPH_AS_TREE_H_
