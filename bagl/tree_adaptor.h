// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_TREE_ADAPTOR_H_
#define BAGL_BAGL_TREE_ADAPTOR_H_

#include <queue>

#include "bagl/detail/container_generators.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/tree_traits.h"

namespace bagl {

// This library provides function templates to adapt a MutableGraph
// such that it has the mutable interface of a tree structure.

//================= Tree =================

// Returns the first vertex-descriptor as the root of the tree.
template <concepts::VertexListGraph G>
auto tree_root(const G& g) {
  auto v_rg = vertices(g);
  return (v_rg.empty()) ? graph_traits<G>::null_vertex() : *v_rg.begin();
}

// Returns the vertex iterator range for all the child-vertices of a given vertex of the tree.
template <concepts::AdjacencyGraph G>
auto children(graph_vertex_descriptor_t<G> v, const G& g) {
  return adjacent_vertices(v, g);
}

//================= BidirectionalTree =================

// Returns the parent vertex of a given vertex descriptor in the tree.
template <typename G>
graph_vertex_descriptor_t<G> parent(graph_vertex_descriptor_t<G> v, const G& g) {
  auto ie_rg = in_edges(v, g);
  return (ie_rg.empty()) ? graph_traits<G>::null_vertex() : source(*ie_rg.begin(), g);
}

//================= MutablePropertyTree / MutableTree =================

// Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
template <concepts::VertexMutableGraph G>
graph_vertex_descriptor_t<G> create_root(G& g);

// Creates a root for the tree (clears it if not empty), and assigns the given vertex-property to it.
template <concepts::VertexMutableGraph G, typename VProp>
graph_vertex_descriptor_t<G> create_root(VProp&& vp, G& g);

// Adds a child vertex to the given parent vertex, and default-initializes the properties of
// the newly created vertex and edge.
template <concepts::MutableGraph G>
auto add_child(graph_vertex_descriptor_t<G> u, G& g) {
  auto v = add_vertex(g);
  auto [e, b] = add_edge(u, v, g);
  return std::tuple(v, e, b);
}

// Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
// vertex to the given property value.
template <concepts::MutableGraph G, typename VProp>
auto add_child(graph_vertex_descriptor_t<G> u, VProp&& vp, G& g) {
  auto v = add_vertex(std::forward<VProp>(vp), g);
  auto [e, b] = add_edge(u, v, g);
  return std::tuple(v, e, b);
}

// Adds a child vertex to the given parent vertex, and initializes the properties of the newly created
// vertex and edge to the given property values.
template <concepts::MutableGraph G, typename VProp, typename EProp>
auto add_child(graph_vertex_descriptor_t<G> u, VProp&& vp, EProp&& ep, G& g) {
  auto v = add_vertex(std::forward<VProp>(vp), g);
  auto [e, b] = add_edge(u, v, std::forward<EProp>(ep), g);
  return std::tuple(v, e, b);
}

namespace tree_adaptor_detail {

template <typename G>
void record_vertex_prop_impl(graph_vertex_descriptor_t<G> /*unused*/, container_detail::ignore_output_iter& /*unused*/,
                             G& /*unused*/) {}

template <typename G, typename VertexOIter>
void record_vertex_prop_impl(graph_vertex_descriptor_t<G> u, VertexOIter& vit_out, G& g) {
  *(vit_out++) = std::move(g[u]);
}

template <typename G>
void record_edge_prop_impl(graph_edge_descriptor_t<G> /*unused*/, container_detail::ignore_output_iter& /*unused*/,
                           G& /*unused*/) {}

template <typename G, typename EdgeOIter>
void record_edge_prop_impl(graph_edge_descriptor_t<G> e, EdgeOIter& eit_out, G& g) {
  *(eit_out++) = std::move(g[e]);
}

template <typename G>
void record_in_edge_prop_impl(graph_vertex_descriptor_t<G> /*unused*/, container_detail::ignore_output_iter& /*unused*/,
                              G& /*unused*/) {}

template <typename G, typename EdgeOIter>
std::enable_if_t<is_bidirectional_graph_v<G>> record_in_edge_prop_impl(graph_vertex_descriptor_t<G> u,
                                                                       EdgeOIter& eit_out, G& g) {
  // depends on bidirectionality of the graph type:
  auto iei_rg = in_edges(u, g);
  if (!iei_rg.empty()) {
    record_edge_prop_impl(*iei_rg.begin(), eit_out, g);
  } else {
    // fill in a dummy edge-property so that number of written vertex and edge properties match.
    *(eit_out++) = edge_bundle_type<G>{};
  }
}

template <typename G, typename EdgeOIter>
std::enable_if_t<!is_bidirectional_graph_v<G>> record_in_edge_prop_impl(graph_vertex_descriptor_t<G> u,
                                                                        EdgeOIter& eit_out, G& g) {
  // fall-back solution for a uni-directional graph type:
  using Vertex = graph_vertex_descriptor_t<G>;
  std::queue<Vertex> v_queue;
  Vertex v = tree_root(g);
  if (v == u) {
    // fill in a dummy edge-property so that number of written vertex and edge properties match.
    *(eit_out++) = edge_bundle_type<G>{};
    return;
  }
  v_queue.push(v);

  while (!v_queue.empty()) {
    v = v_queue.front();
    v_queue.pop();
    for (auto e : out_edges(v, g)) {
      if (target(e, g) == u) {
        record_edge_prop_impl(e, eit_out, g);
        return;
      }
      v_queue.push(target(e, g));
    }
  }
}

}  // namespace tree_adaptor_detail

template <concepts::MutableGraph G, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children_impl(graph_vertex_descriptor_t<G> u, VertexOIter vit_out,
                                                      EdgeOIter eit_out, G& g) {
  using Vertex = graph_vertex_descriptor_t<G>;
  std::queue<Vertex> v_queue;
  v_queue.push(u);
  while (!v_queue.empty()) {
    Vertex v = v_queue.front();
    v_queue.pop();
    for (auto e : out_edges(v, g)) {
      v_queue.push(target(e, g));
      tree_adaptor_detail::record_vertex_prop_impl(target(e, g), vit_out, g);
      tree_adaptor_detail::record_edge_prop_impl(e, eit_out, g);
    }
    if (v != u) {
      clear_vertex(v, g);
      remove_vertex(v, g);
    }
  }
  return {vit_out, eit_out};
}

template <concepts::MutableGraph G, typename VertexOIter, typename EdgeOIter>
requires std::ranges::random_access_range<decltype(vertices(std::declval<G>()))> std::pair<VertexOIter, EdgeOIter>
clear_children_impl(graph_vertex_descriptor_t<G> u, VertexOIter vit_out, EdgeOIter eit_out, G& g) {
  /* TODO: Figure out a way to implement this better! */

  using Vertex = graph_vertex_descriptor_t<G>;
  using std::swap;

  G g_tmp;
  std::queue<Vertex> v_tmp_queue;
  std::queue<Vertex> v_ori_queue;
  Vertex v = tree_root(g);
  v_ori_queue.push(v);
  v_tmp_queue.push(create_root(std::move(g[v]), g_tmp));

  while (!v_ori_queue.empty()) {
    v = v_ori_queue.front();
    v_ori_queue.pop();
    Vertex v_tmp = v_tmp_queue.front();
    v_tmp_queue.pop();
    if (v == u) {
      continue;
    }
    for (auto e : out_edges(v, g)) {
      v_ori_queue.push(target(e, g));
      v_tmp_queue.push(add_child(v_tmp, std::move(g[target(e, g)]), std::move(g[e]), g_tmp).first);
    }
  }

  v_ori_queue.push(u);
  while (!v_ori_queue.empty()) {
    v = v_ori_queue.front();
    v_ori_queue.pop();
    for (auto e : out_edges(v, g)) {
      v_ori_queue.push(target(e, g));
      tree_adaptor_detail::record_vertex_prop_impl(target(e, g), vit_out, g);
      tree_adaptor_detail::record_edge_prop_impl(e, eit_out, g);
    }
  }

  swap(g, g_tmp);

  return {vit_out, eit_out};
}

// Removes a branch (sub-tree) starting from but excluding the given vertex.
template <concepts::MutableGraph G>
void clear_children(graph_vertex_descriptor_t<G> u, G& g) {
  clear_children_impl(u, container_detail::ignore_output_iter(), container_detail::ignore_output_iter(), g);
}

// Removes a branch (sub-tree) starting from but excluding the given vertex, while
// recording the vertex-properties of all the removed vertices into an output-iterator.
template <concepts::MutableGraph G, typename OutputIter>
OutputIter clear_children(graph_vertex_descriptor_t<G> u, OutputIter it_out, G& g) {
  return clear_children_impl(u, it_out, container_detail::ignore_output_iter(), g).first;
}

// Removes a branch (sub-tree) starting from but excluding the given vertex, while
// recording the vertex-properties of all the removed vertices into an output-iterator.
template <concepts::MutableGraph G, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> clear_children(graph_vertex_descriptor_t<G> u, VertexOIter vit_out, EdgeOIter eit_out,
                                                 G& g) {
  return clear_children_impl(u, vit_out, eit_out, g);
}

template <concepts::MutableGraph G, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch_impl(graph_vertex_descriptor_t<G> u, VertexOIter vit_out,
                                                     EdgeOIter eit_out, G& g) {
  using Vertex = graph_vertex_descriptor_t<G>;
  std::queue<Vertex> v_queue;
  v_queue.push(u);
  tree_adaptor_detail::record_vertex_prop_impl(u, vit_out, g);
  tree_adaptor_detail::record_in_edge_prop_impl(u, eit_out, g);

  while (!v_queue.empty()) {
    Vertex v = v_queue.front();
    v_queue.pop();
    for (auto e : out_edges(v, g)) {
      v_queue.push(target(e, g));
      tree_adaptor_detail::record_vertex_prop_impl(target(e, g), vit_out, g);
      tree_adaptor_detail::record_edge_prop_impl(e, eit_out, g);
    }
    clear_vertex(v, g);
    remove_vertex(v, g);
  }

  return {vit_out, eit_out};
}

template <concepts::MutableGraph G, typename VertexOIter, typename EdgeOIter>
requires std::ranges::random_access_range<decltype(vertices(std::declval<G>()))> std::pair<VertexOIter, EdgeOIter> remove_branch_impl(
    graph_vertex_descriptor_t<G> u, VertexOIter vit_out, EdgeOIter eit_out, G& g) {
  /* TODO: Figure out a way to implement this better! */

  using Vertex = graph_vertex_descriptor_t<G>;
  using std::swap;

  G g_tmp;
  std::queue<Vertex> v_tmp_queue;
  std::queue<Vertex> v_ori_queue;
  Vertex v = tree_root(g);
  if (v != u) {
    v_ori_queue.push(v);
    v_tmp_queue.push(create_root(std::move(g[v]), g_tmp));
  } else {
    tree_adaptor_detail::record_vertex_prop_impl(v, vit_out, g);
    // fill in a dummy edge-property so that number of written vertex and edge properties match.
    *(eit_out++) = edge_bundle_type<G>{};
  }

  while (!v_ori_queue.empty()) {
    v = v_ori_queue.front();
    v_ori_queue.pop();
    Vertex v_tmp = v_tmp_queue.front();
    v_tmp_queue.pop();
    for (auto e : out_edges(v, g)) {
      if (target(e, g) == u) {
        tree_adaptor_detail::record_vertex_prop_impl(u, vit_out, g);
        tree_adaptor_detail::record_edge_prop_impl(e, eit_out, g);
        continue;
      }
      v_ori_queue.push(target(e, g));
      v_tmp_queue.push(add_child(v_tmp, std::move(g[target(e, g)]), std::move(g[e]), g_tmp).first);
    }
  }

  v_ori_queue.push(u);
  while (!v_ori_queue.empty()) {
    v = v_ori_queue.front();
    v_ori_queue.pop();
    for (auto e : out_edges(v, g)) {
      v_ori_queue.push(target(e, g));
      tree_adaptor_detail::record_vertex_prop_impl(target(e, g), vit_out, g);
      tree_adaptor_detail::record_edge_prop_impl(e, eit_out, g);
    }
  }

  swap(g, g_tmp);

  return {vit_out, eit_out};
}

// Removes a branch (sub-tree) starting from and including the given vertex.
template <concepts::MutableGraph G>
void remove_branch(graph_vertex_descriptor_t<G> u, G& g) {
  remove_branch_impl(u, container_detail::ignore_output_iter(), container_detail::ignore_output_iter(), g);
}

// Removes a branch (sub-tree) starting from and including the given vertex, while
// recording the vertex-properties of all the removed vertices into an output-iterator.
// The first vertex-property to figure in the output range is that of the vertex v.
template <concepts::MutableGraph G, typename OutputIter>
OutputIter remove_branch(graph_vertex_descriptor_t<G> u, OutputIter it_out, G& g) {
  return remove_branch_impl(u, it_out, container_detail::ignore_output_iter(), g).first;
}

// Removes a branch (sub-tree) starting from and including the given vertex, while
// recording the vertex and edge properties of all the removed vertices and edges into output-ranges.
// The first vertex-property to figure in the output range is that of the vertex v.
template <concepts::MutableGraph G, typename VertexOIter, typename EdgeOIter>
std::pair<VertexOIter, EdgeOIter> remove_branch(graph_vertex_descriptor_t<G> u, VertexOIter vit_out, EdgeOIter eit_out,
                                                G& g) {
  return remove_branch_impl(u, vit_out, eit_out, g);
}

template <concepts::VertexMutableGraph G>
graph_vertex_descriptor_t<G> create_root(G& g) {
  auto old_root = tree_root(g);
  if (old_root != graph_traits<G>::null_vertex()) {
    remove_branch(old_root, g);
  }
  return add_vertex(g);
}

template <concepts::VertexMutableGraph G, typename VProp>
graph_vertex_descriptor_t<G> create_root(VProp&& vp, G& g) {
  auto old_root = tree_root(g);
  if (old_root != graph_traits<G>::null_vertex()) {
    remove_branch(old_root, g);
  }
  return add_vertex(std::forward<VProp>(vp), g);
}

//================= PropertyGraph =================

template <typename Graph>
vertex_bundle_type<Graph>& get_property(graph_vertex_descriptor_t<Graph> v, Graph& g) {
  return g[v];
}

template <typename Graph>
const vertex_bundle_type<Graph>& get_property(graph_vertex_descriptor_t<Graph> v, const Graph& g) {
  return g[v];
}

template <typename Graph>
edge_bundle_type<Graph>& get_property(graph_edge_descriptor_t<Graph> e, Graph& g) {
  return g[e];
}

template <typename Graph>
const edge_bundle_type<Graph>& get_property(graph_edge_descriptor_t<Graph> e, const Graph& g) {
  return g[e];
}

}  // namespace bagl

#endif  // BAGL_BAGL_TREE_ADAPTOR_H_
