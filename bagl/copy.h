// Copyright 1997-2001 University of Notre Dame.
// Authors: Jeremy G. Siek, Lie-Quan Lee, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

/*
  This file implements the following functions (and overloads):

  template <typename GIn, typename GOut>
  void copy_graph(const GIn& g_in, GOut& g_out)

  template <typename GIn, typename GOut>
  graph_vertex_descriptor_t<GOut>
  copy_component(const GIn& g_in, graph_vertex_descriptor_t<GIn> src, GOut& g_out)
 */

#ifndef BAGL_BAGL_COPY_H_
#define BAGL_BAGL_COPY_H_

#include <type_traits>
#include <vector>

#include "bagl/breadth_first_search.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

namespace bagl {

namespace copy_detail {

// Copy all the vertices and edges of graph g_in into graph g_out.
// The copy_vertex and copy_edge function objects control how vertex
// and edge properties are copied.

template <int Version>
struct copy_graph_impl {};

template <>
struct copy_graph_impl<0> {
  template <typename Graph, typename MutableGraph, typename CopyVertex, typename CopyEdge, typename IndexMap,
            typename Orig2CopyVertexIndexMap>
  static void apply(const Graph& g_in, MutableGraph& g_out, CopyVertex copy_vertex, CopyEdge copy_edge,
                    Orig2CopyVertexIndexMap orig2copy, IndexMap /*unused*/) {
    for (auto v : vertices(g_in)) {
      auto new_v = add_vertex(g_out);
      put(orig2copy, v, new_v);
      copy_vertex(v, new_v);
    }
    for (auto e : edges(g_in)) {
      auto [new_e, inserted] = add_edge(get(orig2copy, source(e, g_in)), get(orig2copy, target(e, g_in)), g_out);
      copy_edge(e, new_e);
    }
  }
};

// for directed graphs
template <>
struct copy_graph_impl<1> {
  template <typename Graph, typename MutableGraph, typename CopyVertex, typename CopyEdge, typename IndexMap,
            typename Orig2CopyVertexIndexMap>
  static void apply(const Graph& g_in, MutableGraph& g_out, CopyVertex copy_vertex, CopyEdge copy_edge,
                    Orig2CopyVertexIndexMap orig2copy, IndexMap /*unused*/) {
    for (auto v : vertices(g_in)) {
      auto new_v = add_vertex(g_out);
      put(orig2copy, v, new_v);
      copy_vertex(v, new_v);
    }
    for (auto v : vertices(g_in)) {
      for (auto e : out_edges(v, g_in)) {
        auto [new_e, inserted] = add_edge(get(orig2copy, source(e, g_in)), get(orig2copy, target(e, g_in)), g_out);
        copy_edge(e, new_e);
      }
    }
  }
};

// for undirected graphs
template <>
struct copy_graph_impl<2> {
  template <typename Graph, typename MutableGraph, typename CopyVertex, typename CopyEdge, typename IndexMap,
            typename Orig2CopyVertexIndexMap>
  static void apply(const Graph& g_in, MutableGraph& g_out, CopyVertex copy_vertex, CopyEdge copy_edge,
                    Orig2CopyVertexIndexMap orig2copy, IndexMap index_map) {
    using Color = color_traits<default_color_type>;
    std::vector<default_color_type> color(num_vertices(g_in), Color::white());
    for (auto v : vertices(g_in)) {
      auto new_v = add_vertex(g_out);
      put(orig2copy, v, new_v);
      copy_vertex(v, new_v);
    }
    for (auto v : vertices(g_in)) {
      for (auto e : out_edges(v, g_in)) {
        if (color[get(index_map, target(e, g_in))] == Color::white()) {
          auto [new_e, inserted] = add_edge(get(orig2copy, source(e, g_in)), get(orig2copy, target(e, g_in)), g_out);
          copy_edge(e, new_e);
        }
      }
      color[get(index_map, v)] = Color::black();
    }
  }
};

template <typename Graph>
constexpr int copy_graph_algo_number = (is_vertex_list_graph_v<Graph> && is_edge_list_graph_v<Graph>) ? 0
                                       : is_directed_graph_v<Graph>                                   ? 1
                                                                                                      : 2;

template <class Graph>
using copy_graph_algo = copy_graph_impl<copy_graph_algo_number<Graph>>;

}  // namespace copy_detail

// Default edge and vertex property copiers

template <typename Graph1, typename Graph2>
struct default_edge_copier {
  default_edge_copier(const Graph1& g1, Graph2& g2)
      : edge_all_map1(get(edge_all, g1)), edge_all_map2(get(edge_all, g2)) {}

  template <typename Edge1, typename Edge2>
  void operator()(const Edge1& e1, Edge2& e2) const {
    put(edge_all_map2, e2, get(edge_all_map1, e1));
  }
  property_map_const_t<Graph1, edge_all_t> edge_all_map1;
  mutable property_map_t<Graph2, edge_all_t> edge_all_map2;
};

template <typename Graph1, typename Graph2>
struct default_vertex_copier {
  default_vertex_copier(const Graph1& g1, Graph2& g2)
      : vertex_all_map1(get(vertex_all, g1)), vertex_all_map2(get(vertex_all, g2)) {}

  template <typename Vertex1, typename Vertex2>
  void operator()(const Vertex1& v1, Vertex2& v2) const {
    put(vertex_all_map2, v2, get(vertex_all_map1, v1));
  }
  property_map_const_t<Graph1, vertex_all_t> vertex_all_map1;
  mutable property_map_t<Graph2, vertex_all_t> vertex_all_map2;
};

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          std::invocable<graph_vertex_descriptor_t<GIn>, graph_vertex_descriptor_t<GOut>> VertexCopier,
          std::invocable<graph_edge_descriptor_t<GIn>, graph_edge_descriptor_t<GOut>> EdgeCopier,
          concepts::ReadWriteVertexPropertyMap<GIn> OrigToCopyMap,
          concepts::ReadableVertexIndexMap<GIn> OrigVertexIndexMap>
void copy_graph(const GIn& g_in, GOut& g_out, VertexCopier v_copier, EdgeCopier e_copier, OrigToCopyMap orig_to_copy,
                OrigVertexIndexMap v_index) {
  if (num_vertices(g_in) == 0) {
    return;
  }
  using copy_impl = copy_detail::copy_graph_algo<GIn>;
  copy_impl::apply(g_in, g_out, v_copier, e_copier, orig_to_copy, v_index);
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          std::invocable<graph_vertex_descriptor_t<GIn>, graph_vertex_descriptor_t<GOut>> VertexCopier,
          std::invocable<graph_edge_descriptor_t<GIn>, graph_edge_descriptor_t<GOut>> EdgeCopier>
void copy_graph(const GIn& g_in, GOut& g_out, VertexCopier v_copier, EdgeCopier e_copier) {
  copy_graph(g_in, g_out, v_copier, e_copier,
             vector_property_map(num_vertices(g_in), get(vertex_index, g_in), graph_traits<GOut>::null_vertex()).ref(),
             get(vertex_index, g_in));
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          concepts::ReadWriteVertexPropertyMap<GIn> OrigToCopyMap,
          concepts::ReadableVertexIndexMap<GIn> OrigVertexIndexMap>
void copy_graph(const GIn& g_in, GOut& g_out, OrigToCopyMap orig_to_copy, OrigVertexIndexMap v_index) {
  copy_graph(g_in, g_out, default_vertex_copier(g_in, g_out), default_edge_copier(g_in, g_out), orig_to_copy, v_index);
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut>
void copy_graph(const GIn& g_in, GOut& g_out) {
  copy_graph(g_in, g_out, default_vertex_copier(g_in, g_out), default_edge_copier(g_in, g_out),
             vector_property_map(num_vertices(g_in), get(vertex_index, g_in), graph_traits<GOut>::null_vertex()).ref(),
             get(vertex_index, g_in));
}

namespace copy_detail {

template <class NewGraph, class CopyToOrigIndexMap, class CopyVertex, class CopyEdge>
struct graph_copy_visitor : public bfs_visitor<> {
  graph_copy_visitor(NewGraph& graph, CopyToOrigIndexMap c, CopyVertex cv, CopyEdge ce)
      : g_out_(graph), orig_to_copy_(c), copy_vertex_(cv), copy_edge_(ce) {}

  template <class Vertex>
  auto copy_one_vertex(Vertex u) const {
    auto new_u = add_vertex(g_out_);
    put(orig_to_copy_, u, new_u);
    copy_vertex_(u, new_u);
    return new_u;
  }

  template <class Edge, class Graph>
  void tree_edge(Edge e, const Graph& g_in) const {
    // For a tree edge, the target vertex has not been copied yet.
    auto [new_e, inserted] = add_edge(get(orig_to_copy_, source(e, g_in)), copy_one_vertex(target(e, g_in)), g_out_);
    copy_edge_(e, new_e);
  }

  template <class Edge, class Graph>
  void non_tree_edge(Edge e, const Graph& g_in) const {
    // For a non-tree edge, the target vertex has already been copied.
    auto [new_e, inserted] = add_edge(get(orig_to_copy_, source(e, g_in)), get(orig_to_copy_, target(e, g_in)), g_out_);
    copy_edge_(e, new_e);
  }

 private:
  NewGraph& g_out_;
  CopyToOrigIndexMap orig_to_copy_;
  CopyVertex copy_vertex_;
  CopyEdge copy_edge_;
};

}  // namespace copy_detail

// Copy all the vertices and edges of graph g_in that are reachable
// from the source vertex into graph g_out. Return the vertex
// in g_out that matches the source vertex of g_in.

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          std::invocable<graph_vertex_descriptor_t<GIn>, graph_vertex_descriptor_t<GOut>> VertexCopier,
          std::invocable<graph_edge_descriptor_t<GIn>, graph_edge_descriptor_t<GOut>> EdgeCopier,
          concepts::ReadWriteVertexPropertyMap<GIn> OrigToCopyMap,
          concepts::ReadableVertexIndexMap<GIn> OrigVertexIndexMap>
graph_vertex_descriptor_t<GOut> copy_component(const GIn& g_in, graph_vertex_descriptor_t<GIn> src, GOut& g_out,
                                               VertexCopier v_copier, EdgeCopier e_copier, OrigToCopyMap orig_to_copy,
                                               OrigVertexIndexMap v_index) {
  auto vis = copy_detail::graph_copy_visitor(g_out, orig_to_copy, v_copier, e_copier);
  auto src_copy = vis.copy_one_vertex(src);
  breadth_first_search(g_in, std::ranges::single_view{src}, vis, v_index);
  return src_copy;
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          std::invocable<graph_vertex_descriptor_t<GIn>, graph_vertex_descriptor_t<GOut>> VertexCopier,
          std::invocable<graph_edge_descriptor_t<GIn>, graph_edge_descriptor_t<GOut>> EdgeCopier,
          concepts::ReadWriteVertexPropertyMap<GIn> OrigToCopyMap>
graph_vertex_descriptor_t<GOut> copy_component(const GIn& g_in, graph_vertex_descriptor_t<GIn> src, GOut& g_out,
                                               VertexCopier v_copier, EdgeCopier e_copier, OrigToCopyMap orig_to_copy) {
  return copy_component(g_in, src, g_out, v_copier, e_copier, orig_to_copy, get(vertex_index, g_in));
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          concepts::ReadWriteVertexPropertyMap<GIn> OrigToCopyMap>
graph_vertex_descriptor_t<GOut> copy_component(const GIn& g_in, graph_vertex_descriptor_t<GIn> src, GOut& g_out,
                                               OrigToCopyMap orig_to_copy) {
  return copy_component(g_in, src, g_out, default_vertex_copier(g_in, g_out), default_edge_copier(g_in, g_out),
                        orig_to_copy, get(vertex_index, g_in));
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut,
          std::invocable<graph_vertex_descriptor_t<GIn>, graph_vertex_descriptor_t<GOut>> VertexCopier,
          std::invocable<graph_edge_descriptor_t<GIn>, graph_edge_descriptor_t<GOut>> EdgeCopier>
graph_vertex_descriptor_t<GOut> copy_component(const GIn& g_in, graph_vertex_descriptor_t<GIn> src, GOut& g_out,
                                               VertexCopier v_copier, EdgeCopier e_copier) {
  return copy_component(
      g_in, src, g_out, v_copier, e_copier,
      vector_property_map(num_vertices(g_in), get(vertex_index, g_in), graph_traits<GOut>::null_vertex()).ref(),
      get(vertex_index, g_in));
}

template <concepts::VertexListGraph GIn, concepts::MutableGraph GOut>
graph_vertex_descriptor_t<GOut> copy_component(const GIn& g_in, graph_vertex_descriptor_t<GIn> src, GOut& g_out) {
  return copy_component(
      g_in, src, g_out, default_vertex_copier(g_in, g_out), default_edge_copier(g_in, g_out),
      vector_property_map(num_vertices(g_in), get(vertex_index, g_in), graph_traits<GOut>::null_vertex()).ref(),
      get(vertex_index, g_in));
}

}  // namespace bagl

#endif  // BAGL_BAGL_COPY_H_
