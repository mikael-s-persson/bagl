// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_NEIGHBOR_BFS_H_
#define BAGL_BAGL_NEIGHBOR_BFS_H_

#include <concepts>
#include <queue>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"

namespace bagl {

// Neighbor Breadth First Search
// Like BFS, but traverses in-edges as well as out-edges.
// (for directed graphs only. use normal BFS for undirected graphs)

namespace concepts {

template <typename Visitor, typename Graph>
concept NeighborBFSVisitor = std::copy_constructible<Visitor> &&
    requires(Visitor& vis, Graph& g, graph_vertex_descriptor_t<Graph> u, graph_edge_descriptor_t<Graph> e) {
  vis.initialize_vertex(u, g);
  vis.discover_vertex(u, g);
  vis.examine_vertex(u, g);
  vis.examine_out_edge(e, g);
  vis.examine_in_edge(e, g);
  vis.tree_out_edge(e, g);
  vis.tree_in_edge(e, g);
  vis.non_tree_out_edge(e, g);
  vis.non_tree_in_edge(e, g);
  vis.gray_target(e, g);
  vis.black_target(e, g);
  vis.gray_source(e, g);
  vis.black_source(e, g);
  vis.finish_vertex(u, g);
};

}  // namespace concepts

template <typename Visitors = null_visitor>
class neighbor_bfs_visitor {
 public:
  neighbor_bfs_visitor() = default;
  explicit neighbor_bfs_visitor(Visitors vis) : vis_(vis) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_initialize_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_discover_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void examine_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_examine_vertex_on_all(vis_, u, g);
  }
  template <typename Edge, typename Graph>
  void examine_out_edge(Edge e, Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void tree_out_edge(Edge e, Graph& g) {
    visitors_detail::invoke_tree_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void non_tree_out_edge(Edge e, Graph& g) {
    visitors_detail::invoke_non_tree_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void gray_target(Edge e, Graph& g) {
    visitors_detail::invoke_gray_target_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void black_target(Edge e, Graph& g) {
    visitors_detail::invoke_black_target_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void examine_in_edge(Edge e, Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void tree_in_edge(Edge e, Graph& g) {
    visitors_detail::invoke_tree_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void non_tree_in_edge(Edge e, Graph& g) {
    visitors_detail::invoke_non_tree_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void gray_source(Edge e, Graph& g) {
    visitors_detail::invoke_gray_target_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void black_source(Edge e, Graph& g) {
    visitors_detail::invoke_black_target_on_all(vis_, e, g);
  }
  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_finish_vertex_on_all(vis_, u, g);
  }

 protected:
  Visitors vis_;
};

template <typename Visitors>
auto make_neighbor_bfs_visitor(Visitors vis) {
  return neighbor_bfs_visitor<Visitors>(vis);
}

namespace neighbor_bfs_detail {

template <concepts::BidirectionalGraph G, class Buffer, concepts::NeighborBFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void neighbor_bfs_impl(const G& g, graph_vertex_descriptor_t<G> start, Buffer& q, V vis, ColorMap color) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using Edge = graph_edge_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  put(color, start, Color::gray());
  vis.discover_vertex(start, g);
  q.push(start);
  while (!q.empty()) {
    Vertex u = q.top();
    q.pop();  // pop before push to avoid problem if Q is priority_queue.
    vis.examine_vertex(u, g);

    for (auto e : out_edges(u, g)) {
      vis.examine_out_edge(e, g);
      Vertex v = target(e, g);
      ColorValue v_color = get(color, v);
      if (v_color == Color::white()) {
        vis.tree_out_edge(e, g);
        put(color, v, Color::gray());
        vis.discover_vertex(v, g);
        q.push(v);
      } else {
        vis.non_tree_out_edge(e, g);
        if (v_color == Color::gray()) {
          vis.gray_target(e, g);
        } else {
          vis.black_target(e, g);
        }
      }
    }  // for out-edges

    for (auto e : in_edges(u, g)) {
      vis.examine_in_edge(e, g);
      Vertex v = source(e, g);
      ColorValue v_color = get(color, v);
      if (v_color == Color::white()) {
        vis.tree_in_edge(e, g);
        put(color, v, Color::gray());
        vis.discover_vertex(v, g);
        q.push(v);
      } else {
        vis.non_tree_in_edge(e, g);
        if (v_color == Color::gray()) {
          vis.gray_source(e, g);
        } else {
          vis.black_source(e, g);
        }
      }
    }  // for in-edges

    put(color, u, Color::black());
    vis.finish_vertex(u, g);
  }  // while
}

}  // namespace neighbor_bfs_detail
template <concepts::BidirectionalGraph G, concepts::NeighborBFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
requires concepts::VertexListGraph<G>
void neighbor_breadth_first_search(G& g, graph_vertex_descriptor_t<G> start, V vis, ColorMap color) {
  // Buffer default
  using Vertex = graph_vertex_descriptor_t<G>;
  std::queue<Vertex> q;
  // Initialization
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  for (auto v : vertices(g)) {
    put(color, v, Color::white());
    vis.initialize_vertex(v, g);
  }
  neighbor_bfs_detail::neighbor_bfs_impl(g, start, std::ref(q), vis, color);
}

template <concepts::BidirectionalGraph G, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
requires concepts::VertexListGraph<G>
void neighbor_breadth_first_search(G& g, graph_vertex_descriptor_t<G> start, ColorMap color) {
  neighbor_breadth_first_search(g, start, make_neighbor_bfs_visitor(null_visitor()), color);
}

template <concepts::BidirectionalGraph G, concepts::NeighborBFSVisitor<G> V>
requires concepts::VertexListGraph<G>
void neighbor_breadth_first_search(G& g, graph_vertex_descriptor_t<G> start, V vis) {
  neighbor_breadth_first_search(g, start, vis,
                                make_vector_property_map<default_color_type>(num_vertices(g), get(vertex_index, g)));
}

template <concepts::BidirectionalGraph G>
requires concepts::VertexListGraph<G>
void neighbor_breadth_first_search(G& g, graph_vertex_descriptor_t<G> start) {
  neighbor_breadth_first_search(g, start, make_neighbor_bfs_visitor(null_visitor()),
                                make_vector_property_map<default_color_type>(num_vertices(g), get(vertex_index, g)));
}
}  // namespace bagl

#endif  // BAGL_BAGL_NEIGHBOR_BFS_H_
