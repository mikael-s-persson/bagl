// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BREADTH_FIRST_SEARCH_H_
#define BAGL_BAGL_BREADTH_FIRST_SEARCH_H_

#include <cassert>
#include <concepts>
#include <ranges>
#include <vector>
#include <queue>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/visitors.h"

namespace bagl {

namespace concepts {

template <typename Visitor, typename Graph>
concept BFSVisitor = std::copy_constructible<Visitor> &&
    requires(Visitor& vis, Graph& g, graph_vertex_descriptor_t<Graph> u, graph_edge_descriptor_t<Graph> e) {
  vis.initialize_vertex(u, g);
  vis.discover_vertex(u, g);
  vis.examine_vertex(u, g);
  vis.examine_edge(e, g);
  vis.tree_edge(e, g);
  vis.non_tree_edge(e, g);
  vis.gray_target(e, g);
  vis.black_target(e, g);
  vis.finish_vertex(u, g);
};

}  // namespace concepts

// Multiple-see version
template <concepts::IncidenceGraph G, class Buffer, concepts::BFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, std::ranges::input_range Seeds>
void breadth_first_visit(const G& g, Seeds seeds, Buffer& q, V vis,
                         ColorMap color) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  for (Vertex s : seeds) {
    put(color, s, Color::gray());
    vis.discover_vertex(s, g);
    q.push(s);
  }
  while (!q.empty()) {
    Vertex u = q.top();
    q.pop();
    vis.examine_vertex(u, g);
    for (auto e : out_edges(u, g)) {
      Vertex v = target(e, g);
      vis.examine_edge(e, g);
      ColorValue v_color = get(color, v);
      if (v_color == Color::white()) {
        vis.tree_edge(e, g);
        put(color, v, Color::gray());
        vis.discover_vertex(v, g);
        q.push(v);
      } else {
        vis.non_tree_edge(e, g);
        if (v_color == Color::gray()) {
          vis.gray_target(e, g);
        } else {
          vis.black_target(e, g);
        }
      }
    }  // end for
    put(color, u, Color::black());
    vis.finish_vertex(u, g);
  }  // end while
}  // breadth_first_visit

// Single-seed version
template <concepts::IncidenceGraph G, concepts::Buffer Buffer, concepts::BFSVisitor<G> V, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void breadth_first_visit(const G& g, graph_vertex_descriptor_t<G> s, Buffer& q, V vis, ColorMap color) {
  breadth_first_visit(g, std::ranges::single_view(s), q, vis, color);
}

template <concepts::VertexListGraph G, std::ranges::input_range Seeds, concepts::Buffer Buffer, concepts::BFSVisitor<G> V, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void breadth_first_search(const G& g, Seeds seeds, Buffer& q, V vis, ColorMap color) {
  // Initialization
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  for (auto v : vertices(g)) {
    vis.initialize_vertex(v, g);
    put(color, v, Color::white());
  }
  breadth_first_visit(g, seeds, q, vis, color);
}

template <concepts::VertexListGraph G, concepts::Buffer Buffer, concepts::BFSVisitor<G> V, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void breadth_first_search(const G& g, graph_vertex_descriptor_t<G> s, Buffer& q, V vis, ColorMap color) {
  breadth_first_search(g, std::ranges::single_view(s), q, vis, color);
}

template <concepts::VertexListGraph G, std::ranges::input_range Seeds, concepts::BFSVisitor<G> V, concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void breadth_first_search(const G& g, Seeds seeds, V vis, ColorMap color) {
  buffer_queue<graph_vertex_descriptor_t<G>> q;
  breadth_first_search(g, seeds, q, vis, color);
}

template <concepts::VertexListGraph G, std::ranges::input_range Seeds, concepts::BFSVisitor<G> V, concepts::ReadableVertexIndexMap<G> IndexMap>
void breadth_first_search(const G& g, Seeds seeds, V vis, IndexMap v_index) {
  buffer_queue<graph_vertex_descriptor_t<G>> q;
  breadth_first_search(g, seeds, q, vis, two_bit_color_map(num_vertices(g), v_index));
}

template <concepts::VertexListGraph G, std::ranges::input_range Seeds, concepts::BFSVisitor<G> V>
void breadth_first_search(const G& g, Seeds seeds, V vis) {
  buffer_queue<graph_vertex_descriptor_t<G>> q;
  breadth_first_search(g, seeds, q, vis, two_bit_color_map(num_vertices(g), get(vertex_index, g)));
}

template <class Visitors = null_visitors>
class bfs_visitor {
 public:
  bfs_visitor() = default;
  explicit bfs_visitor(Visitors vis) : vis_(vis) {}

  template <class Vertex, class Graph>
  void initialize_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_initialize_vertex_on_all(vis_, u, g);
  }

  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_discover_vertex_on_all(vis_, u, g);
  }

  template <class Vertex, class Graph>
  void examine_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_examine_vertex_on_all(vis_, u, g);
  }

  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }

  template <class Edge, class Graph>
  void tree_edge(Edge e, Graph& g) {
    visitors_detail::invoke_tree_edge_on_all(vis_, e, g);
  }

  template <class Edge, class Graph>
  void non_tree_edge(Edge e, Graph& g) {
    visitors_detail::invoke_non_tree_edge_on_all(vis_, e, g);
  }

  template <class Edge, class Graph>
  void gray_target(Edge e, Graph& g) {
    visitors_detail::invoke_gray_target_on_all(vis_, e, g);
  }

  template <class Edge, class Graph>
  void black_target(Edge e, Graph& g) {
    visitors_detail::invoke_black_target_on_all(vis_, e, g);
  }

  template <class Vertex, class Graph>
  void finish_vertex(Vertex u, Graph& g) {
    visitors_detail::invoke_finish_vertex_on_all(vis_, u, g);
  }

 protected:
  Visitors vis_;
};

template <typename... Visitors>
auto make_bfs_visitor(Visitors&&... vis) {
  if constexpr (sizeof...(Visitors) == 0) {
    return bfs_visitor<>();
  } else {
    return bfs_visitor<std::tuple<std::decay_t<Visitors>...>>(std::forward<Visitors>(vis)...);
  }
}
using default_bfs_visitor = bfs_visitor<>;

}  // namespace bagl

#endif  // BAGL_BAGL_BREADTH_FIRST_SEARCH_H_
