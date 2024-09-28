// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2003 Bruce Barr
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DEPTH_FIRST_SEARCH_H_
#define BAGL_BAGL_DEPTH_FIRST_SEARCH_H_

#include <concepts>
#include <ranges>
#include <utility>
#include <vector>
#include <tuple>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/has_trait_member.h"
#include "bagl/partial_range.h"
#include "bagl/properties.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/visitors.h"

namespace bagl {

namespace concepts {

template <typename Visitor, typename Graph>
concept DFSVisitor = std::copy_constructible<Visitor> &&
    requires(Visitor& vis, Graph& g, graph_vertex_descriptor_t<Graph> u, graph_edge_descriptor_t<Graph> e) {
  vis.initialize_vertex(u, g);
  vis.start_vertex(u, g);
  vis.discover_vertex(u, g);
  vis.examine_edge(e, g);
  vis.tree_edge(e, g);
  vis.back_edge(e, g);
  vis.forward_or_cross_edge(e, g);
  visitors_detail::invoke_finish_edge(vis, e, g); // Optional for user
  vis.finish_vertex(u, g);
};

}  // namespace concepts

namespace depth_first_search_detail {

struct nontruth2 {
  template <typename T, typename T2>
  bool operator()(const T&, const T2&) const {
    return false;
  }
};

// If the vertex u and the iterators ei and ei_end are thought of as the
// context of the algorithm, each push and pop from the stack could
// be thought of as a context shift.
// Each pass through "while (ei != ei_end)" may refer to the out-edges of
// an entirely different vertex, because the context of the algorithm
// shifts every time a white adjacent vertex is discovered.
// The corresponding context shift back from the adjacent vertex occurs
// after all of its out-edges have been examined.
//
// See https://lists.boost.org/Archives/boost/2003/06/49265.php for FAQ.

template <concepts::IncidenceGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, typename TerminatorFunc>
    requires concepts::ColorValue<property_traits_value_t<ColorMap>>
void depth_first_visit_impl(const G& g, graph_vertex_descriptor_t<G> start, V& vis, ColorMap color,
                            TerminatorFunc func = TerminatorFunc()) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using Edge = graph_edge_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  using OutEdgeRange = graph_out_edge_range_t<G>;
  using VertexInfo = std::tuple<Vertex, std::optional<Edge>, partial_view<OutEdgeRange>>;

  std::vector<VertexInfo> stack;

  // Possible optimization for vertex list graphs
  // if constexpr (is_vertex_list_graph_v<G> && std::ranges::random_access_range<graph_vertex_range_t<G>>) {
  //   stack.reserve(vertices(g).size());
  // }

  put(color, start, Color::gray());
  vis.discover_vertex(start, g);
  auto start_oe = partial_view(out_edges(start, g));
  if (func(start, g)) {
    // If this vertex terminates the search, we push empty range
    // This will run the loop once leading to "finish_vertex" on start.
    start_oe.move_begin_to_end();
    stack.emplace_back(start, std::nullopt, std::move(start_oe));
  } else {
    stack.emplace_back(start, std::nullopt, std::move(start_oe));
  }
  while (!stack.empty()) {
    auto [u, u_src_e, u_er] = std::move(stack.back());
    stack.pop_back();
    // finish_edge has to be called here, not after the
    // loop. Think of the pop as the return from a recursive call.
    if (u_src_e.has_value()) {
      visitors_detail::invoke_finish_edge(vis, u_src_e.value(), g);
    }
    while (!u_er.empty()) {
      auto u_ei = u_er.begin();
      Vertex v = target(*u_ei, g);
      vis.examine_edge(*u_ei, g);
      ColorValue v_color = get(color, v);
      if (v_color == Color::white()) {
        vis.tree_edge(*u_ei, g);
        // Re-stack the remaining edges of u.
        u_er.move_begin_to(std::next(u_ei));
        stack.emplace_back(u, *u_ei, std::move(u_er));
        // v becomes the new u.
        u = v;
        put(color, u, Color::gray());
        vis.discover_vertex(u, g);
        u_er = partial_view(out_edges(u, g));
        if (func(u, g)) {
          u_er.move_begin_to_end(); // break
        }
      } else {
        if (v_color == Color::gray()) {
          vis.back_edge(*u_ei, g);
        } else {
          vis.forward_or_cross_edge(*u_ei, g);
        }
        visitors_detail::invoke_finish_edge(vis, *u_ei, g);
        u_er.move_begin_to(std::next(u_ei));
      }
    }
    put(color, u, Color::black());
    vis.finish_vertex(u, g);
  }
}

}  // namespace depth_first_search_detail

template <concepts::VertexListGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void depth_first_search(const G& g, V vis, ColorMap color, graph_vertex_descriptor_t<G> start_vertex) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  for (Vertex u : vertices(g)) {
    put(color, u, Color::white());
    vis.initialize_vertex(u, g);
  }

  if (start_vertex != get_default_starting_vertex(g)) {
    vis.start_vertex(start_vertex, g);
    depth_first_search_detail::depth_first_visit_impl(g, start_vertex, vis, color,
                                                      depth_first_search_detail::nontruth2());
  }

  for (Vertex u : vertices(g)) {
    ColorValue u_color = get(color, u);
    if (u_color == Color::white()) {
      vis.start_vertex(u, g);
      depth_first_search_detail::depth_first_visit_impl(g, u, vis, color, depth_first_search_detail::nontruth2());
    }
  }
}

template <concepts::VertexListGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void depth_first_search(const G& g, V vis, ColorMap color) {
  if (vertices(g).empty()) {
    return;
  }

  depth_first_search(g, vis, color, get_default_starting_vertex(g));
}

template <concepts::VertexListGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void depth_first_search(const G& g, V vis) {
  if (vertices(g).empty()) {
    return;
  }

  depth_first_search(g, vis, two_bit_color_map(num_vertices(g), get(vertex_index, g)), get_default_starting_vertex(g));
}

template <typename Visitors = null_visitors>
class dfs_visitor {
 public:
  dfs_visitor() = default;
  explicit dfs_visitor(Visitors vis) : vis_(vis) {}

  template <typename Vertex, typename Graph>
  void initialize_vertex(Vertex u, const Graph& g) {
    visitors_detail::invoke_initialize_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void start_vertex(Vertex u, const Graph& g) {
    visitors_detail::invoke_start_vertex_on_all(vis_, u, g);
  }
  template <typename Vertex, typename Graph>
  void discover_vertex(Vertex u, const Graph& g) {
    visitors_detail::invoke_discover_vertex_on_all(vis_, u, g);
  }
  template <typename Edge, typename Graph>
  void examine_edge(Edge e, const Graph& g) {
    visitors_detail::invoke_examine_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph& g) {
    visitors_detail::invoke_tree_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void back_edge(Edge e, const Graph& g) {
    visitors_detail::invoke_back_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void forward_or_cross_edge(Edge e, const Graph& g) {
    visitors_detail::invoke_forward_or_cross_edge_on_all(vis_, e, g);
  }
  template <typename Edge, typename Graph>
  void finish_edge(Edge e, const Graph& g) {
    visitors_detail::invoke_finish_edge_on_all(vis_, e, g);
  }
  template <typename Vertex, typename Graph>
  void finish_vertex(Vertex u, const Graph& g) {
    visitors_detail::invoke_finish_vertex_on_all(vis_, u, g);
  }

 protected:
  Visitors vis_;
};

template <typename Visitors>
dfs_visitor<Visitors> make_dfs_visitor(Visitors vis) {
  return dfs_visitor<Visitors>(vis);
}
using default_dfs_visitor = dfs_visitor<>;

template <concepts::IncidenceGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap>
void depth_first_visit(const G& g, graph_vertex_descriptor_t<G> u, V vis, ColorMap color) {
  vis.start_vertex(u, g);
  depth_first_search_detail::depth_first_visit_impl(g, u, vis, color, depth_first_search_detail::nontruth2());
}

template <concepts::IncidenceGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWriteVertexPropertyMap<G> ColorMap, typename TerminatorFunc>
void depth_first_visit(const G& g, graph_vertex_descriptor_t<G> u, V vis, ColorMap color,
                       TerminatorFunc func = TerminatorFunc()) {
  vis.start_vertex(u, g);
  depth_first_search_detail::depth_first_visit_impl(g, u, vis, color, func);
}
}  // namespace bagl

#endif  // BAGL_BAGL_DEPTH_FIRST_SEARCH_H_
