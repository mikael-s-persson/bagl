// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_UNDIRECTED_DFS_H_
#define BAGL_BAGL_UNDIRECTED_DFS_H_

#include <vector>

#include "bagl/depth_first_search.h"
#include "bagl/graph_traits.h"
#include "bagl/partial_range.h"

namespace bagl {

namespace udfs_detail {

template <concepts::IncidenceGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> VertexColorMap,
          concepts::ReadWritePropertyMap<graph_edge_descriptor_t<G>> EdgeColorMap>
requires concepts::ColorValue<property_traits_value_t<VertexColorMap>> &&
    concepts::ColorValue<property_traits_value_t<EdgeColorMap>>
void undir_dfv_impl(const G& g, graph_vertex_descriptor_t<G> start, V& vis, VertexColorMap vertex_color,
                    EdgeColorMap edge_color) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using Edge = graph_edge_descriptor_t<G>;
  using ColorValue = property_traits_value_t<VertexColorMap>;
  using EColorValue = property_traits_value_t<EdgeColorMap>;
  using Color = color_traits<ColorValue>;
  using EColor = color_traits<EColorValue>;
  using OutEdgeRange = graph_out_edge_range_t<G>;
  using VertexInfo = std::tuple<Vertex, std::optional<Edge>, partial_view<OutEdgeRange>>;

  std::vector<VertexInfo> stack;

  put(vertex_color, start, Color::gray());
  vis.discover_vertex(start, g);
  stack.emplace_back(start, std::nullopt, partial_view(out_edges(start, g)));
  while (!stack.empty()) {
    auto [u, u_src_e, u_er] = std::move(stack.back());
    stack.pop_back();
    while (!u_er.empty()) {
      auto u_ei = u_er.begin();
      Vertex v = target(*u_ei, g);
      vis.examine_edge(*u_ei, g);
      ColorValue v_color = get(vertex_color, v);
      EColorValue uv_color = get(edge_color, *u_ei);
      put(edge_color, *u_ei, EColor::black());
      if (v_color == Color::white()) {
        vis.tree_edge(*u_ei, g);
        u_src_e = *u_ei;
        // Re-stack the remaining edges of u.
        u_er.move_begin_to(std::next(u_ei));
        stack.emplace_back(u, u_src_e, u_er);
        // v becomes the new u.
        u = v;
        put(vertex_color, u, Color::gray());
        vis.discover_vertex(u, g);
        u_er = partial_view(out_edges(u, g));
      } else if (v_color == Color::gray()) {
        if (uv_color == EColor::white()) {
          vis.back_edge(*u_ei, g);
        }
        visitors_detail::invoke_finish_edge(vis, *u_ei, g);
        u_er.move_begin_to(std::next(u_ei));
      } else {  // if (v_color == Color::black())
        visitors_detail::invoke_finish_edge(vis, *u_ei, g);
        u_er.move_begin_to(std::next(u_ei));
      }
    }
    put(vertex_color, u, Color::black());
    vis.finish_vertex(u, g);
    if (u_src_e.has_value()) {
      visitors_detail::invoke_finish_edge(vis, u_src_e.value(), g);
    }
  }
}

}  // namespace udfs_detail

template <concepts::IncidenceGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> VertexColorMap,
          concepts::ReadWritePropertyMap<graph_edge_descriptor_t<G>> EdgeColorMap>
void undirected_depth_first_visit(const G& g, graph_vertex_descriptor_t<G> u, V vis, VertexColorMap vcolor,
                                  EdgeColorMap ecolor) {
  udfs_detail::undir_dfv_impl(g, u, vis, vcolor, ecolor);
}

template <concepts::VertexAndEdgeListGraph G, concepts::DFSVisitor<G> V,
          concepts::ReadWritePropertyMap<graph_vertex_descriptor_t<G>> VertexColorMap,
          concepts::ReadWritePropertyMap<graph_edge_descriptor_t<G>> EdgeColorMap>
void undirected_dfs(const G& g, V vis, VertexColorMap vcolor, EdgeColorMap ecolor, graph_vertex_descriptor_t<G> start) {
  using ColorValue = property_traits_value_t<VertexColorMap>;
  using Color = color_traits<ColorValue>;

  for (auto u : vertices(g)) {
    put(vcolor, u, Color::white());
    vis.initialize_vertex(u, g);
  }
  for (auto e : edges(g)) {
    put(ecolor, e, Color::white());
  }

  if (start != *vertices(g).begin()) {
    vis.start_vertex(start, g);
    udfs_detail::undir_dfv_impl(g, start, vis, vcolor, ecolor);
  }

  for (auto u : vertices(g)) {
    ColorValue u_color = get(vcolor, u);
    if (u_color == Color::white()) {
      vis.start_vertex(u, g);
      udfs_detail::undir_dfv_impl(g, u, vis, vcolor, ecolor);
    }
  }
}

template <typename Graph, typename DFSVisitor, typename VertexColorMap, typename EdgeColorMap>
void undirected_dfs(const Graph& g, DFSVisitor vis, VertexColorMap vcolor, EdgeColorMap ecolor) {
  undirected_dfs(g, vis, vcolor, ecolor, *vertices(g).begin());
}

}  // namespace bagl

#endif  // BAGL_BAGL_UNDIRECTED_DFS_H_
