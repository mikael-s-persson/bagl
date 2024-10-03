// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_STRONG_COMPONENTS_H_
#define BAGL_BAGL_STRONG_COMPONENTS_H_

#include <queue>
#include <stack>
#include <type_traits>
#include <vector>

#include "bagl/connected_components.h"  // for components_recorder
#include "bagl/depth_first_search.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/transpose_graph.h"
#include "bagl/zip_range.h"

namespace bagl {

//==========================================================================
// This is Tarjan's algorithm for strongly connected components
// from his paper "Depth first search and linear graph algorithms".
// It calculates the components in a single application of DFS.
// We implement the algorithm as a dfs-visitor.

namespace scc_detail {

template <typename ComponentMap, typename RootMap, typename DiscoverTime, typename Stack>
class tarjan_scc_visitor : public dfs_visitor<> {
  using comp_type = property_traits_value_t<ComponentMap>;
  using time_type = property_traits_value_t<DiscoverTime>;

 public:
  tarjan_scc_visitor(ComponentMap comp_map, RootMap r, DiscoverTime d, comp_type& c, Stack& s)
      : c_(c), comp_(comp_map), root_(r), discover_time_(d), dfs_time_(time_type()), s_(s) {}

  template <typename Graph>
  void discover_vertex(graph_vertex_descriptor_t<Graph> v, const Graph& /*unused*/) {
    put(root_, v, v);
    put(comp_, v, (std::numeric_limits<comp_type>::max)());
    put(discover_time_, v, dfs_time_++);
    s_.push(v);
  }
  template <typename Graph>
  void finish_vertex(graph_vertex_descriptor_t<Graph> v, const Graph& g) {
    auto w = graph_traits<Graph>::null_vertex();
    for (auto e : out_edges(v, g)) {
      w = target(e, g);
      if (get(comp_, w) == std::numeric_limits<comp_type>::max()) {
        put(root_, v, min_discover_time(get(root_, v), get(root_, w)));
      }
    }
    if (get(root_, v) == v) {
      do {
        w = s_.top();
        s_.pop();
        put(comp_, w, c_);
        put(root_, w, v);
      } while (w != v);
      ++c_;
    }
  }

 private:
  template <typename Vertex>
  Vertex min_discover_time(Vertex u, Vertex v) {
    return get(discover_time_, u) < get(discover_time_, v) ? u : v;
  }

  comp_type& c_;
  ComponentMap comp_;
  RootMap root_;
  DiscoverTime discover_time_;
  time_type dfs_time_;
  Stack& s_;
};

}  // namespace scc_detail

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> ComponentMap,
          concepts::ReadWriteVertexPropertyMap<Graph> RootMap, concepts::ReadWriteVertexPropertyMap<Graph> DiscoverTime>
requires std::convertible_to<property_traits_value_t<RootMap>, graph_vertex_descriptor_t<Graph>>
auto strong_components(const Graph& g,     // Input
                       ComponentMap comp,  // Output
                       // Internal record keeping
                       RootMap root, DiscoverTime discover_time) {
  static_assert(is_directed_graph_v<Graph>);
  property_traits_value_t<ComponentMap> total = {};

  std::stack<graph_vertex_descriptor_t<Graph>> s;
  depth_first_search(g, scc_detail::tarjan_scc_visitor(comp, root, discover_time, total, s));
  return total;
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> ComponentMap,
          concepts::ReadWriteVertexPropertyMap<Graph> RootMap>
requires std::convertible_to<property_traits_value_t<RootMap>, graph_vertex_descriptor_t<Graph>>
auto strong_components(const Graph& g,     // Input
                       ComponentMap comp,  // Output
                       // Internal record keeping
                       RootMap root) {
  return strong_components(g, comp, root, vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0}));
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> ComponentMap>
auto strong_components(const Graph& g,       // Input
                       ComponentMap comp) {  // Output
  return strong_components(
      g, comp, vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex()));
}

template <typename Graph, typename ComponentMap, typename ComponentLists>
void build_component_lists(const Graph& g, std::size_t num_scc, ComponentMap component_number,
                           ComponentLists& components) {
  components.resize(num_scc);
  for (auto v : vertices(g)) {
    components[component_number[v]].push_back(v);
  }
}

//==========================================================================
// This is the version of strongly connected components from
// "Intro. to Algorithms" by Cormen, Leiserson, Rivest, which was
// adapted from "Data Structure and Algorithms" by Aho, Hopcroft,
// and Ullman, who credit the algorithm to S.R. Kosaraju and M. Sharir.
// The algorithm is based on computing DFS forests the graph
// and its transpose.

// This algorithm is slower than Tarjan's by a constant factor, uses
// more memory, and puts more requirements on the graph type.

template <concepts::MutableGraph Graph, class DFSVisitor, class ComponentsMap, class DiscoverTime, class FinishTime,
          class ColorMap>
auto kosaraju_strong_components(Graph& g, ComponentsMap c, FinishTime finish_time, ColorMap color) {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  property_traits_value_t<FinishTime> time = 0;
  depth_first_search(g, make_dfs_visitor(time_stamper_on_finish_vertex(finish_time, time)), color);

  Graph g_t(num_vertices(g));
  transpose_graph(g, g_t);

  property_traits_value_t<ComponentsMap> c_count(0);
  connected_components_detail::components_recorder vis(c, c_count);

  // initialize g_t
  for (auto v : vertices(g_t)) {
    put(color, v, Color::white());
  }

  auto fl = [&finish_time](const auto& lhs, const auto& rhs) { return get(finish_time, lhs) < get(finish_time, rhs); };
  std::priority_queue<Vertex, std::vector<Vertex>, decltype(fl)> q(fl);

  for (auto [i, j] : zip_range(vertices(g_t), vertices(g))) {
    put(finish_time, i, get(finish_time, j));
    q.push(i);
  }

  while (!q.empty()) {
    Vertex u = q.top();
    q.pop();
    if (get(color, u) == Color::white()) {
      depth_first_visit(g_t, u, vis, color);
      ++c_count;
    }
  }
  return c_count;
}

}  // namespace bagl

#endif  // BAGL_BAGL_STRONG_COMPONENTS_H_
