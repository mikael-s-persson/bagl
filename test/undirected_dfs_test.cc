// Copyright 2001 University of Notre Dame.
// Author: Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/undirected_dfs.h"

#include <cstdlib>
#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename ColorMap, typename ParentMap, typename DiscoverTimeMap, typename FinishTimeMap>
class dfs_test_visitor {
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

 public:
  dfs_test_visitor(ColorMap color, ParentMap p, DiscoverTimeMap d, FinishTimeMap f)
      : color_(color), parent_(p), discover_time_(d), finish_time_(f), time_(0) {}

  template <class Vertex, class Graph>
  void initialize_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(color_, u), Color::white());
  }
  template <class Vertex, class Graph>
  void start_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(color_, u), Color::white());
  }
  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(color_, u), Color::gray());
    EXPECT_EQ(get(color_, get(parent_, u)), Color::gray());

    put(discover_time_, u, time_++);
  }
  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(color_, source(e, g)), Color::gray());
  }
  template <class Edge, class Graph>
  void tree_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(color_, target(e, g)), Color::white());
    put(parent_, target(e, g), source(e, g));
  }
  template <class Edge, class Graph>
  void back_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(color_, target(e, g)), Color::gray());
  }
  template <class Edge, class Graph>
  void forward_or_cross_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(color_, target(e, g)), Color::black());
  }
  template <class Edge, class Graph>
  void finish_edge(Edge e, Graph& g) {
    EXPECT_THAT(get(color_, target(e, g)), ::testing::AnyOf(Color::gray(), Color::black()));
  }
  template <class Vertex, class Graph>
  void finish_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(color_, u), Color::black());
    put(finish_time_, u, time_++);
  }

 private:
  ColorMap color_;
  ParentMap parent_;
  DiscoverTimeMap discover_time_;
  FinishTimeMap finish_time_;
  property_traits_value_t<DiscoverTimeMap> time_;
};

TEST(UndirectedDFSTest, Basic) {
  using Graph = adjacency_list<vec_s, vec_s, undirected_s, property<vertex_color_t, default_color_type>,
                               property<edge_color_t, default_color_type> >;
  using ColorMap = property_map_t<Graph, vertex_color_t>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  using EColorMap = property_map_t<Graph, edge_color_t>;
  using EColorValue = property_traits_value_t<EColorMap>;
  using EColor = color_traits<EColorValue>;

  std::mt19937 gen(42);
  const std::size_t max_v = 7;

  for (std::size_t i = 0; i < max_v; ++i) {
    for (std::size_t j = 0; j < i * i; ++j) {
      Graph g;
      generate_random_graph(g, i, j, gen);

      // Get vertex index map
      auto idx = get(vertex_index, g);

      auto color = get(vertex_color, g);
      auto e_color = get(edge_color, g);
      auto parent_pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
      for (auto v : vertices(g)) {
        parent_pm[v] = v;
      }
      auto discover_time_pm = vector_property_map(num_vertices(g), idx, int{0});
      auto finish_time_pm = vector_property_map(num_vertices(g), idx, int{0});

      dfs_test_visitor vis(color, parent_pm.ref(), discover_time_pm.ref(), finish_time_pm.ref());

      undirected_dfs(g, vis, color, e_color);

      // all vertices should be black
      for (auto v : vertices(g)) {
        EXPECT_EQ(get(color, v), Color::black());
      }

      // all edges should be black
      for (auto e : edges(g)) {
        EXPECT_EQ(get(e_color, e), EColor::black());
      }

      // check parenthesis structure of discover/finish times
      // See CLR p.480
      for (auto u : vertices(g)) {
        for (auto v : vertices(g)) {
          if (u == v) {
            continue;
          }
          EXPECT_THAT(
              true, ::testing::AnyOf(
                        finish_time_pm[u] < discover_time_pm[v], finish_time_pm[v] < discover_time_pm[u],
                        ::testing::AllOf(discover_time_pm[v] < discover_time_pm[u],
                                         finish_time_pm[u] < finish_time_pm[v], is_descendant(u, v, parent_pm.ref())),
                        ::testing::AllOf(discover_time_pm[u] < discover_time_pm[v],
                                         finish_time_pm[v] < finish_time_pm[u], is_descendant(v, u, parent_pm.ref()))));
        }
      }
    }
  }
}

}  // namespace
}  // namespace bagl
