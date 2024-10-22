// Copyright 2001 University of Notre Dame.
// Author: Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/depth_first_search.h"

#include <cstdlib>
#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(DepthFirstSearchConcepts, Visitor) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag>;
  using vlist_graph_t = vertex_list_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag, incidence_graph_t>;
  using graph_t = property_graph_archetype<vlist_graph_t, vertex_index_t, std::size_t>;
  graph_t g;
  dfs_visitor<> v;
  depth_first_search(g, v);
}

TEST(DepthFirstSearchConcepts, ColorMap) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag>;
  using graph_t = vertex_list_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag, incidence_graph_t>;
  graph_t g;
  dfs_visitor<> v;
  read_write_property_map_archetype<vertex_t, color_value_archetype> color;
  depth_first_search(g, v, color);
}

template <typename ColorMap, typename ParentMap, typename DiscoverTimeMap, typename FinishTimeMap>
class dfs_test_visitor {
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

 public:
  dfs_test_visitor(ColorMap color, ParentMap p, DiscoverTimeMap d, FinishTimeMap f)
      : m_color(color), m_parent(p), m_discover_time(d), m_finish_time(f), m_time(0) {}

  template <class Vertex, class Graph>
  void initialize_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(m_color, u), Color::white());
  }
  template <class Vertex, class Graph>
  void start_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(m_color, u), Color::white());
  }
  template <class Vertex, class Graph>
  void discover_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(m_color, u), Color::gray());
    EXPECT_EQ(get(m_color, get(m_parent, u)), Color::gray());

    put(m_discover_time, u, m_time++);
  }
  template <class Edge, class Graph>
  void examine_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(m_color, source(e, g)), Color::gray());
  }
  template <class Edge, class Graph>
  void tree_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(m_color, target(e, g)), Color::white());

    put(m_parent, target(e, g), source(e, g));
  }
  template <class Edge, class Graph>
  void back_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(m_color, target(e, g)), Color::gray());
  }
  template <class Edge, class Graph>
  void forward_or_cross_edge(Edge e, Graph& g) {
    EXPECT_EQ(get(m_color, target(e, g)), Color::black());
  }
  template <class Edge, class Graph>
  void finish_edge(Edge e, Graph& g) {
    EXPECT_THAT(get(m_color, target(e, g)), ::testing::AnyOf(Color::gray(), Color::black()));
  }
  template <class Vertex, class Graph>
  void finish_vertex(Vertex u, Graph&) {
    EXPECT_EQ(get(m_color, u), Color::black());

    put(m_finish_time, u, m_time++);
  }

 private:
  ColorMap m_color;
  ParentMap m_parent;
  DiscoverTimeMap m_discover_time;
  FinishTimeMap m_finish_time;
  property_traits_value_t<DiscoverTimeMap> m_time;
};

template <typename T>
class DepthFirstSearchTest : public ::testing::Test {
  static_assert(concepts::IncidenceGraph<T>);
  static_assert(concepts::VertexListGraph<T>);
  static_assert(concepts::MutableGraph<T>);
};
TYPED_TEST_SUITE_P(DepthFirstSearchTest);

TYPED_TEST_P(DepthFirstSearchTest, Search) {
  using Graph = TypeParam;
  using ColorMap = property_map_t<Graph, vertex_color_t>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  std::mt19937 gen(42);
  const std::size_t max_v = 7;

  for (std::size_t i = 0; i < max_v; ++i) {
    for (std::size_t j = 0; j < i * i; ++j) {
      Graph g;
      generate_random_graph(g, i, j, gen);

      // Get vertex index map
      auto idx = get(vertex_index, g);

      ColorMap color = get(vertex_color, g);
      auto parent_pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
      for (std::size_t k = 0; k < num_vertices(g); ++k) {
        parent_pm[k] = k;
      }
      auto discover_time_pm = vector_property_map(num_vertices(g), idx, 0);
      auto finish_time_pm = vector_property_map(num_vertices(g), idx, 0);

      dfs_test_visitor vis(color, parent_pm.ref(), discover_time_pm.ref(), finish_time_pm.ref());

      depth_first_search(g, vis, color);

      // all vertices should be black
      for (auto v : vertices(g)) {
        EXPECT_EQ(get(color, v), Color::black());
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

REGISTER_TYPED_TEST_SUITE_P(DepthFirstSearchTest, Search);

using DirectedGraph = adjacency_list<vec_s, vec_s, directed_s, property<vertex_color_t, default_color_type> >;
using UndirectedGraph = adjacency_list<vec_s, vec_s, undirected_s, property<vertex_color_t, default_color_type> >;
using TestGraphs = ::testing::Types<DirectedGraph, UndirectedGraph>;

INSTANTIATE_TYPED_TEST_SUITE_P(DepthFirstSearch, DepthFirstSearchTest, TestGraphs);

}  // namespace
}  // namespace bagl
