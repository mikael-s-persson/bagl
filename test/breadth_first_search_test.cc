// Copyright 2001 University of Notre Dame.
// Author: Andrew Janiszewski, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/breadth_first_search.h"

#include <random>
#include <ranges>

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

TEST(BreadthFirstSearchConcepts, Visitor) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, undirected_tag, allow_parallel_edge_tag>;
  using vlist_graph_t =
      vertex_list_graph_archetype<vertex_t, undirected_tag, allow_parallel_edge_tag, incidence_graph_t>;
  using graph_t = property_graph_archetype<vlist_graph_t, vertex_index_t, std::size_t>;
  graph_t g;
  bfs_visitor<> v;
  breadth_first_search(g, std::ranges::empty_view<vertex_t>(), v);
}

TEST(BreadthFirstSearchConcepts, ColorMap) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag>;
  using graph_t = vertex_list_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag, incidence_graph_t>;
  graph_t g;
  bfs_visitor<> v;
  read_write_property_map_archetype<vertex_t, color_value_archetype> color;
  breadth_first_search(g, std::ranges::empty_view<vertex_t>(), v, color);
}

TEST(BreadthFirstSearchConcepts, Buffer) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag>;
  using graph_t = vertex_list_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag, incidence_graph_t>;
  graph_t g;
  bfs_visitor<> v;
  buffer_archetype<vertex_t> b;
  read_write_property_map_archetype<vertex_t, color_value_archetype> color;
  breadth_first_search(g, std::ranges::empty_view<vertex_t>(), b, v, color);
}

TEST(BreadthFirstSearchConcepts, VertexIndexMap) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag>;
  using graph_t = vertex_list_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag, incidence_graph_t>;
  graph_t g;
  bfs_visitor<> v;
  readable_property_map_archetype<vertex_t, std::size_t> v_index;
  breadth_first_search(g, std::ranges::empty_view<vertex_t>(), v, v_index);
}

template <typename Vertex, typename DistanceMap, typename ParentMap, typename ColorMap>
class bfs_testing_visitor {
  using Color = color_traits<property_traits_value_t<ColorMap>>;

 public:
  bfs_testing_visitor(Vertex s, DistanceMap d, ParentMap p, ColorMap c)
      : current_distance(0), distance(d), parent(p), color(c), src(s) {}

  template <typename Graph>
  void initialize_vertex(const Vertex& u, const Graph&) const {
    EXPECT_EQ(get(color, u), Color::white());
  }
  template <typename Graph>
  void examine_vertex(const Vertex& u, const Graph&) const {
    current_vertex = u;
    // Ensure that the distances monotonically increase.
    EXPECT_THAT(distance[u], ::testing::AnyOf(current_distance, current_distance + 1));
    if (distance[u] == current_distance + 1) {
      ++current_distance;
    }
  }
  template <typename Graph>
  void discover_vertex(const Vertex& u, const Graph&) const {
    EXPECT_EQ(get(color, u), Color::gray());
    if (u == src) {
      current_vertex = src;
    } else {
      EXPECT_EQ(parent[u], current_vertex);
      EXPECT_EQ(distance[u], current_distance + 1);
      EXPECT_EQ(distance[u], distance[parent[u]] + 1);
    }
  }
  template <typename Graph>
  void examine_edge(const graph_edge_descriptor_t<Graph>& e, const Graph& g) const {
    EXPECT_EQ(source(e, g), current_vertex);
  }
  template <typename Graph>
  void tree_edge(const graph_edge_descriptor_t<Graph>& e, const Graph& g) const {
    auto u = source(e, g);
    auto v = target(e, g);
    EXPECT_EQ(get(color, target(e, g)), Color::white());
    EXPECT_EQ(distance[u], current_distance);
    parent[v] = u;
    distance[v] = distance[u] + 1;
  }
  template <typename Graph>
  void non_tree_edge(const graph_edge_descriptor_t<Graph>& e, const Graph& g) const {
    auto u = source(e, g);
    auto v = target(e, g);
    EXPECT_NE(color[v], Color::white());

    if constexpr (is_directed_graph_v<Graph>) {
      // cross or back edge
      EXPECT_LE(distance[v], distance[u] + 1);
    } else {
      // cross edge (or going backwards on a tree edge)
      EXPECT_THAT(distance[v], ::testing::AnyOf(distance[u], distance[u] + 1, distance[u] - 1));
    }
  }

  template <typename Graph>
  void gray_target(const graph_edge_descriptor_t<Graph>& e, const Graph& g) const {
    EXPECT_EQ(color[target(e, g)], Color::gray());
  }

  template <typename Graph>
  void black_target(const graph_edge_descriptor_t<Graph>& e, const Graph& g) const {
    EXPECT_EQ(color[target(e, g)], Color::black());

    // All vertices adjacent to a black vertex must already be discovered
    for (auto u : adjacent_vertices(target(e, g), g)) {
      EXPECT_NE(color[u], Color::white());
    }
  }
  template <typename Graph>
  void finish_vertex(const Vertex& u, const Graph&) const {
    EXPECT_EQ(color[u], Color::black());
  }

 private:
  mutable Vertex current_vertex;
  mutable property_traits_value_t<DistanceMap> current_distance;
  DistanceMap distance;
  ParentMap parent;
  ColorMap color;
  Vertex src;
};

template <typename T>
class BreadthFirstSearchTest : public ::testing::Test {
  static_assert(concepts::IncidenceGraph<T>);
  static_assert(concepts::VertexListGraph<T>);
  static_assert(concepts::MutableGraph<T>);
};
TYPED_TEST_SUITE_P(BreadthFirstSearchTest);

TYPED_TEST_P(BreadthFirstSearchTest, Search) {
  using Graph = TypeParam;
  using vertex_descriptor = graph_vertex_descriptor_t<Graph>;
  using Color = color_traits<default_color_type>;

  std::mt19937 gen(42);
  const std::size_t max_v = 7;

  for (std::size_t i = 0; i < max_v; ++i) {
    for (std::size_t j = 0; j < i * i; ++j) {
      Graph g;
      generate_random_graph(g, i, j, gen);

      // declare the "start" variable
      vertex_descriptor start = random_vertex(g, gen);

      // Get vertex index map
      auto idx = get(vertex_index, g);

      // Vertex properties
      auto distance_pm = vector_property_map(i, idx, std::numeric_limits<int>::max());
      distance_pm[start] = 0;
      auto parent_pm = vector_property_map(i, idx, graph_traits<Graph>::null_vertex());
      for (auto u : vertices(g)) {
        parent_pm[u] = u;
      }
      auto color_pm = vector_property_map(i, idx, Color::white());

      // Create the testing visitor.
      bfs_testing_visitor vis(start, distance_pm, parent_pm, color_pm);

      breadth_first_search(g, std::ranges::single_view(start), vis, color_pm);

      // All white vertices should be unreachable from the source.
      for (auto u : vertices(g)) {
        if (color_pm[u] == Color::white()) {
          auto color2_pm = vector_property_map(i, idx, Color::white());
          EXPECT_FALSE(is_reachable(start, u, g, color2_pm));
        }
      }

      // The shortest path to a child should be one longer than
      // shortest path to the parent.
      for (auto u : vertices(g)) {
        if (parent_pm[u] != u) {
          // *ui not the root of the bfs tree
          EXPECT_EQ(distance_pm[u], distance_pm[parent_pm[u]] + 1);
        }
      }
    }
  }
}

REGISTER_TYPED_TEST_SUITE_P(BreadthFirstSearchTest, Search);

using DirectedGraph = adjacency_list<vec_s, vec_s, directed_s>;
using UndirectedGraph = adjacency_list<vec_s, vec_s, undirected_s>;
using TestGraphs = ::testing::Types<DirectedGraph, UndirectedGraph>;

INSTANTIATE_TYPED_TEST_SUITE_P(BreadthFirstSearch, BreadthFirstSearchTest, TestGraphs);

}  // namespace
}  // namespace bagl
