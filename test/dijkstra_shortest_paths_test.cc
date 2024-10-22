// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson

#include "bagl/dijkstra_shortest_paths.h"

#include <cmath>
#include <iostream>
#include <list>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/dijkstra_shortest_paths_no_color_map.h"
#include "bagl/exception.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/random.h"
#include "bagl/relax.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test/new_york_state_map.h"

namespace bagl {
namespace {

TEST(DijkstraSearchConcepts, Visitor) {
  using vertex_t = semiregular_archetype<1>;
  using incidence_graph_t = incidence_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag>;
  using vlist_graph_t = vertex_list_graph_archetype<vertex_t, directed_tag, allow_parallel_edge_tag, incidence_graph_t>;
  using graph_t = property_graph_archetype<vlist_graph_t, vertex_index_t, std::size_t>;
  using edge_t = graph_edge_descriptor_t<graph_t>;
  graph_t g;
  readable_property_map_archetype<edge_t, std::size_t> weight;
  readable_property_map_archetype<vertex_t, int> index;
  read_write_property_map_archetype<vertex_t, vertex_t> pred;
  read_write_property_map_archetype<vertex_t, std::size_t> distance;
  dijkstra_shortest_paths(g, std::ranges::empty_view<vertex_t>(), pred, distance, weight, index,
                          make_dijkstra_visitor());
  dijkstra_shortest_paths_no_color_map(g, std::ranges::empty_view<vertex_t>(), pred, distance, weight, index,
                                       make_dijkstra_visitor());
}

template <bool WithoutColorMap>
void test_dijkstra_search() {
  // specify some types
  using Graph = test::state_map_type;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::mt19937 gen(42);

  Graph g = test::make_new_york_state_map();

  auto name = get(vertex_name, g);
  auto idx = get(vertex_index, g);
  auto weight = get(edge_weight, g);

  auto pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
  auto dm = vector_property_map(num_vertices(g), idx, std::numeric_limits<float>::max());

  for (int i = 0; i < 5; ++i) {
    // pick random start/goal
    Vertex start = gen() % num_vertices(g);
    Vertex goal = gen() % num_vertices(g);

    auto print_problem_def = [&]() {
      std::stringstream out;
      print_graph(g, get(vertex_bundle, g), weight, out);
      out << "\nStart vertex: " << name[start];
      out << "\nGoal vertex: " << name[goal];
      return out.str();
    };
    // Print problem definition if any checks fail.
    SCOPED_TRACE(print_problem_def());

    // call astar named parameter interface

    if constexpr (WithoutColorMap) {
      EXPECT_THROW(dijkstra_shortest_paths_no_color_map(g, std::ranges::single_view{start}, pm.ref(), dm.ref(), weight,
                                                        idx, make_dijkstra_visitor(bfs_stop_at_goal(goal))),
                   search_succeeded)
          << "Did not find a path to the goal at all!";
    } else {
      EXPECT_THROW(dijkstra_shortest_paths(g, start, pm.ref(), dm.ref(), weight, idx,
                                           make_dijkstra_visitor(bfs_stop_at_goal(goal))),
                   search_succeeded)
          << "Did not find a path to the goal at all!";
    }
    std::vector<Vertex> shortest_path;
    float path_weight_sum = 0.0F;
    for (Vertex v = goal;; v = pm[v]) {
      shortest_path.push_back(v);
      if (pm[v] == v) {
        // Expect to end up back at start.
        EXPECT_EQ(v, start);
        break;
      }
      auto [e, e_found] = edge(pm[v], v, g);
      ASSERT_TRUE(e_found) << "Edge from " << name[pm[v]] << " to " << name[v] << " does not exist!";
      EXPECT_THAT(dm[pm[v]] + get(weight, e), ::testing::FloatNear(dm[v], 1.0e-5F));
      path_weight_sum = path_weight_sum + get(weight, e);
    }
    std::ranges::reverse(shortest_path);
    auto print_solution = [&]() {
      std::stringstream out;
      out << "Shortest path from " << name[start] << " to " << name[goal] << ": ";
      bool first_u = true;
      for (auto spi : shortest_path) {
        if (!first_u) {
          out << " -> ";
        }
        out << name[spi];
        first_u = false;
      }
      out << "\nTotal travel time: " << dm[goal];
      return out.str();
    };
    // Print solution if any checks fail.
    SCOPED_TRACE(print_solution());

    EXPECT_THAT(dm[goal], ::testing::FloatNear(path_weight_sum, 1.0e-5F));
    float early_stop_cost = dm[goal];
    // Run astar exhaustively, to be sure.
    if constexpr (WithoutColorMap) {
      dijkstra_shortest_paths_no_color_map(g, std::ranges::single_view{start}, pm.ref(), dm.ref(), weight, idx,
                                           default_dijkstra_visitor());
    } else {
      dijkstra_shortest_paths(g, start, pm.ref(), dm.ref(), weight, idx, default_dijkstra_visitor());
    }
    EXPECT_THAT(early_stop_cost, ::testing::FloatNear(dm[goal], 1.0e-5F));
    for (auto v : vertices(g)) {
      if (pm[v] == v) {
        continue;
      }
      // Recorded predecessor is the best.
      for (auto e : out_edges(v, g)) {
        EXPECT_GE(dm[target(e, g)] + get(weight, e), dm[v] - 1e-5F);
      }
    }
  }
}

TEST(DijkstraSearchTest, Search) {
  SCOPED_TRACE("While testing dijkstra_shortest_paths");
  test_dijkstra_search<false>();
}

TEST(DijkstraSearchTest, NoColorMap) {
  SCOPED_TRACE("While testing dijkstra_shortest_paths_no_color_map");
  test_dijkstra_search<true>();
}

TEST(DijkstraSearchTest, CompareAlgorithms) {
  // specify some types
  using Graph = adjacency_list<list_s, vec_s, directed_s, no_property, property<edge_weight_t, float>>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::mt19937 gen(42);

  for (int i = 0; i < 5; ++i) {
    Graph g;
    generate_random_graph(g, 10, 50, gen);

    auto idx = get(vertex_index, g);
    auto weight = get(edge_weight, g);
    for (auto e : edges(g)) {
      std::uniform_real_distribution<float> uni_dist;
      put(weight, e, uni_dist(gen) + 1e-5F);
    }

    auto pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
    auto dm = vector_property_map(num_vertices(g), idx, default_inf_v<float>);
    auto ncpm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
    auto ncdm = vector_property_map(num_vertices(g), idx, default_inf_v<float>);

    // pick random start/goal
    Vertex start = gen() % num_vertices(g);

    auto print_problem_def = [&]() {
      std::stringstream out;
      print_graph(g, get(vertex_index, g), weight, out);
      out << "\nStart vertex: " << start;
      return out.str();
    };
    // Print problem definition if any checks fail.
    SCOPED_TRACE(print_problem_def());

    // call astar named parameter interface

    dijkstra_shortest_paths_no_color_map(g, std::ranges::single_view{start}, ncpm.ref(), ncdm.ref(), weight, idx,
                                         default_dijkstra_visitor());
    dijkstra_shortest_paths(g, start, pm.ref(), dm.ref(), weight, idx, default_dijkstra_visitor());

    for (auto u : vertices(g)) {
      EXPECT_EQ(pm[u], ncpm[u]);
      EXPECT_THAT(dm[u], ::testing::FloatNear(ncdm[u], 1.0e-5F));
    }
  }
}

}  // namespace
}  // namespace bagl
