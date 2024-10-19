// Copyright (c) 2004 Kristopher Beevers
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/astar_search.h"

#include <cmath>
#include <iostream>
#include <list>
#include <random>
#include <tuple>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test/new_york_state_map.h"

namespace bagl {
namespace {

TEST(AStarSearchTest, Search) {
  // specify some types
  using Graph = test::state_map_type;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  std::mt19937 gen(42);

  Graph g = test::make_new_york_state_map();

  auto name = get(vertex_name, g);
  auto idx = get(vertex_index, g);
  auto weight = get(edge_weight, g);

  auto pm = vector_property_map(num_vertices(g), idx, graph_traits<Graph>::null_vertex());
  auto cm = vector_property_map(num_vertices(g), idx, std::numeric_limits<float>::max());
  auto dm = vector_property_map(num_vertices(g), idx, std::numeric_limits<float>::max());
  auto colorm = vector_property_map(num_vertices(g), idx, default_color_type::white_color);

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
    EXPECT_THROW(astar_search(g, start, test::travel_time_heuristic(g, goal),
                              make_astar_visitor(bfs_stop_at_goal(goal)), pm, cm, dm, weight, idx, colorm),
                 search_succeeded)
        << "Did not find a path to the goal at all!";
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
    astar_search(g, start, test::travel_time_heuristic(g, goal), default_astar_visitor(), pm, cm, dm, weight, idx,
                 colorm);
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

}  // namespace
}  // namespace bagl
