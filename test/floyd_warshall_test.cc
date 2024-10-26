// Copyright 2002 Rensselaer Polytechnic Institute
//  Authors: Lauren Foutz
//           Scott Hill
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <algorithm>
#include <iostream>
#include <map>
#include <random>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/bellman_ford_shortest_paths.h"
#include "bagl/floyd_warshall_shortest.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/matrix_property_map.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(FloydWarshallTest, UndirectedAdjMatrix) {
  using Graph = adjacency_matrix<undirected_s, property<vertex_distance_t, int, property<vertex_name_t, int> >,
                                 property<edge_weight_t, int> >;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  int num_v = 20;
  Graph g(num_v);
  int num_e = 400;
  std::mt19937 ran(46);

  auto index = get(vertex_name, g);
  int x = 0;
  for (auto v : vertices(g)) {
    put(index, v, x++);
  }

  for (int i = 0; i < num_e; i++) {
    add_edge(index[ran() % num_v], index[ran() % num_v], g);
  }

  auto local_edge_map = get(edge_weight, g);
  for (auto e : edges(g)) {
    if (ran() % 100 != 0) {
      put(local_edge_map, e, ran() % 100);
    } else {
      put(local_edge_map, e, 0 - (ran() % 100));
    }
  }

  int int_inf = std::numeric_limits<int>::max();
  auto matrix = matrix_property_store(num_vertices(g), index, int_inf);
  bool floyd = floyd_warshall_all_pairs_shortest_paths(g, matrix.ref(), local_edge_map);

  bool bellman = true;
  auto distance_row = get(vertex_distance, g);
  for (auto v : vertices(g)) {
    put(distance_row, v, int_inf);
  }
  auto matrix2 = matrix_property_store(num_vertices(g), index, int_inf);
  for (auto v : vertices(g)) {
    put(distance_row, v, 0);
    bellman = bellman_ford_shortest_paths(g, default_bellman_visitor(), num_v, get(edge_weight, g),
                                          null_property_map<Vertex>(), get(vertex_distance, g));
    for (auto u : vertices(g)) {
      matrix2[v][u] = get(distance_row, u);
      put(distance_row, u, int_inf);
    }
    if (!bellman) {
      break;
    }
  }

  EXPECT_EQ(floyd, bellman);
  ASSERT_TRUE(bellman) << "Must pick a suitable test seed";

  if (!bellman) {
    return;
  }

  for (auto u : vertices(g)) {
    for (auto v : vertices(g)) {
      EXPECT_EQ(matrix2[u][v], matrix[u][v]) << "At u = " << index[u] << " v = " << index[v];
    }
  }
}

TEST(FloydWarshallTest, DirectedAdjList) {
  using Graph =
      adjacency_list<list_s, list_s, directed_s, property<vertex_distance_t, int, property<vertex_name_t, int> >,
                     property<edge_weight_t, int> >;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  int num_v = 20;
  int num_e = 400;

  std::mt19937 ran(44);
  Graph g;
  generate_random_graph(g, num_v, num_e, ran, true);

  auto index = get(vertex_name, g);
  int x = 0;
  for (auto v : vertices(g)) {
    put(index, v, x++);
  }

  auto local_edge_map = get(edge_weight, g);
  for (auto e : edges(g)) {
    if (ran() % 100 != 0) {
      put(local_edge_map, e, ran() % 100);
    } else {
      put(local_edge_map, e, 0 - (ran() % 100));
    }
  }

  int int_inf = std::numeric_limits<int>::max();
  auto matrix = matrix_property_store(num_vertices(g), index, int_inf);
  bool floyd = floyd_warshall_all_pairs_shortest_paths(g, matrix.ref(), local_edge_map);

  bool bellman = true;
  auto distance_row = get(vertex_distance, g);
  for (auto v : vertices(g)) {
    put(distance_row, v, int_inf);
  }
  auto matrix2 = matrix_property_store(num_vertices(g), index, int_inf);
  for (auto v : vertices(g)) {
    put(distance_row, v, 0);
    bellman = bellman_ford_shortest_paths(g, default_bellman_visitor(), num_v, local_edge_map,
                                          null_property_map<Vertex>(), distance_row);
    for (auto u : vertices(g)) {
      matrix2[v][u] = get(distance_row, u);
      put(distance_row, u, int_inf);
    }
    if (!bellman) {
      break;
    }
  }

  EXPECT_EQ(floyd, bellman);
  ASSERT_TRUE(bellman) << "Must pick a suitable test seed";

  for (auto u : vertices(g)) {
    for (auto v : vertices(g)) {
      EXPECT_EQ(matrix2[u][v], matrix[u][v]) << "At u = " << index[u] << " v = " << index[v];
    }
  }
}

}  // namespace
}  // namespace bagl
