// Copyright 2004 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/gursoy_atun_layout.h"

#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(GursoyAtunLayoutTest, GridInHeart) {
  // Generate a graph structured like a grid
  // lay it out in a heart topology

  using Graph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_weight_t, double> >;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  // Make grid, like Gursoy and Atun used
  std::vector<std::vector<Vertex> > verts;
  const int grid_size = 20;
  std::mt19937 gen(42);
  std::uniform_real_distribution<> random_edge_weight;
  verts.resize(grid_size);
  for (auto& v_row : verts) {
    v_row.resize(grid_size);
    for (auto& v : v_row) {
      v = add_vertex(g);
    }
  }
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      if (i != 0) {
        add_edge(verts[i][j], verts[i - 1][j], g, random_edge_weight(gen));
      }
      if (j != 0) {
        add_edge(verts[i][j], verts[i][j - 1], g, random_edge_weight(gen));
      }
    }
  }

  auto vindex = get(vertex_index, g);

  using Topology = heart_topology<>;
  using Point = Topology::point_type;
  Topology space;

  auto position = vector_property_map(num_vertices(g), vindex, Point{});

  gursoy_atun_layout(g, space, position.ref());
  for (auto v : vertices(g)) {
    EXPECT_TRUE(space.in_heart(get(position, v)));
  }

  gursoy_atun_layout(g, space, position.ref(), get(edge_weight, g));
  for (auto v : vertices(g)) {
    EXPECT_TRUE(space.in_heart(get(position, v)));
  }
}

TEST(GursoyAtunLayoutTest, UnweightedGridInHeart) {
  // Generate a graph structured like a grid
  // lay it out in a heart topology

  using Graph = adjacency_list<vec_s, vec_s, undirected_s>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Graph g;

  // Make grid, like Gursoy and Atun used
  std::vector<std::vector<Vertex> > verts;
  const int grid_size = 40;
  std::mt19937 gen(42);
  std::uniform_real_distribution<> random_edge_weight;
  verts.resize(grid_size);
  for (auto& v_row : verts) {
    v_row.resize(grid_size);
    for (auto& v : v_row) {
      v = add_vertex(g);
    }
  }
  for (int i = 0; i < grid_size; ++i) {
    for (int j = 0; j < grid_size; ++j) {
      if (i != 0) {
        add_edge(verts[i][j], verts[i - 1][j], g);
      }
      if (j != 0) {
        add_edge(verts[i][j], verts[i][j - 1], g);
      }
    }
  }

  auto vindex = get(vertex_index, g);

  using Topology = heart_topology<>;
  using Point = Topology::point_type;
  Topology space;

  auto position = vector_property_map(num_vertices(g), vindex, Point{});

  gursoy_atun_layout(g, space, position.ref());
  for (auto v : vertices(g)) {
    auto pos = get(position, v);
    EXPECT_TRUE(space.in_heart(pos));
  }
  // Heart topology is about 2000 x 2000 in size, so expect edges to be reasonably
  // distributed.
  int below_20_count = 0;
  int above_100_count = 0;
  int above_200_count = 0;
  for (auto e : edges(g)) {
    auto d = space.distance(get(position, source(e, g)), get(position, target(e, g)));
    if (d < 20.0) {
      ++below_20_count;
    }
    if (d > 100.0) {
      ++above_100_count;
    }
    if (d > 200.0) {
      ++above_200_count;
    }
  }
  EXPECT_LT(above_200_count, num_edges(g) / 100);
  EXPECT_LT(above_100_count, num_edges(g) / 10);
  EXPECT_LT(below_20_count, num_edges(g) / 10);
}

}  // namespace
}  // namespace bagl
