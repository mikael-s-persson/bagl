// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/biconnected_components.h"

#include <algorithm>
#include <iterator>
#include <random>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/connected_components.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct EdgeProperty {
  std::size_t component;
};

using Graph = adjacency_list<list_s, vec_s, undirected_s, no_property, EdgeProperty>;
using Vertex = graph_vertex_descriptor_t<Graph>;

void check_articulation_points(const Graph& g, std::vector<Vertex> art_points) {
  auto components = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  int basic_comps = connected_components(g, components.ref());

  std::vector<Vertex> art_points_check;
  for (auto v : vertices(g)) {
    Graph g_copy(g);
    Vertex victim = *std::next(vertices(g_copy).begin(), get(vertex_index, g, v));
    clear_vertex(victim, g_copy);
    remove_vertex(victim, g_copy);

    int copy_comps = connected_components(g_copy, components.ref());
    if (copy_comps > basic_comps) {
      art_points_check.push_back(v);
    }
  }

  EXPECT_THAT(art_points, ::testing::UnorderedElementsAreArray(art_points_check)) << "Failed on graph:\n"
                                                                                  << (print_graph(g), "");
}

TEST(BiconnectedCompTest, Basic) {
  Graph g(4);
  add_edge(2, 3, g);
  add_edge(0, 3, g);
  add_edge(0, 2, g);
  add_edge(1, 0, g);

  std::vector<Vertex> art_points;
  biconnected_components(g, get(&EdgeProperty::component, g), std::back_inserter(art_points));

  check_articulation_points(g, art_points);
}

TEST(BiconnectedCompTest, RandomGraph) {
  std::size_t n = 100;
  std::size_t m = 500;
  std::size_t seed = 5453;

  Graph g(n);
  std::mt19937 gen(seed);
  generate_random_graph(g, n, m, gen);

  std::vector<Vertex> art_points;
  biconnected_components(g, get(&EdgeProperty::component, g), std::back_inserter(art_points));

  check_articulation_points(g, art_points);
}

}  // namespace
}  // namespace bagl
