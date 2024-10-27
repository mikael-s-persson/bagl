// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/incremental_components.h"

#include <ctime>
#include <iostream>
#include <map>
#include <random>
#include <set>

#include "bagl/adjacency_list.h"
#include "bagl/disjoint_sets.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "bagl/reverse_graph.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using VectorGraph = adjacency_list<vec_s, vec_s, undirected_s>;
using ListGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, unsigned int> >;

template <typename Graph>
void test_graph(const Graph& g) {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  auto vindex = get(vertex_index, g);
  auto parent_map = vector_property_map(num_vertices(g), vindex, graph_traits<Graph>::null_vertex());
  auto rank_map = vector_property_map(num_vertices(g), vindex, std::size_t{0});

  // Create disjoint sets of vertices from the g
  disjoint_sets vertex_sets(rank_map.ref(), parent_map.ref());

  initialize_incremental_components(g, vertex_sets);
  incremental_components(g, vertex_sets);

  // Build component index from the g's vertices, its index map,
  // and the disjoint sets.
  component_index vertex_components(parent_map.ref(), vindex, vertices(g));

  // Create a reverse-lookup map for vertex indices
  std::vector<Vertex> rindex(num_vertices(g));
  for (Vertex v : vertices(g)) {
    rindex[get(vindex, v)] = v;
  }

  const auto by_index = [&](Vertex u, Vertex v) { return get(vindex, u) < get(vindex, v); };

  // Verify that components are really connected
  for (std::size_t cindex : vertex_components.components()) {
    std::set<Vertex, decltype(by_index)> component_vertices(by_index);

    for (std::size_t child_index : vertex_components.component(cindex)) {
      Vertex child_vertex = rindex[child_index];
      component_vertices.insert(child_vertex);
    }

    // Verify that children are connected to each other in some
    // manner, but not to vertices outside their component.
    for (Vertex child_vertex : component_vertices) {
      // Skip orphan vertices
      if (out_degree(child_vertex, g) == 0) {
        continue;
      }

      // Make sure at least one edge exists between [child_vertex] and
      // another vertex in the component.
      bool edge_exists = false;

      for (Edge child_edge : out_edges(child_vertex, g)) {
        if (component_vertices.count(target(child_edge, g)) > 0) {
          edge_exists = true;
          break;
        }
      }

      EXPECT_TRUE(edge_exists);
    }
  }
}

TEST(IncrementalComponentsTest, Basic) {
  std::size_t vertices_to_generate = 100;
  std::size_t edges_to_generate = 50;
  std::size_t random_seed = 42;

  std::mt19937 generator(random_seed);

  // Generate random vector and list graphs
  VectorGraph vector_graph;
  generate_random_graph(vector_graph, vertices_to_generate, edges_to_generate, generator, false);

  test_graph(vector_graph);

  ListGraph list_graph;
  generate_random_graph(list_graph, vertices_to_generate, edges_to_generate, generator, false);

  // Assign indices to list_graph's vertices
  std::size_t index = 0;
  for (auto v : vertices(list_graph)) {
    put(get(vertex_index, list_graph), v, index++);
  }

  test_graph(list_graph);
}

}  // namespace
}  // namespace bagl
