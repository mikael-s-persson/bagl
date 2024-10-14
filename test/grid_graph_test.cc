// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/grid_graph.h"

#include <ctime>
#include <iostream>
#include <random>
#include <set>
#include <type_traits>

#include "bagl/graph_traits.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// Function that prints a vertex to std::cout
template <typename Vertex>
void print_vertex(Vertex vertex_to_print) {
  std::cout << "(";
  for (std::size_t i = 0; i < vertex_to_print.size(); ++i) {
    std::cout << vertex_to_print[i];
    if (i + 1 != vertex_to_print.size()) {
      std::cout << ", ";
    }
  }
  std::cout << ")";
}

template <typename T>
class GridGraphTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(GridGraphTest);

TYPED_TEST_P(GridGraphTest, Basic) {
  constexpr int Dims = TypeParam::value;
  std::mt19937 rng_gen(42);
  using Graph = grid_graph<Dims>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;

  std::cout << "Dimensions: " << Dims << ", lengths: ";

  // Randomly generate the dimension lengths (3-10) and wrapping
  std::array<std::size_t, Dims> lengths;
  std::array<bool, Dims> wrapped;

  for (std::size_t i = 0; i < Dims; ++i) {
    lengths[i] = 3 + (rng_gen() % 8);
    wrapped[i] = ((rng_gen() % 2) == 0);

    std::cout << lengths[i] << (wrapped[i] ? " [W]" : " [U]") << ", ";
  }

  std::cout << std::endl;

  Graph graph(lengths, wrapped);

  // Verify dimension lengths and wrapping
  for (std::size_t i = 0; i < Dims; ++i) {
    EXPECT_EQ(graph.length(i), lengths[i]);
    EXPECT_EQ(graph.wrapped(i), wrapped[i]);
  }

  // Verify matching indices
  for (std::size_t vertex_i = 0; vertex_i < num_vertices(graph); ++vertex_i) {
    EXPECT_EQ(get(vertex_index, graph, vertex(vertex_i, graph)), vertex_i);
  }

  for (std::size_t edge_i = 0; edge_i < num_edges(graph); ++edge_i) {
    Edge current_edge = edge_at(edge_i, graph);
    EXPECT_EQ(get(edge_index, graph, current_edge), edge_i);
  }

  // Verify all vertices are within bounds
  std::size_t vertex_count = 0;
  for (Vertex current_vertex : vertices(graph)) {
    std::size_t current_index = get(vertex_index, graph, current_vertex);

    for (std::size_t i = 0; i < Dims; ++i) {
      EXPECT_LT(current_vertex[i], lengths[i]);
    }

    // Verify out-edges of this vertex
    std::size_t out_edge_count = 0;
    std::set<std::size_t> target_vertices;

    for (Edge out_edge : out_edges(current_vertex, graph)) {
      target_vertices.insert(get(vertex_index, graph, target(out_edge, graph)));
      ++out_edge_count;
    }

    EXPECT_EQ(out_edge_count, out_degree(current_vertex, graph));

    // Verify in-edges of this vertex
    std::size_t in_edge_count = 0;
    for (Edge in_edge : in_edges(current_vertex, graph)) {
      EXPECT_GT(target_vertices.count(get(vertex_index, graph, source(in_edge, graph))), 0);
      ++in_edge_count;
    }

    EXPECT_EQ(in_edge_count, in_degree(current_vertex, graph));

    // The number of out-edges and in-edges should be the same
    EXPECT_EQ(degree(current_vertex, graph), out_degree(current_vertex, graph) + in_degree(current_vertex, graph));

    // Verify adjacent vertices to this vertex
    std::size_t adjacent_count = 0;

    for (Vertex adjacent_vertex : adjacent_vertices(current_vertex, graph)) {
      EXPECT_GT(target_vertices.count(get(vertex_index, graph, adjacent_vertex)), 0);
      ++adjacent_count;
    }

    EXPECT_EQ(adjacent_count, out_degree(current_vertex, graph));

    // Verify that this vertex is not listed as connected to any
    // vertices outside of its adjacent vertices.
    for (Vertex unconnected_vertex : vertices(graph)) {
      std::size_t unconnected_index = get(vertex_index, graph, unconnected_vertex);

      if ((unconnected_index == current_index) || (target_vertices.count(unconnected_index) > 0)) {
        continue;
      }

      EXPECT_FALSE(edge(current_vertex, unconnected_vertex, graph).second);
      EXPECT_FALSE(edge(unconnected_vertex, current_vertex, graph).second);
    }

    ++vertex_count;
  }

  EXPECT_EQ(vertex_count, num_vertices(graph));

  // Verify all edges are within bounds
  std::size_t edge_count = 0;
  for (Edge current_edge : edges(graph)) {
    std::size_t source_index = get(vertex_index, graph, source(current_edge, graph));
    std::size_t target_index = get(vertex_index, graph, target(current_edge, graph));

    EXPECT_NE(source_index, target_index);
    EXPECT_LT(source_index, num_vertices(graph));
    EXPECT_LT(target_index, num_vertices(graph));

    // Verify that the edge is listed as existing in both directions
    EXPECT_TRUE(edge(source(current_edge, graph), target(current_edge, graph), graph).second);
    EXPECT_TRUE(edge(target(current_edge, graph), source(current_edge, graph), graph).second);

    ++edge_count;
  }

  EXPECT_EQ(edge_count, num_edges(graph));
}

REGISTER_TYPED_TEST_SUITE_P(GridGraphTest, Basic);

using GridGraphDimensions = ::testing::Types<std::integral_constant<int, 1>, std::integral_constant<int, 2>,
                                             std::integral_constant<int, 3>, std::integral_constant<int, 4>>;

INSTANTIATE_TYPED_TEST_SUITE_P(GridGraph, GridGraphTest, GridGraphDimensions);

}  // namespace
}  // namespace bagl
