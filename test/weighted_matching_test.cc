// Copyright (c) 2018 Yi Ji
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/graph_traits.h"
#include "bagl/max_cardinality_matching.h"
#include "bagl/maximum_weighted_matching.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using EdgeProperty = property<edge_weight_t, float, property<edge_index_t, std::size_t> >;

using UndirectedGraph = adjacency_list<vec_s, vec_s, undirected_s, no_property, EdgeProperty>;
using UndirectedListGraph =
    adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, std::size_t>, EdgeProperty>;
using UndirectedMatrixGraph = adjacency_matrix<undirected_s, no_property, EdgeProperty>;

template <typename Graph>
struct vertex_index_installer {
  static void install(Graph&) {}
};

template <>
struct vertex_index_installer<UndirectedListGraph> {
  static void install(UndirectedListGraph& g) {
    std::size_t i = 0;
    for (auto v : vertices(g)) {
      put(vertex_index, g, v, i++);
    }
  }
};

template <typename Graph>
void print_graph(const Graph& g) {
  std::cout << std::endl << "The graph is:\n";
  for (auto e : edges(g)) {
    std::cout << "(" << get(vertex_index, g, source(e, g)) << ", " << get(vertex_index, g, target(e, g))
              << "): " << get(edge_weight, g, e) << std::endl;
  }
}

template <typename Graph>
void weighted_matching_test(const Graph& g, property_traits_value_t<property_map_t<Graph, edge_weight_t> > answer) {
  auto mate = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
  maximum_weighted_matching(g, mate.ref());
  EXPECT_EQ(matching_weight_sum(g, mate.ref()), answer);
  if (matching_weight_sum(g, mate.ref()) == answer) {
    return;
  }

  // Print problem information in case of a failure.

  auto max_mate = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
  brute_force_maximum_weighted_matching(g, max_mate.ref());

  std::cout << std::endl
            << "Found a weighted matching of weight sum " << matching_weight_sum(g, mate.ref()) << std::endl
            << "While brute-force search found a weighted matching of weight sum "
            << matching_weight_sum(g, max_mate.ref()) << std::endl;

  print_graph(g);

  std::cout << std::endl << "The algorithmic matching is:" << std::endl;
  for (auto v : vertices(g)) {
    if (mate[v] != graph_traits<Graph>::null_vertex()) {
      auto v_id = get(vertex_index, g, v);
      auto mate_id = get(vertex_index, g, mate[v]);
      if (v_id < mate_id) {
        std::cout << "{" << v_id << ", " << mate_id << "}" << std::endl;
      }
    }
  }

  std::cout << std::endl << "The brute-force matching is:" << std::endl;
  for (auto v : vertices(g)) {
    if (max_mate[v] != graph_traits<Graph>::null_vertex()) {
      auto v_id = get(vertex_index, g, v);
      auto mate_id = get(vertex_index, g, max_mate[v]);
      if (v_id < mate_id) {
        std::cout << "{" << v_id << ", " << mate_id << "}" << std::endl;
      }
    }
  }

  std::cout << std::endl;
}

template <typename Graph>
Graph make_graph(std::size_t num_v, std::size_t num_e, std::deque<std::size_t> input_edges) {
  Graph g(num_v);
  vertex_index_installer<Graph>::install(g);
  for (std::size_t i = 0; i < num_e; ++i) {
    auto src_v = input_edges.front();
    input_edges.pop_front();
    auto tgt_v = input_edges.front();
    input_edges.pop_front();
    auto edge_weight = input_edges.front();
    input_edges.pop_front();
    add_edge(*std::next(vertices(g).begin(), src_v), *std::next(vertices(g).begin(), tgt_v), g,
             EdgeProperty(edge_weight));
  }
  return g;
}

TEST(WeightedMatchingTest, Basic) {
  // Change this to "weighted_matching.dat" for exhaustive test suite (50k).
  std::ifstream in_file("test/testdata/weighted_matching_small.dat");
  ASSERT_TRUE(in_file);
  std::string line;
  while (std::getline(in_file, line)) {
    std::istringstream in_graph(line);
    std::size_t answer, num_v, num_e;
    in_graph >> answer >> num_v >> num_e;

    std::deque<std::size_t> input_edges;
    std::size_t i;
    while (in_graph >> i) {
      input_edges.push_back(i);
    }

    weighted_matching_test(make_graph<UndirectedGraph>(num_v, num_e, input_edges), answer);
    weighted_matching_test(make_graph<UndirectedListGraph>(num_v, num_e, input_edges), answer);
    weighted_matching_test(make_graph<UndirectedMatrixGraph>(num_v, num_e, input_edges), answer);
  }
}

}  // namespace
}  // namespace bagl
