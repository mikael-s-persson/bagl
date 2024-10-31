// Copyright Daniel Trebbien 2010.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/connected_components.h"
#include "bagl/exception.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/read_dimacs.h"
#include "bagl/single_property_map.h"
#include "bagl/stoer_wagner_min_cut.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

typedef adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_weight_t, int> > undirected_graph;
typedef property_map_t<undirected_graph, edge_weight_t> weight_map_type;
typedef property_traits_value_t<weight_map_type> weight_type;

typedef adjacency_list<vec_s, vec_s, undirected_s> undirected_unweighted_graph;

std::string test_dir;

using edge_t = std::pair<std::size_t, std::size_t>;

TEST(StoerWagnerTest, OneVertexGraph) {
  // if only one vertex, can't run `stoer_wagner_min_cut`
  undirected_graph g;
  add_vertex(g);

  EXPECT_THROW(stoer_wagner_min_cut(g, get(edge_weight, g)), bad_graph);
}

// the example from Stoer & Wagner (1997)
TEST(StoerWagnerTest, Test0) {
  std::vector<edge_t> edges = {{0, 1}, {1, 2}, {2, 3}, {0, 4}, {1, 4}, {1, 5},
                               {2, 6}, {3, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}};
  std::vector<int> ws = {2, 3, 4, 3, 2, 2, 2, 2, 2, 3, 1, 3};
  undirected_graph g(8, edges, ws);

  auto weights = get(edge_weight, g);
  auto parities = vector_property_map(num_vertices(g), get(vertex_index, g), false);
  int w = stoer_wagner_min_cut(g, weights, parities.ref(), get(vertex_index, g));
  EXPECT_EQ(w, 4);
  const bool parity0 = get(parities, 0);
  EXPECT_EQ(parity0, get(parities, 1));
  EXPECT_EQ(parity0, get(parities, 4));
  EXPECT_EQ(parity0, get(parities, 5));
  const bool parity2 = get(parities, 2);
  EXPECT_NE(parity0, parity2);
  EXPECT_EQ(parity2, get(parities, 3));
  EXPECT_EQ(parity2, get(parities, 6));
  EXPECT_EQ(parity2, get(parities, 7));
}

// three vertices with one multi-edge
TEST(StoerWagnerTest, Test1) {
  std::vector<edge_t> edges = {{0, 1}, {1, 2}, {1, 2}, {2, 0}};
  std::vector<int> ws = {3, 1, 1, 1};
  undirected_graph g(3, edges, ws);

  auto weights = get(edge_weight, g);
  auto parities = vector_property_map(num_vertices(g), get(vertex_index, g), false);
  int w = stoer_wagner_min_cut(g, weights, parities.ref(), get(vertex_index, g));
  EXPECT_EQ(w, 3);
  const bool parity2 = get(parities, 2), parity0 = get(parities, 0);
  EXPECT_NE(parity2, parity0);
  EXPECT_EQ(parity0, get(parities, 1));
}

// example by Daniel Trebbien
TEST(StoerWagnerTest, Test2) {
  std::vector<edge_t> edges = {{5, 2}, {0, 6}, {5, 6}, {3, 1}, {0, 1}, {6, 3}, {4, 6}, {2, 4}, {5, 3}};
  std::vector<int> ws = {1, 3, 4, 6, 4, 1, 2, 5, 2};
  undirected_graph g(7, edges, ws);

  auto weights = get(edge_weight, g);
  auto parities = vector_property_map(num_vertices(g), get(vertex_index, g), false);
  int w = stoer_wagner_min_cut(g, weights, parities.ref(), get(vertex_index, g));
  EXPECT_EQ(w, 3);
  const bool parity2 = get(parities, 2);
  EXPECT_EQ(parity2, get(parities, 4));
  const bool parity5 = get(parities, 5);
  EXPECT_NE(parity2, parity5);
  EXPECT_EQ(parity5, get(parities, 3));
  EXPECT_EQ(parity5, get(parities, 6));
  EXPECT_EQ(parity5, get(parities, 1));
  EXPECT_EQ(parity5, get(parities, 0));
}

// example by Daniel Trebbien
TEST(StoerWagnerTest, Test3) {
  std::vector<edge_t> edges = {{3, 4}, {3, 6}, {3, 5}, {0, 4}, {0, 1}, {0, 6}, {0, 7}, {0, 5},
                               {0, 2}, {4, 1}, {1, 6}, {1, 5}, {6, 7}, {7, 5}, {5, 2}, {3, 4}};
  std::vector<int> ws = {0, 3, 1, 3, 1, 2, 6, 1, 8, 1, 1, 80, 2, 1, 1, 4};
  undirected_graph g(8, edges, ws);

  auto weights = get(edge_weight, g);
  auto parities = vector_property_map(num_vertices(g), get(vertex_index, g), false);
  int w = stoer_wagner_min_cut(g, weights, parities.ref(), get(vertex_index, g));
  EXPECT_EQ(w, 7);
  const bool parity1 = get(parities, 1);
  EXPECT_EQ(parity1, get(parities, 5));
  const bool parity0 = get(parities, 0);
  EXPECT_NE(parity1, parity0);
  EXPECT_EQ(parity0, get(parities, 2));
  EXPECT_EQ(parity0, get(parities, 3));
  EXPECT_EQ(parity0, get(parities, 4));
  EXPECT_EQ(parity0, get(parities, 6));
  EXPECT_EQ(parity0, get(parities, 7));
}

TEST(StoerWagnerTest, Test4) {
  std::vector<edge_t> edges = {{0, 1}, {1, 2}, {2, 3}, {0, 4}, {1, 4}, {1, 5}, {2, 6},
                               {3, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}, {0, 4}, {6, 7}};
  undirected_unweighted_graph g(8, edges);

  auto parities = vector_property_map(num_vertices(g), get(vertex_index, g), false);
  int w = stoer_wagner_min_cut(g, single_property_map(1), parities.ref(), get(vertex_index, g));
  EXPECT_EQ(w, 2);
  const bool parity0 = get(parities, 0);
  EXPECT_EQ(parity0, get(parities, 1));
  EXPECT_EQ(parity0, get(parities, 4));
  EXPECT_EQ(parity0, get(parities, 5));
  const bool parity2 = get(parities, 2);
  EXPECT_NE(parity0, parity2);
  EXPECT_EQ(parity2, get(parities, 3));
  EXPECT_EQ(parity2, get(parities, 6));
  EXPECT_EQ(parity2, get(parities, 7));
}

// Non regression test for github.com/boostorg/graph/issues/286
TEST(StoerWagnerTest, Test5) {
  std::vector<edge_t> edges = {{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}, {4, 5},
                               {4, 6}, {4, 7}, {5, 6}, {5, 7}, {6, 7}, {0, 4}};
  std::vector<int> ws = {3, 3, 3, 2, 2, 2, 3, 3, 3, 2, 2, 2, 6};
  undirected_graph g(8, edges, ws);

  auto weights = get(edge_weight, g);
  auto parities = vector_property_map(num_vertices(g), get(vertex_index, g), false);
  int w = stoer_wagner_min_cut(g, weights, parities.ref(), get(vertex_index, g));
  EXPECT_EQ(w, 6);
  const bool parity0 = get(parities, 0);
  EXPECT_EQ(parity0, get(parities, 1));
  EXPECT_EQ(parity0, get(parities, 2));
  EXPECT_EQ(parity0, get(parities, 3));
  const bool parity4 = get(parities, 4);
  EXPECT_NE(parity0, parity4);
  EXPECT_EQ(parity4, get(parities, 5));
  EXPECT_EQ(parity4, get(parities, 6));
  EXPECT_EQ(parity4, get(parities, 7));
}

// The input for the `test_prgen` family of tests comes from a program, named
// `prgen`, that comes with a package of min-cut solvers by Chandra Chekuri,
// Andrew Goldberg, David Karger, Matthew Levine, and Cliff Stein. `prgen` was
// used to generate input graphs and the solvers were used to verify the return
// value of `stoer_wagner_min_cut` on the input graphs.
//
// http://www.columbia.edu/~cs2035/code.html
//
// Note that it is somewhat more difficult to verify the parities because
// "`prgen` graphs" often have several min-cuts. This is why only the cut
// weight of the min-cut is verified.

// 3 min-cuts
TEST(StoerWagnerTest, PRGen_20_70_2) {
  std::ifstream ifs("test/testdata/prgen_20_70_2.net");
  undirected_graph g;
  read_dimacs_min_cut(g, get(edge_weight, g), dummy_property_map(), ifs);

  // verify the connectedness assumption
  auto components = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  EXPECT_EQ(connected_components(g, components.ref()), 1U);

  auto distances = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  auto indices_in_heap =
      vector_property_map(num_vertices(g), get(vertex_index, g), std::numeric_limits<std::size_t>::max());
  auto pq = make_d_ary_heap_indirect<graph_vertex_descriptor_t<undirected_graph>, 22>(
      distances.ref(), indices_in_heap.ref(), std::greater<>());

  int w = stoer_wagner_min_cut(g, get(edge_weight, g), pq);
  EXPECT_EQ(w, 3407);
}

// 7 min-cuts
TEST(StoerWagnerTest, PRGen_50_40_2) {
  std::ifstream ifs("test/testdata/prgen_50_40_2.net");
  undirected_graph g;
  read_dimacs_min_cut(g, get(edge_weight, g), dummy_property_map(), ifs);

  // verify the connectedness assumption
  auto components = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  EXPECT_EQ(connected_components(g, components.ref()), 1U);

  int w = stoer_wagner_min_cut(g, get(edge_weight, g));
  EXPECT_EQ(w, 10056);
}

// 6 min-cuts
TEST(StoerWagnerTest, PRGen_50_70_2) {
  std::ifstream ifs("test/testdata/prgen_50_70_2.net");
  undirected_graph g;
  read_dimacs_min_cut(g, get(edge_weight, g), dummy_property_map(), ifs);

  // verify the connectedness assumption
  auto components = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  EXPECT_EQ(connected_components(g, components.ref()), 1U);

  int w = stoer_wagner_min_cut(g, get(edge_weight, g));
  EXPECT_EQ(w, 21755);
}

}  // namespace
}  // namespace bagl
