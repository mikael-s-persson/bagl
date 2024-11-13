// Copyright (C) 2004-2008 The Trustees of Indiana University.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/metis.h"

#include <array>
#include <fstream>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/breadth_first_search.h"
#include "bagl/graph_utility.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

/* An undirected graph with distance values stored on the vertices. */
using Graph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_weight_t, std::size_t> >;

TEST(MetisTest, WeightedGraph) {
  const char* filename = "test/testdata/weighted_graph.gr";

  // Open the METIS input file
  std::ifstream in(filename);
  ASSERT_TRUE(in.good());
  metis_reader reader(in);

  // Load the graph using the default distribution
  Graph g(reader.num_vertices(), reader.edges(), reader.weights());
  EXPECT_EQ(num_vertices(g), 5);
  EXPECT_EQ(num_edges(g), 9);
  using ETrip = std::array<std::size_t, 3>;
  std::vector<ETrip> actual_edges;
  for (auto e : edges(g)) {
    actual_edges.push_back(
        ETrip{{std::min(source(e, g), target(e, g)), std::max(source(e, g), target(e, g)), get(edge_weight, g, e)}});
  }
  EXPECT_THAT(actual_edges, ::testing::UnorderedElementsAreArray({
                                ETrip{{0, 2, 1}},
                                ETrip{{1, 1, 2}},
                                ETrip{{1, 3, 1}},
                                ETrip{{1, 4, 2}},
                                ETrip{{1, 2, 7}},
                                ETrip{{2, 3, 3}},
                                ETrip{{3, 4, 1}},
                                ETrip{{0, 4, 1}},
                                ETrip{{1, 4, 1}},
                            }));
  if (::testing::Test::HasFailure()) {
    print_graph(g);
    for (auto e : edges(g)) {
      std::cout << source(e, g) << "-" << target(e, g) << ": " << get(edge_weight, g, e) << std::endl;
    }
  }
}

}  // namespace
}  // namespace bagl
