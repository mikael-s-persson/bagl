// Copyright 2024 Mikael Persson

#include <cstddef>
#include <fstream>
#include <sstream>

#include "bagl/adjacency_list.h"
#include "bagl/connected_components.h"
#include "bagl/directed_graph.h"
#include "bagl/read_dimacs.h"
#include "bagl/vector_property_map.h"
#include "bagl/write_dimacs.h"
#include "gtest/gtest.h"
#include "test/min_cost_max_flow_utils.h"

namespace bagl {
namespace {

// Min-cut problem
TEST(ReadDimacsTest, MinCutProblem) {
  using Graph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_weight_t, int> >;

  std::ifstream ifs("test/testdata/prgen_20_70_2.net");
  Graph g;
  EXPECT_TRUE(read_dimacs_min_cut(g, get(edge_weight, g), dummy_property_map(), ifs));

  // verify the connectedness assumption
  auto components = vector_property_map(num_vertices(g), get(vertex_index, g), std::size_t{0});
  EXPECT_EQ(connected_components(g, components.ref()), 1U);

  EXPECT_EQ(num_vertices(g), 20);
  EXPECT_EQ(num_edges(g), 270);
}

// Min-cut problem
TEST(ReadWriteDimacsTest, MaxFlowProblem) {
  using Graph = sample_graph::Graph;
  using Edge = graph_edge_descriptor_t<Graph>;
  Graph g;
  auto rev = vector_property_map(0, get(edge_index, g), Edge{});
  auto [s, t] = sample_graph::get_sample_graph_2(g, rev.ref());

  std::stringstream dimacs_stream;
  write_dimacs_max_flow(g, get(edge_capacity, g), get(vertex_index, g), s, t, dimacs_stream);

  auto s2 = graph_traits<Graph>::null_vertex();
  auto t2 = graph_traits<Graph>::null_vertex();
  Graph g2;
  auto rev2 = vector_property_map(0, get(edge_index, g2), Edge{});
  EXPECT_TRUE(read_dimacs_max_flow(g2, get(edge_capacity, g2), rev2.ref(), s2, t2, dimacs_stream));

  EXPECT_EQ(num_vertices(g), num_vertices(g2));
  EXPECT_EQ(num_edges(g) * 2, num_edges(g2));
  EXPECT_EQ(get(vertex_index, g, s), get(vertex_index, g2, s2));
  EXPECT_EQ(get(vertex_index, g, t), get(vertex_index, g2, t2));
}

}  // namespace
}  // namespace bagl
