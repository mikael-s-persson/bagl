// Copyright 2010 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/random_spanning_tree.h"

#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "bagl/dynamic_property_map.h"
#include "bagl/graph_traits.h"
#include "bagl/graphviz.h"
#include "bagl/grid_graph.h"
#include "bagl/properties.h"
#include "bagl/random.h"
#include "bagl/single_property_map.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph, typename PredMap, typename WeightMap>
void write_spanning_tree(const Graph& g, PredMap pred, WeightMap weight, std::string filename) {
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using Edge = graph_edge_descriptor_t<Graph>;
  auto edge_style = vector_property_map(num_edges(g), get(edge_index, g), std::string{});
  auto vertex_pos = vector_property_map(num_vertices(g), get(vertex_index, g), std::string{});
  for (auto e : edges(g)) {
    put(edge_style, e,
        (get(pred, target(e, g)) == source(e, g) || get(pred, source(e, g)) == target(e, g)) ? "bold" : "dotted");
  }
  for (auto v : vertices(g)) {
    put(vertex_pos, v, std::to_string(v[0]) + "," + std::to_string(v[1]));
  }
  dynamic_properties dp;
  dp.property<Edge>("style", edge_style.ref());
  dp.property<Vertex>("pos", vertex_pos.ref());
  dp.property<Edge>("len", weight);
  dp.property("node_id", vertex_index, g);
  std::ofstream out(filename.c_str());
  write_graphviz(out, g, dp);
}

bool print_out_test_graphs = false;

TEST(RandomSpanningTreeTest, Basic) {
  using Graph = grid_graph<2>;

  std::array<std::size_t, 2> sizes = {{5, 5}};
  Graph g(sizes);

  auto pred = vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<Graph>::null_vertex());
  auto weight = vector_property_map(num_edges(g), get(edge_index, g), double{1.0});

  for (auto e : edges(g)) {
    put(weight, e, (1. + get(edge_index, g, e)) / num_edges(g));
  }

  std::mt19937 gen{42};
  random_spanning_tree(g, gen, random_vertex(g, gen), pred.ref());
  if (print_out_test_graphs) {
    // write_spanning_tree(g, pred.ref(), single_property_map(1.), "unweight_random_st.dot");
  }
  random_spanning_tree(g, gen, random_vertex(g, gen), pred.ref());
  if (print_out_test_graphs) {
    // write_spanning_tree(g, pred.ref(), single_property_map(1.), "unweight_random_st2.dot");
  }

  random_spanning_tree(g, gen, random_vertex(g, gen), pred.ref(), weight.ref());
  if (print_out_test_graphs) {
    // write_spanning_tree(g, pred.ref(), weight.ref(), "weight_random_st.dot");
  }
}

}  // namespace
}  // namespace bagl
