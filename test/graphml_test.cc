// Copyright (C) 2006  Tiago de Paula Peixoto <tiago@forked.de>
//  Authors: Tiago de Paula Peixoto
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/graphml.h"

#include <cmath>
#include <fstream>
#include <string>

#include "bagl/adjacency_list.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(GraphMLTest, RoundTrip) {
  using Graph =
      adjacency_list<vec_s, vec_s, directed_s, property<vertex_color_t, int, property<vertex_name_t, std::string> >,
                     property<edge_weight_t, double>, property<graph_name_t, std::string> >;
  Graph g;
  dynamic_properties dp;
  dp.property("foo", vertex_color, g);
  dp.property("weight", edge_weight, g);
  dp.property("name", vertex_name, g);
  dp.property("description", graph_name, g);

  std::ifstream ifile("test/testdata/graphml_test.xml");
  read_graphml(ifile, g, dp);
  ifile.close();

  EXPECT_EQ(num_vertices(g), 9);
  EXPECT_EQ(num_edges(g), 9);
  EXPECT_EQ(get(vertex_color, g, vertex(2, g)), 100);
  EXPECT_EQ(get(vertex_color, g, vertex(3, g)), 42);
  EXPECT_THAT(get(edge_weight, g, edge(vertex(0, g), vertex(1, g), g).first), ::testing::DoubleNear(0.0, 0.00001));
  EXPECT_THAT(get(edge_weight, g, edge(vertex(1, g), vertex(2, g), g).first), ::testing::DoubleNear(0.8, 0.00001));
  EXPECT_EQ(get("description", dp, &g), "Root graph.");

  std::stringstream outstr;
  write_graphml(outstr, g, dp);

  Graph g2;
  dynamic_properties dp2;
  dp2.property("foo", vertex_color, g2);
  dp2.property("weight", edge_weight, g2);
  dp2.property("name", vertex_name, g2);
  dp2.property("description", graph_name, g2);

  read_graphml(outstr, g2, dp2);

  EXPECT_EQ(num_vertices(g), num_vertices(g2));
  EXPECT_EQ(num_edges(g), num_edges(g2));
  EXPECT_EQ(get("description", dp, &g), get("description", dp2, &g2));

  for (auto v : vertices(g)) {
    EXPECT_EQ(get(vertex_color, g, v), get(vertex_color, g2, v));
  }

  for (auto e : edges(g)) {
    EXPECT_THAT(get(edge_weight, g, e), ::testing::DoubleNear(get(edge_weight, g2, e), 0.00001));
  }
}

}  // namespace
}  // namespace bagl
