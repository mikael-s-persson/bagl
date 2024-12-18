// Copyright (C) 2006-2009 Dmitry Bufistov and Andrey Parfenov
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/howard_cycle_ratio.h"

#include <sstream>
#include <string_view>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/graphviz.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// The graph has two equal cycles with ratio 2/3
const std::string_view test_graph1 = R"""(digraph G {
  edge [w1=1, w2=1];
  1->2;
  2->3 [w1=0];
  3->4;
  4->2;
  1->5;
  5->6;
  6->7 [w1=0];
  7->5;
})""";

// The graph has no cycles
const std::string_view test_graph2 = R"""(digraph G {
  edge [w1=1];
  1->3 [w2=1];
  1->2 [w2=2];
  1->4 [w2=7];
})""";

// Example from the paper "Nunerical computation of spectral elements"
// Maximum cycle ratio is 5.5
const std::string_view test_graph3 = R"""(digraph article {
  edge [w2 =2];
  1->1 [w1 = 1];
  1->2 [w1 = 2];
  1->4 [w1 = 7];
  2->2 [w1 = 3];
  2->3 [w1 = 5];
  3->2 [w1 = 4];
  3->4 [w1 = 3];
  4->2 [w1 = 2];
  4->3 [w1 = 8];
})""";

// Simple multigraph.
// Maximum cycle ratio is 2.5, minimum  0.5
const std::string_view test_graph4 = R"""(digraph G {
  edge [w2=1];
  a->b  [w1=1];
  b->a  [w1=0];
  a->b [w1=2];
  b->a [w1=3];
})""";

// The big graph with two equal cycles
const std::string_view test_graph5 = R"""(digraph G {
  edge [w2=1, w1=1];
  n94->n8;
  n95->n8;
  n93->n8;
  n93->n9;
  n42->n9;
  n23->n13;
  n29->n13;
  n95->n14;
  n37->n14;
  n95->n19;
  n37->n19;
  n94->n23;
  n60->n26;
  n76->n26;
  n94->n29;
  n9->n33 [w1=0];
  n80->n33;
  n14->n34 [w1=0];
  n19->n34;
  n94->n37;
  n94->n42;
  n95->n42;
  n8->n60;
  n26->n60;
  n26->n76;
  n106->n76;
  n93->n80;
  n42->n80;
  n33->n93;
  n60->n93;
  n13->n94;
  n60->n94;
  n34->n95;
  n60->n95;
  n94->n106;
  n95->n106;
  n93->n106;
})""";

// Random graph generated by hands.
// Maximum cycle ratio is 3.58, minimum is 0.294118
const std::string_view test_graph6 = R"""(digraph test_graph6 {
  16;
  17;

  1->2 [w1=1, w2=0.1];
  2->3 [w1 = 2, w2=3.6];
  3->4 [w1=7, w2=8];
  4->5 [w1=3.1,w2=0.8];
  4->5 [w1 = 4.2, w2=0.6];
  4->5 [w1 = 5.3, w2=0.4];
  5->6 [w1=-10, w2 = 34.75];
  6->1 [w1=100, w2 = 20];

  1->7 [w1=10, w2 = 20];
  7->8 [w1=3.75, w2 = 1.25];
  7->8 [w1=30, w2 = 22.2];
  8->9 [w1=10, w2 = 20];
  9->10 [w1=-2.1, w2 = 30];
  10->6 [w1=10, w2 = 20];

  11->12 [w1 = 5, w2=0.45];
  12->13 [w1 = 4, w2=0.2];
  13->14 [w1 = 3, w2=15.75];
  14->11 [w1 = -2.5, w2=0.6];
  11->10 [w1 = -8, w2=0.9];
  11->10 [w1 = -15, w2=2.9];

  18 -> 19 [w1=18, w2=6];
  18 -> 20 [w1=16.3, w2=8.2];
  18 -> 21 [w1=-3, w2=15];
  18 -> 18 [w1=2, w2=1];
  19 -> 18 [w1=0.06, w2=0.01];
  19 -> 19 [w1=1, w2=1.2];
  19 -> 20 [w1=5, w2=2];
  19 -> 21 [w1=3, w2=0.1];
  20 -> 18 [w1=4, w2=0.2];
  20 -> 19 [w1=11, w2=21];
  20 -> 20 [w1=6, w2=5];
  20 -> 21 [w1=7, w2=1];
  21 -> 18 [w1=8, w2=2];
  21 -> 19 [w1=12, w2=6];
  21 -> 20 [w1=7.5, w2=4.3];
  21 -> 21 [w1=1.25, w2=2.15];
})""";

constexpr double epsilon = 0.005;

template <typename TW>
struct Graph {
  using VertexProps = property<vertex_index_t, std::size_t, property<vertex_name_t, std::string> >;
  using EdgeProps = property<edge_weight_t, TW, property<edge_weight2_t, TW, property<edge_index_t, std::size_t> > >;
  using type = adjacency_list<list_s, list_s, directed_s, VertexProps, EdgeProps>;
};
using DiGraphInt = Graph<int>::type;
using DiGraphReal = Graph<double>::type;

using CCInt = std::vector<graph_edge_descriptor_t<DiGraphInt> >;
using CCReal = std::vector<graph_edge_descriptor_t<DiGraphReal> >;

template <typename TW>
struct CEdgeProps {
  explicit CEdgeProps(TW a_w1 = 1, TW a_w2 = 2) : w1(a_w1), w2(a_w2), eindex(std::numeric_limits<std::size_t>::max()) {}
  TW w1;
  TW w2;
  std::size_t eindex;
};
using GraphMInt = adjacency_matrix<directed_s, no_property, CEdgeProps<int> >;

/// Create "tokens_map" for reading graph properties from .dot file
template <typename TG>
void make_dynamic_properties(TG& g, dynamic_properties& p) {
  p.property("node_id", vertex_name, g);
  p.property("label", edge_weight, g);
  p.property("w1", edge_weight, g);
  p.property("w2", edge_weight2, g);
}

template <typename TG>
void read_data(std::istream& is, TG& g) {
  dynamic_properties p;
  make_dynamic_properties(g, p);
  read_graphviz(is, g, p);
  int i = 0;
  for (auto vd : vertices(g)) {
    put(vertex_index, g, vd, i++);
  }
  i = 0;
  for (auto ed : edges(g)) {
    put(edge_index, g, ed, i++);
  }
}

struct my_float : mcr_float<> {
  static double infinity() { return 1000; }
};

TEST(HowardCycleRatio, TestGraph1) {
  DiGraphInt tg;
  auto vim = get(vertex_index, tg);
  auto ew1m = get(edge_weight, tg);
  auto ew2m = get(edge_weight2, tg);

  std::stringstream iss;
  iss << test_graph1;
  read_data(iss, tg);
  EXPECT_EQ(num_vertices(tg), 7);
  EXPECT_EQ(num_edges(tg), 8);
  auto max_cr = maximum_cycle_ratio(tg, vim, ew1m, ew2m);
  EXPECT_THAT(max_cr, ::testing::DoubleNear(0.666666666, epsilon));
}

TEST(HowardCycleRatio, TestGraph2) {
  DiGraphInt tg;
  auto vim = get(vertex_index, tg);
  auto ew1m = get(edge_weight, tg);
  auto ew2m = get(edge_weight2, tg);

  std::stringstream iss;
  iss << test_graph2;
  read_data(iss, tg);
  EXPECT_EQ(num_vertices(tg), 4);
  EXPECT_EQ(num_edges(tg), 3);
  auto max_cr = maximum_cycle_ratio(tg, vim, ew1m, ew2m, nullptr, my_float());
  EXPECT_THAT(max_cr, ::testing::DoubleNear(-1000.0, epsilon));
}

TEST(HowardCycleRatio, TestGraph3) {
  DiGraphInt tg;
  auto vim = get(vertex_index, tg);
  auto ew1m = get(edge_weight, tg);
  auto ew2m = get(edge_weight2, tg);

  std::stringstream iss;
  iss << test_graph3;
  read_data(iss, tg);
  EXPECT_EQ(num_vertices(tg), 4);
  EXPECT_EQ(num_edges(tg), 9);
  auto max_cr = maximum_cycle_ratio(tg, vim, ew1m, ew2m);
  EXPECT_THAT(max_cr, ::testing::DoubleNear(2.75, epsilon));
  double maxmc = maximum_cycle_mean(tg, vim, ew1m);
  EXPECT_THAT(maxmc, ::testing::DoubleNear(5.5, epsilon));
}

TEST(HowardCycleRatio, TestGraph4) {
  DiGraphInt tg;
  auto vim = get(vertex_index, tg);
  auto ew1m = get(edge_weight, tg);
  auto ew2m = get(edge_weight2, tg);

  std::stringstream iss;
  iss << test_graph4;
  read_data(iss, tg);
  EXPECT_EQ(num_vertices(tg), 2);
  EXPECT_EQ(num_edges(tg), 4);
  auto max_cr = maximum_cycle_ratio(tg, vim, ew1m, ew2m);
  EXPECT_THAT(max_cr, ::testing::DoubleNear(2.5, epsilon));
  auto min_cr = minimum_cycle_ratio(tg, vim, ew1m, ew2m);
  EXPECT_THAT(min_cr, ::testing::DoubleNear(0.5, epsilon));
}

TEST(HowardCycleRatio, TestGraph5) {
  DiGraphInt tg;
  auto vim = get(vertex_index, tg);
  auto ew1m = get(edge_weight, tg);
  auto ew2m = get(edge_weight2, tg);

  std::stringstream iss;
  iss << test_graph5;
  read_data(iss, tg);
  EXPECT_EQ(num_vertices(tg), 19);
  EXPECT_EQ(num_edges(tg), 37);
  CCInt cc;  /// critical cycle
  auto min_cr = minimum_cycle_ratio(tg, vim, ew1m, ew2m, &cc, my_float());
  EXPECT_THAT(min_cr, ::testing::DoubleNear(0.6666666666, epsilon));
  EXPECT_THAT(cc, ::testing::SizeIs(3));
  std::pair<double, double> cr(.0, .0);
  for (auto c : cc) {
    cr.first += get(ew1m, c);
    cr.second += get(ew2m, c);
  }
  EXPECT_THAT(cr.first / cr.second, ::testing::DoubleNear(min_cr, epsilon))
      << "From cycle: " << cr.first << " / " << cr.second;
}

TEST(HowardCycleRatio, TestGraph6) {
  DiGraphReal tg;
  auto vim = get(vertex_index, tg);
  auto ew1m = get(edge_weight, tg);
  auto ew2m = get(edge_weight2, tg);
  CCReal cc;
  std::stringstream iss;
  iss << test_graph6;
  read_data(iss, tg);
  EXPECT_EQ(num_vertices(tg), 20);
  EXPECT_EQ(num_edges(tg), 36);
  auto max_cr = maximum_cycle_ratio(tg, vim, ew1m, ew2m);
  EXPECT_THAT(max_cr, ::testing::DoubleNear(3.58025, epsilon));
  auto min_cr = minimum_cycle_ratio(tg, vim, ew1m, ew2m, &cc);
  EXPECT_THAT(min_cr, ::testing::DoubleNear(0.294118, epsilon));
  EXPECT_THAT(cc, ::testing::SizeIs(2));
  std::pair<double, double> cr(.0, .0);
  for (auto c : cc) {
    cr.first += get(ew1m, c);
    cr.second += get(ew2m, c);
  }
  EXPECT_THAT(cr.first / cr.second, ::testing::DoubleNear(min_cr, epsilon))
      << "From cycle: " << cr.first << " / " << cr.second;
}

TEST(HowardCycleRatio, TestMatrix) {
  GraphMInt gm(10);
  auto v_rg = vertices(gm);
  for (auto vi1 = v_rg.begin(), vi_end = v_rg.end(); vi1 != vi_end; ++vi1) {
    for (auto vi2 = v_rg.begin(); vi2 != vi_end; ++vi2) {
      add_edge(*vi1, *vi2, gm);
    }
  }

  auto vim = get(vertex_index, gm);
  auto ew1m = get(&CEdgeProps<int>::w1, gm);
  auto ew2m = get(&CEdgeProps<int>::w2, gm);
  auto max_cr = maximum_cycle_ratio(gm, vim, ew1m, ew2m);
  EXPECT_THAT(max_cr, ::testing::DoubleNear(0.5, epsilon));
}

}  // namespace
}  // namespace bagl
