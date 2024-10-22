//  Copyright (c) 2006, Stephan Diederich
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/boykov_kolmogorov_max_flow.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <type_traits>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph, typename CapacityMap, typename ReverseEdgeMap, typename EdgeIndexMap, typename VertexIndexMap>
auto fill_random_max_flow_graph(Graph& g, CapacityMap cap, ReverseEdgeMap rev, EdgeIndexMap eindex,
                                VertexIndexMap vindex, std::size_t n_verts, std::size_t n_edges, std::size_t seed) {
  using Edge = graph_edge_descriptor_t<Graph>;
  const int cap_low = 1;
  const int cap_high = 1000;

  // init random numer generator
  std::mt19937 gen(seed);
  // generate graph
  generate_random_graph(g, n_verts, n_edges, gen);

  // initialize vertex indices, if needed.
  if constexpr (!std::is_integral_v<graph_vertex_descriptor_t<Graph>>) {
    std::size_t index = 0;
    for (auto v : vertices(g)) {
      put(vindex, v, index++);
    }
  }

  // init an uniform distribution int generator
  std::uniform_int_distribution<int> int_gen(cap_low, cap_high);
  // randomize edge-capacities
  std::size_t e_id = 0;
  for (auto e : edges(g)) {
    put(eindex, e, e_id++);
    put(cap, e, int_gen(gen));
  }

  // get source and sink node
  auto s = random_vertex(g, gen);
  auto t = random_vertex(g, gen);
  while (t == s) {
    t = random_vertex(g, gen);
  }

  // add reverse edges (ugly... how to do better?!)
  std::deque<Edge> edges_copy;
  std::ranges::copy(edges(g), std::back_inserter(edges_copy));
  while (!edges_copy.empty()) {
    Edge old_edge = edges_copy.front();
    edges_copy.pop_front();
    auto [new_edge, inserted] = add_edge(target(old_edge, g), source(old_edge, g), g);
    assert(inserted);
    put(eindex, new_edge, e_id++);
    put(rev, old_edge, new_edge);
    put(rev, new_edge, old_edge);
    put(cap, new_edge, 0);
  }
  return std::pair(s, t);
}

template <typename Graph, typename EdgeCapacityMap, typename ResidualCapacityEdgeMap, typename ReverseEdgeMap,
          typename PredecessorMap, typename ColorMap, typename DistanceMap, typename IndexMap>
class bk_max_flow_tester
    : public max_flow_detail::bk_max_flow<Graph, EdgeCapacityMap, ResidualCapacityEdgeMap, ReverseEdgeMap,
                                          PredecessorMap, ColorMap, DistanceMap, IndexMap> {
  using Edge = graph_edge_descriptor_t<Graph>;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using EdgeValue = property_traits_value_t<EdgeCapacityMap>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using ColorTraits = color_traits<ColorValue>;
  using DistanceValue = property_traits_value_t<DistanceMap>;
  using Super = max_flow_detail::bk_max_flow<Graph, EdgeCapacityMap, ResidualCapacityEdgeMap, ReverseEdgeMap,
                                             PredecessorMap, ColorMap, DistanceMap, IndexMap>;

 public:
  bk_max_flow_tester(Graph& g, EdgeCapacityMap cap, ResidualCapacityEdgeMap res_cap, ReverseEdgeMap rev,
                     PredecessorMap pred, ColorMap color, DistanceMap dist, IndexMap vindex, Vertex src, Vertex sink)
      : Super(g, cap, res_cap, rev, pred, color, dist, vindex, src, sink) {}

  void invariant_four(Vertex v) {
    // passive nodes in S or T
    if (v == this->source_ || v == this->sink_) {
      return;
    }
    auto it = std::find(this->orphans_.begin(), this->orphans_.end(), v);
    // a node is active, if its in the active_list AND (is has_a_parent, or
    // its already in the orphans_list or its the sink, or its the source)
    bool is_active = (this->in_active_list_map_[v] && (this->has_parent(v) || it != this->orphans_.end()));
    if (this->get_tree(v) != ColorTraits::gray() && !is_active) {
      for (auto e : out_edges(v, this->g_)) {
        auto other_node = target(e, this->g_);
        if (this->get_tree(other_node) != this->get_tree(v)) {
          if (this->get_tree(v) == ColorTraits::black()) {
            EXPECT_EQ(this->res_cap_map_[e], 0);
          } else {
            EXPECT_EQ(this->res_cap_map_[this->rev_edge_map_[e]], 0);
          }
        }
      }
    }
  }

  void invariant_five(const Vertex& v) {
    auto tree_color = this->get_tree(v);
    if (tree_color == ColorTraits::gray()) {
      EXPECT_LE(this->time_map_[v], this->time_);
    }
  }

  void invariant_six(const Vertex& v) {
    if (this->get_tree(v) == ColorTraits::gray() || this->time_map_[v] != this->time_) {
      return;
    }

    Vertex current_node = v;
    DistanceValue distance = 0;
    ColorValue color = this->get_tree(v);
    Vertex terminal = (color == ColorTraits::black()) ? this->source_ : this->sink_;
    while (current_node != terminal) {
      EXPECT_TRUE(this->has_parent(current_node));
      Edge e = this->get_edge_to_parent(current_node);
      ++distance;
      current_node = (color == ColorTraits::black()) ? source(e, this->g_) : target(e, this->g_);
      if (distance > this->dist_map_[v]) {
        break;
      }
    }
    EXPECT_EQ(distance, this->dist_map_[v]);
  }

  void invariant_seven(const Vertex& v) {
    if (this->get_tree(v) == ColorTraits::gray()) {
      return;
    }
    ColorValue color = this->get_tree(v);
    auto time = this->time_map_[v];
    Vertex current_node = v;
    while (this->has_parent(current_node)) {
      Edge e = this->get_edge_to_parent(current_node);
      current_node = (color == ColorTraits::black()) ? source(e, this->g_) : target(e, this->g_);
      EXPECT_GE(this->time_map_[current_node], time);
    }
  }

  void invariant_eight(const Vertex& v) {
    if (this->get_tree(v) == ColorTraits::gray()) {
      return;
    }
    ColorValue color = this->get_tree(v);
    auto time = this->time_map_[v];
    DistanceValue distance = this->dist_map_[v];
    Vertex current_node = v;
    while (this->has_parent(current_node)) {
      Edge e = this->get_edge_to_parent(current_node);
      current_node = (color == ColorTraits::black()) ? source(e, this->g_) : target(e, this->g_);
      if (this->time_map_[current_node] == time) {
        EXPECT_LT(this->dist_map_[current_node], distance);
      }
    }
  }

  void check_invariants() {
    for (auto v : vertices(this->g_)) {
      invariant_four(v);
      invariant_five(v);
      invariant_six(v);
      invariant_seven(v);
      invariant_eight(v);
    }
  }

  EdgeValue test() {
    this->add_active_node(this->sink_);
    this->augment_direct_paths();
    check_invariants();
    // start the main-loop
    while (true) {
      // find a path from source to sink
      auto [connecting_edge, path_found] = this->grow();
      if (!path_found) {
        // we're finished, no more paths were found
        break;
      }
      check_invariants();
      this->time_++;
      this->augment(connecting_edge);  // augment that path
      check_invariants();
      this->adopt();  // rebuild search tree structure
      check_invariants();
    }

    // check if flow is the sum of outgoing edges of src
    EdgeValue src_sum = 0;
    for (auto e : out_edges(this->source_, this->g_)) {
      src_sum += this->cap_map_[e] - this->res_cap_map_[e];
    }
    EXPECT_EQ(this->flow_, src_sum);
    // check if flow is the sum of ingoing edges of sink
    EdgeValue sink_sum = 0;
    for (auto e : out_edges(this->sink_, this->g_)) {
      Edge in_edge = this->rev_edge_map_[e];
      sink_sum += this->cap_map_[in_edge] - this->res_cap_map_[in_edge];
    }
    EXPECT_EQ(this->flow_, sink_sum);
    return this->flow_;
  }
};

class BoykoyKolmogorovTest : public ::testing::Test {
 public:
  static constexpr std::size_t n_verts = 10;
  static constexpr std::size_t n_edges = 500;
  static constexpr std::size_t seed = 3453;

  using AdjVecGraph = adjacency_list<
      vec_s, vec_s, directed_s, property<vertex_color_t, default_color_type, property<vertex_distance_t, int>>,
      property<edge_index_t, std::size_t, property<edge_capacity_t, int, property<edge_residual_capacity_t, int>>>>;

  using AdjListGraph = adjacency_list<
      list_s, list_s, directed_s,
      property<vertex_index_t, std::size_t,
               property<vertex_color_t, default_color_type, property<vertex_distance_t, int>>>,
      property<edge_index_t, std::size_t, property<edge_capacity_t, int, property<edge_residual_capacity_t, int>>>>;

  struct Node {
    default_color_type color;
    int distance;
  };
  struct Link {
    std::size_t index;
    int capacity;
    int residual_capacity;
  };
  using BundleGraph = adjacency_list<vec_s, vec_s, directed_s, Node, Link>;

  void SetUp() override {
    AdjVecGraph g;
    auto rev = vector_property_map(2 * n_edges, get(edge_index, g), graph_edge_descriptor_t<AdjVecGraph>{});

    auto [src, sink] = fill_random_max_flow_graph(g, get(edge_capacity, g), rev.ref(), get(edge_index, g),
                                                  get(vertex_index, g), n_verts, n_edges, seed);

    auto pred = vector_property_map(n_verts, get(vertex_index, g), graph_edge_descriptor_t<AdjVecGraph>{});
    expected_max_flow_ =
        boykov_kolmogorov_max_flow(g, get(edge_capacity, g), get(edge_residual_capacity, g), rev.ref(), pred.ref(),
                                   get(vertex_color, g), get(vertex_distance, g), get(vertex_index, g), src, sink);
  }

 protected:
  int expected_max_flow_ = 0;
};

TEST_F(BoykoyKolmogorovTest, AdjListGraphTest) {
  AdjListGraph g;
  auto rev = vector_property_map(2 * n_edges, get(edge_index, g), graph_edge_descriptor_t<AdjListGraph>{});

  auto [src, sink] = fill_random_max_flow_graph(g, get(edge_capacity, g), rev.ref(), get(edge_index, g),
                                                get(vertex_index, g), n_verts, n_edges, seed);

  auto pred = vector_property_map(n_verts, get(vertex_index, g), graph_edge_descriptor_t<AdjListGraph>{});
  EXPECT_EQ(boykov_kolmogorov_max_flow(g, get(edge_capacity, g), get(edge_residual_capacity, g), rev.ref(), pred.ref(),
                                       get(vertex_color, g), get(vertex_distance, g), get(vertex_index, g), src, sink),
            expected_max_flow_);
}

TEST_F(BoykoyKolmogorovTest, BundleGraphTest) {
  BundleGraph g;
  auto rev = vector_property_map(2 * n_edges, get(&Link::index, g), graph_edge_descriptor_t<BundleGraph>{});

  auto [src, sink] = fill_random_max_flow_graph(g, get(&Link::capacity, g), rev.ref(), get(&Link::index, g),
                                                get(vertex_index, g), n_verts, n_edges, seed);

  auto pred = vector_property_map(n_verts, get(vertex_index, g), graph_edge_descriptor_t<BundleGraph>{});
  EXPECT_EQ(
      boykov_kolmogorov_max_flow(g, get(&Link::capacity, g), get(&Link::residual_capacity, g), rev.ref(), pred.ref(),
                                 get(&Node::color, g), get(&Node::distance, g), get(vertex_index, g), src, sink),
      expected_max_flow_);
}

TEST_F(BoykoyKolmogorovTest, InvariantsTest) {
  AdjVecGraph g;
  auto rev = vector_property_map(2 * n_edges, get(edge_index, g), graph_edge_descriptor_t<AdjVecGraph>{});

  auto [src, sink] = fill_random_max_flow_graph(g, get(edge_capacity, g), rev.ref(), get(edge_index, g),
                                                get(vertex_index, g), n_verts, n_edges, seed);

  auto pred = vector_property_map(n_verts, get(vertex_index, g), graph_edge_descriptor_t<AdjVecGraph>{});
  bk_max_flow_tester tester(g, get(edge_capacity, g), get(edge_residual_capacity, g), rev.ref(), pred.ref(),
                            get(vertex_color, g), get(vertex_distance, g), get(vertex_index, g), src, sink);
  EXPECT_EQ(tester.test(), expected_max_flow_);
}

}  // namespace
}  // namespace bagl
