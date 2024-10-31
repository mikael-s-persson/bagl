// Copyright Fernando Vilas 2012.
// Based on stoer_wagner_test.cpp by Daniel Trebbien.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/maximum_adjacency_search.h"

#include <array>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/connected_components.h"
#include "bagl/d_ary_heap.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/single_property_map.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using undirected_graph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_weight_t, int> >;
using undirected_unweighted_graph = adjacency_list<vec_s, vec_s, undirected_s>;
using vertex_descriptor = int;

std::string test_dir;

using edge_t = std::pair<std::size_t, std::size_t>;

template <typename KeyedUpdatablePriorityQueue>
class mas_test_visitor {
 public:
  typedef typename KeyedUpdatablePriorityQueue::key_type weight_type;

  explicit mas_test_visitor(KeyedUpdatablePriorityQueue& pq)
      : m_pq_(pq),
        vertex_visit_order_(std::make_shared<std::vector<vertex_descriptor> >()),
        vertex_weights_when_visited_(std::make_shared<std::vector<weight_type> >()) {}

  void clear() {
    vertex_visit_order_->clear();
    vertex_weights_when_visited_->clear();
  }

  template <typename Graph>
  void start_vertex(vertex_descriptor u, const Graph& /*g*/) {
    vertex_visit_order_->push_back(u);

    const weight_type u_weight = get(m_pq_.keys(), u);
    vertex_weights_when_visited_->push_back(u_weight);
  }

  const std::vector<vertex_descriptor>& vertex_visit_order() const { return *vertex_visit_order_; }

  const std::vector<weight_type>& vertex_weights_when_visited() const { return *vertex_weights_when_visited_; }

 private:
  const KeyedUpdatablePriorityQueue& m_pq_;
  std::shared_ptr<std::vector<vertex_descriptor> > vertex_visit_order_;
  std::shared_ptr<std::vector<weight_type> > vertex_weights_when_visited_;
};

// the example from Stoer & Wagner (1997)
// Check various implementations of the ArgPack where
// the weights are provided in it, and one case where
// they are not.
TEST(MaximumAdjacencySearchTest, Test0) {
  std::vector<edge_t> edge_list = {{0, 1}, {1, 2}, {2, 3}, {0, 4}, {1, 4}, {1, 5},
                                   {2, 6}, {3, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}};
  const std::vector<int> ws = {2, 3, 4, 3, 2, 2, 2, 2, 2, 3, 1, 3};
  const std::size_t vertices_count = 8;

  undirected_graph g(vertices_count, edge_list, ws);

  auto weights = get(edge_weight, g);
  auto assignments =
      vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<undirected_graph>::null_vertex());
  auto distances = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  auto indices_in_heap =
      vector_property_map(num_vertices(g), get(vertex_index, g), std::numeric_limits<std::size_t>::max());
  auto pq = make_d_ary_heap_indirect<vertex_descriptor, 22>(distances.ref(), indices_in_heap.ref(), std::greater<>());

  mas_test_visitor test_vis(pq);

  maximum_adjacency_search(g, weights, make_mas_visitor(test_vis), *vertices(g).begin(), assignments.ref(), pq);

  const std::vector<vertex_descriptor> expected_vertex_order1 = {0, 4, 1, 5, 2, 3, 6, 7};
  const std::vector<int> expected_weights_when_visited1 = {9, 3, 4, 5, 3, 4, 5, 5};

  EXPECT_THAT(test_vis.vertex_visit_order(), ::testing::ElementsAreArray(expected_vertex_order1));
  EXPECT_THAT(test_vis.vertex_weights_when_visited(), ::testing::ElementsAreArray(expected_weights_when_visited1));

  test_vis.clear();

  maximum_adjacency_search(g, weights, make_mas_visitor(test_vis), *vertices(g).begin(), assignments.ref(), pq);

  EXPECT_THAT(test_vis.vertex_visit_order(), ::testing::ElementsAreArray(expected_vertex_order1));
  EXPECT_THAT(test_vis.vertex_weights_when_visited(), ::testing::ElementsAreArray(expected_weights_when_visited1));

  test_vis.clear();

  maximum_adjacency_search(g, single_property_map(1), make_mas_visitor(test_vis), *vertices(g).begin(),
                           assignments.ref(), pq);

  const std::vector<vertex_descriptor> expected_vertex_order2 = {0, 1, 4, 5, 2, 6, 3, 7};
  const std::vector<int> expected_weights_when_visited2 = {9, 1, 2, 2, 1, 2, 2, 2};

  EXPECT_THAT(test_vis.vertex_visit_order(), ::testing::ElementsAreArray(expected_vertex_order2));
  EXPECT_THAT(test_vis.vertex_weights_when_visited(), ::testing::ElementsAreArray(expected_weights_when_visited2));
}

// Check the unweighted case
// with and without providing a weight_map
TEST(MaximumAdjacencySearchTest, Test1) {
  typedef graph_traits<undirected_unweighted_graph>::edge_descriptor edge_descriptor;

  std::vector<edge_t> edge_list = {{0, 1}, {0, 4}, {1, 2}, {1, 4}, {1, 5}, {2, 3},
                                   {2, 6}, {3, 6}, {3, 7}, {4, 5}, {5, 6}, {6, 7}};
  const std::size_t vertices_count = 8;

  undirected_unweighted_graph g(vertices_count, edge_list);

  auto assignments =
      vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<undirected_graph>::null_vertex());
  auto distances = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  auto indices_in_heap =
      vector_property_map(num_vertices(g), get(vertex_index, g), std::numeric_limits<std::size_t>::max());
  auto pq = make_d_ary_heap_indirect<vertex_descriptor, 22>(distances.ref(), indices_in_heap.ref(), std::greater<>());

  mas_test_visitor test_vis(pq);

  maximum_adjacency_search(g, single_property_map(1), make_mas_visitor(test_vis), *vertices(g).begin(),
                           assignments.ref(), pq);

  const std::vector<vertex_descriptor> expected_vertex_order1 = {0, 1, 4, 5, 2, 6, 3, 7};
  const std::vector<int> expected_weights_when_visited1 = {9, 1, 2, 2, 1, 2, 2, 2};

  EXPECT_THAT(test_vis.vertex_visit_order(), ::testing::ElementsAreArray(expected_vertex_order1));
  EXPECT_THAT(test_vis.vertex_weights_when_visited(), ::testing::ElementsAreArray(expected_weights_when_visited1));

  test_vis.clear();

  const std::vector<int> ws = {2, 3, 3, 2, 2, 4, 2, 2, 2, 3, 1, 3};
  std::map<edge_descriptor, int> wm;
  int i = 0;
  for (auto e : edges(g)) {
    wm[e] = ws[i];
    ++i;
  }
  associative_property_map ws_map(wm);

  maximum_adjacency_search(g, ws_map, make_mas_visitor(test_vis), *vertices(g).begin(), assignments.ref(), pq);

  const std::vector<vertex_descriptor> expected_vertex_order2 = {0, 4, 1, 5, 2, 3, 6, 7};
  const std::vector<int> expected_weights_when_visited2 = {9, 3, 4, 5, 3, 4, 5, 5};

  EXPECT_THAT(test_vis.vertex_visit_order(), ::testing::ElementsAreArray(expected_vertex_order2));
  EXPECT_THAT(test_vis.vertex_weights_when_visited(), ::testing::ElementsAreArray(expected_weights_when_visited2));
}

void test_weighted(const std::vector<edge_t>& edge_list, const std::vector<int> weights_list,
                   const std::vector<vertex_descriptor>& expected_vertex_order,
                   const std::vector<int>& expected_weights_when_visited, vertex_descriptor start_vertex = 0) {
  const undirected_unweighted_graph g(expected_vertex_order.size(), edge_list);

  auto assignments =
      vector_property_map(num_vertices(g), get(vertex_index, g), graph_traits<undirected_graph>::null_vertex());
  auto distances = vector_property_map(num_vertices(g), get(vertex_index, g), int{0});
  auto indices_in_heap =
      vector_property_map(num_vertices(g), get(vertex_index, g), std::numeric_limits<std::size_t>::max());
  auto pq = make_d_ary_heap_indirect<vertex_descriptor, 4>(distances.ref(), indices_in_heap.ref(), std::greater<>());
  mas_test_visitor test_vis(pq);

  std::map<graph_edge_descriptor_t<undirected_unweighted_graph>, int> wm;
  int i = 0;
  for (auto e : edges(g)) {
    wm[e] = weights_list[i];
    ++i;
  }
  associative_property_map ws_map(wm);

  maximum_adjacency_search(g, ws_map, make_mas_visitor(test_vis), start_vertex, assignments.ref(), pq);

  EXPECT_THAT(test_vis.vertex_visit_order(), ::testing::ElementsAreArray(expected_vertex_order));
  EXPECT_THAT(test_vis.vertex_weights_when_visited(), ::testing::ElementsAreArray(expected_weights_when_visited));
}

void test_unweighted(const std::vector<edge_t>& edge_list, const std::vector<vertex_descriptor>& expected_vertex_order,
                     const std::vector<int>& expected_weights_when_visited, vertex_descriptor start_vertex = 0) {
  std::vector<int> weights_list;
  weights_list.resize(edge_list.size());
  for (std::size_t i = 0; i < edge_list.size(); i++) {
    weights_list[i] = 1;
  }

  test_weighted(edge_list, weights_list, expected_vertex_order, expected_weights_when_visited, start_vertex);
}

TEST(MaximumAdjacencySearchTest, Test2) {
  const std::size_t vertices_count = 2;

  const std::vector<edge_t> edge_list = {{{0, 1}}};
  const std::vector<vertex_descriptor> expected_vertex_order = {0, 1};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 1};

  test_unweighted(edge_list, expected_vertex_order, expected_weights_when_visited);
}

TEST(MaximumAdjacencySearchTest, Test3) {
  const std::size_t vertices_count = 3;

  const std::vector<edge_t> edge_list = {{{0, 1}, {1, 2}}};
  const std::vector<vertex_descriptor> expected_vertex_order = {0, 1, 2};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 1, 1};

  test_unweighted(edge_list, expected_vertex_order, expected_weights_when_visited);
}

TEST(MaximumAdjacencySearchTest, Test4) {
  const std::size_t vertices_count = 3;

  const std::vector<edge_t> edge_list = {{{0, 1}, {0, 2}, {1, 2}}};
  const std::vector<vertex_descriptor> expected_vertex_order = {0, 1, 2};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 1, 2};

  test_unweighted(edge_list, expected_vertex_order, expected_weights_when_visited);
}

// The example graph from Matula (1993)
TEST(MaximumAdjacencySearchTest, Test5Matula1993) {
  const std::size_t vertices_count = 12;

  const std::vector<edge_t> edge_list = {{{0, 1},  {0, 2}, {0, 3}, {0, 9},  {1, 2},  {1, 4},  {1, 10}, {2, 5},
                                          {2, 11}, {3, 4}, {3, 5}, {3, 6},  {4, 5},  {4, 7},  {5, 8},  {6, 7},
                                          {6, 8},  {6, 9}, {7, 8}, {7, 10}, {8, 11}, {9, 10}, {9, 11}, {10, 11}}};
  const std::vector<vertex_descriptor> expected_vertex_order = {0, 1, 2, 10, 9, 11, 6, 3, 7, 4, 8, 5};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 1, 2, 1, 2, 3, 1, 2, 2, 3, 3, 4};

  test_unweighted(edge_list, expected_vertex_order, expected_weights_when_visited);
}

// Testing with a different start vertex
TEST(MaximumAdjacencySearchTest, Test6StartVertex) {
  const std::size_t vertices_count = 3;

  const std::vector<edge_t> edge_list = {{{0, 1}, {1, 2}}};
  const std::vector<vertex_descriptor> expected_vertex_order = {1, 0, 2};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 1, 1};

  test_unweighted(edge_list, expected_vertex_order, expected_weights_when_visited, 1);
}

TEST(MaximumAdjacencySearchTest, Test7) {
  const std::size_t vertices_count = 3;

  const std::vector<edge_t> edge_list = {{{0, 1}, {1, 2}}};
  const std::vector<int> weights_list = {2, 6};
  const std::vector<vertex_descriptor> expected_vertex_order = {0, 1, 2};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 2, 6};

  test_weighted(edge_list, weights_list, expected_vertex_order, expected_weights_when_visited);
}

TEST(MaximumAdjacencySearchTest, Test8) {
  const std::size_t vertices_count = 3;

  const std::vector<edge_t> edge_list = {{{0, 1}, {0, 2}, {1, 2}}};
  const std::vector<int> weights_list = {2, 6, 7};
  const std::vector<vertex_descriptor> expected_vertex_order = {0, 2, 1};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 6, 9};

  test_weighted(edge_list, weights_list, expected_vertex_order, expected_weights_when_visited);
}

TEST(MaximumAdjacencySearchTest, Test9StartVertex) {
  const std::size_t vertices_count = 3;

  const std::vector<edge_t> edge_list = {{{0, 1}, {0, 2}, {1, 2}}};
  const std::vector<int> weights_list = {2, 6, 7};
  const std::vector<vertex_descriptor> expected_vertex_order = {1, 2, 0};
  const std::vector<int> expected_weights_when_visited = {vertices_count + 1, 7, 8};

  test_weighted(edge_list, weights_list, expected_vertex_order, expected_weights_when_visited, 1);
}

}  // namespace
}  // namespace bagl
