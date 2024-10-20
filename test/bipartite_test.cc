/// Copyright (c) 2010 Matthias Walter (xammy@xammy.homelinux.net)
// Authors: Matthias Walter
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/bipartite.h"

#include <ranges>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/lookup_edge.h"
#include "bagl/property.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename Graph, typename ColorMap>
void check_two_coloring(const Graph& g, const ColorMap color_map) {
  for (auto e : edges(g)) {
    EXPECT_NE(get(color_map, source(e, g)), get(color_map, target(e, g)));
  }
}

/// Tests for a vertex sequence to define an odd cycle

template <typename Graph, typename RandomAccessIterator>
void check_odd_cycle(const Graph& g, RandomAccessIterator first, RandomAccessIterator beyond) {
  EXPECT_EQ((beyond - first) % 2, 1);

  auto first_vertex = *first;
  auto current_vertex = *first++;
  for (; first != beyond; ++first) {
    auto last_vertex = current_vertex;
    current_vertex = *first;

    EXPECT_TRUE(lookup_edge(current_vertex, last_vertex, g).second);
  }

  EXPECT_TRUE(lookup_edge(first_vertex, current_vertex, g).second);
}

/// Call the is_bipartite and find_odd_cycle functions and verify their results.

template <typename Graph, typename IndexMap>
void check_bipartite(const Graph& g, IndexMap index_map, bool expect_bipartite) {
  auto partition_map = vector_property_map(num_vertices(g), index_map, default_color_type::white_color);

  std::vector<graph_vertex_descriptor_t<Graph> > odd_cycle;
  odd_cycle.resize(num_vertices(g));

  bool first_result = is_bipartite(g, index_map, partition_map);

  EXPECT_EQ(first_result, is_bipartite(g, index_map));
  EXPECT_EQ(first_result, is_bipartite(g));

  if (first_result) {
    check_two_coloring(g, partition_map);
  }

  EXPECT_EQ(first_result, expect_bipartite);

  auto second_first = odd_cycle.begin();
  auto second_beyond = find_odd_cycle(g, index_map, partition_map, second_first);

  if (expect_bipartite) {
    EXPECT_EQ(second_beyond, second_first);
    check_two_coloring(g, partition_map);
  } else {
    check_odd_cycle(g, second_first, second_beyond);
  }

  second_beyond = find_odd_cycle(g, index_map, second_first);
  if (expect_bipartite) {
    EXPECT_EQ(second_beyond, second_first);
  } else {
    check_odd_cycle(g, second_first, second_beyond);
  }
}

using VectorGraph = adjacency_list<vec_s, vec_s, undirected_s>;
using ListGraph = adjacency_list<list_s, list_s, undirected_s, property<vertex_index_t, std::size_t> >;

template <typename Graph>
Graph make_bipartite_graph() {
  using E = std::pair<std::size_t, std::size_t>;

  /**
   * Create the graph drawn below.
   *
   *       0 - 1 - 2
   *       |       |
   *   3 - 4 - 5 - 6
   *  /      \   /
   *  |        7
   *  |        |
   *  8 - 9 - 10
   **/

  std::array bipartite_edges{E(0, 1), E(0, 4), E(1, 2), E(2, 6),  E(3, 4), E(3, 8), E(4, 5),
                             E(4, 7), E(5, 6), E(6, 7), E(7, 10), E(8, 9), E(9, 10)};

  if constexpr (std::is_same_v<vertex_property_type<Graph>, no_property>) {
    return Graph{11, bipartite_edges};
  } else {
    return Graph{11, std::ranges::iota_view(0, 11), bipartite_edges};
  }
}

template <typename Graph>
Graph make_non_bipartite_graph() {
  using E = std::pair<std::size_t, std::size_t>;

  /**
   * Create the graph drawn below.
   *
   *       2 - 1 - 0
   *       |       |
   *   3 - 6 - 5 - 4
   *  /      \   /
   *  |        7
   *  |       /
   *  8 ---- 9
   *
   **/

  std::array non_bipartite_edges{E(0, 1), E(0, 4), E(1, 2), E(2, 6), E(3, 4), E(3, 8),
                                 E(4, 5), E(4, 7), E(5, 6), E(6, 7), E(7, 9), E(8, 9)};

  if constexpr (std::is_same_v<vertex_property_type<Graph>, no_property>) {
    return Graph{10, non_bipartite_edges};
  } else {
    return Graph{10, std::ranges::iota_view(0, 10), non_bipartite_edges};
  }
}

TEST(BipartiteTest, BipartiteVectorGraph) {
  auto v_graph = make_bipartite_graph<VectorGraph>();
  check_bipartite(v_graph, get(vertex_index, v_graph), true);
}

TEST(BipartiteTest, BipartiteListGraph) {
  auto l_graph = make_bipartite_graph<ListGraph>();
  check_bipartite(l_graph, get(vertex_index, l_graph), true);
}

TEST(BipartiteTest, NonBipartiteVectorGraph) {
  auto v_graph = make_non_bipartite_graph<VectorGraph>();
  check_bipartite(v_graph, get(vertex_index, v_graph), false);
}

TEST(BipartiteTest, NonBipartiteListGraph) {
  auto l_graph = make_non_bipartite_graph<ListGraph>();
  check_bipartite(l_graph, get(vertex_index, l_graph), false);
}

}  // namespace
}  // namespace bagl
