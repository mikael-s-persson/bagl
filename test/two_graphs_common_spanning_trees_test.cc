// Copyright (C) 2012, Michele Caini.
// Copyright 2024 Mikael Persson - Modernized to C++20

// Two Graphs Common Spanning Trees Algorithm
// Based on academic article of Mint, Read and Tarjan
// Efficient Algorithm for Common Spanning Tree Problem
// Electron. Lett., 28 April 1983, Volume 19, Issue 9, p.346-347

#include "bagl/two_graphs_common_spanning_trees.h"

#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/property.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

using Graph = adjacency_list<vec_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t>>;

template <typename Coll, typename Seq>
void edge_checker(Coll& coll, Seq& seq) {
  bool found = false;
  for (auto& coll_seq : coll) {
    EXPECT_EQ(coll_seq.size(), seq.size());
    if (std::ranges::equal(coll_seq, seq)) {
      found = true;
      break;
    }
  }
  EXPECT_TRUE(found);
}

TEST(TwoGraphsCommonSpanningTreesTest, Basic) {
  std::vector ig_edges = {
      std::pair{0, 1}, std::pair{1, 3}, std::pair{3, 2}, std::pair{1, 5}, std::pair{5, 4},
      std::pair{5, 6}, std::pair{5, 3}, std::pair{3, 1}, std::pair{1, 3},
  };
  Graph ig(7, ig_edges, std::ranges::iota_view(0));

  std::vector vg_edges = {
      std::pair{0, 2}, std::pair{0, 4}, std::pair{0, 5}, std::pair{5, 1}, std::pair{5, 3},
      std::pair{5, 6}, std::pair{5, 4}, std::pair{5, 2}, std::pair{2, 6},
  };
  Graph vg(7, vg_edges, std::ranges::iota_view(0));

  std::vector<std::vector<bool>> coll;
  tree_collector<std::vector<std::vector<bool>>, std::vector<bool>> collector(coll);

  two_graphs_common_spanning_trees(ig, get(vertex_index, ig), get(edge_index, ig), vg, get(vertex_index, vg),
                                   get(edge_index, vg), collector);

  std::vector<bool> check;
  check = {true, true, true, true, true, true, false, false, false};
  edge_checker(coll, check);

  check = {true, true, true, true, true, true, false, false, false};
  edge_checker(coll, check);
}

}  // namespace
}  // namespace bagl
