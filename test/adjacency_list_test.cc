// (C) Copyright Andrew Sutton 2009
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"

#include <cstddef>
#include <iostream>
#include <random>
#include <string_view>

#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/random.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include "test/graph_mutation_test_suite.h"

namespace bagl {
namespace {

template <typename Kind>
struct test_graphs {
  using vprop = test_vertex_property;
  using eprop = test_edge_property;

  using VVL = adjacency_list<vec_s, vec_s, Kind, vprop, eprop>;
  using PVL = adjacency_list<pool_s, vec_s, Kind, vprop, eprop>;
  using LVL = adjacency_list<list_s, vec_s, Kind, vprop, eprop>;
  using MVL = adjacency_list<multiset_s, vec_s, Kind, vprop, eprop>;
  using UMVL = adjacency_list<unordered_multiset_s, vec_s, Kind, vprop, eprop>;

  using VLL = adjacency_list<vec_s, list_s, Kind, vprop, eprop>;
  using PLL = adjacency_list<pool_s, list_s, Kind, vprop, eprop>;
  using LLL = adjacency_list<list_s, list_s, Kind, vprop, eprop>;
  using MLL = adjacency_list<multiset_s, list_s, Kind, vprop, eprop>;
  using UMLL = adjacency_list<unordered_multiset_s, list_s, Kind, vprop, eprop>;

  using VPL = adjacency_list<vec_s, pool_s, Kind, vprop, eprop>;
  using PPL = adjacency_list<pool_s, pool_s, Kind, vprop, eprop>;
  using LPL = adjacency_list<list_s, pool_s, Kind, vprop, eprop>;
  using MPL = adjacency_list<multiset_s, pool_s, Kind, vprop, eprop>;
  using UMPL = adjacency_list<unordered_multiset_s, pool_s, Kind, vprop, eprop>;

  using SVL = adjacency_list<set_s, vec_s, Kind, vprop, eprop>;
  using USVL = adjacency_list<unordered_set_s, vec_s, Kind, vprop, eprop>;

  using SLL = adjacency_list<set_s, list_s, Kind, vprop, eprop>;
  using USLL = adjacency_list<unordered_set_s, list_s, Kind, vprop, eprop>;

  using SPL = adjacency_list<set_s, pool_s, Kind, vprop, eprop>;
  using USPL = adjacency_list<unordered_set_s, pool_s, Kind, vprop, eprop>;

  using allgraphs = ::testing::Types<VVL, PVL, LVL, MVL, UMVL, VLL, PLL, LLL, MLL, UMLL, VPL, PPL, LPL, MPL, UMPL, SVL,
                                     USVL, SLL, USLL, SPL, USPL>;
};

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMulti, GraphMutationTest, test_graphs<undirected_s>::allgraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMulti, GraphMutationTest, test_graphs<directed_s>::allgraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMulti, GraphMutationTest, test_graphs<bidirectional_s>::allgraphs);

}  // namespace
}  // namespace bagl
