// Copyright 2024 Mikael Persson

#include "bagl/adjacency_list.h"
#include "gtest/gtest.h"
#include "test/graph_mutation_test_suite.h"

namespace bagl {
namespace {

template <typename Kind>
struct test_graphs {
  using vprop = test_vertex_property;
  using eprop = test_edge_property;

  using VVL = adjacency_list<vec_s, vec_s, Kind, vprop, eprop>;
  using VLL = adjacency_list<vec_s, list_s, Kind, vprop, eprop>;
  using VPL = adjacency_list<vec_s, pool_s, Kind, vprop, eprop>;

  using allgraphs = ::testing::Types<VVL, VLL, VPL>;
};

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMulti, GraphMutationTest, test_graphs<undirected_s>::allgraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMulti, GraphMutationTest, test_graphs<directed_s>::allgraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMulti, GraphMutationTest, test_graphs<bidirectional_s>::allgraphs);

}  // namespace
}  // namespace bagl
