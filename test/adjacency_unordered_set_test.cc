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

  using SVL = adjacency_list<unordered_set_s, vec_s, Kind, vprop, eprop>;
  using SLL = adjacency_list<unordered_set_s, list_s, Kind, vprop, eprop>;
  using SPL = adjacency_list<unordered_set_s, pool_s, Kind, vprop, eprop>;

  using allgraphs = ::testing::Types<SVL, SLL, SPL>;
};

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMulti, GraphMutationTest, test_graphs<undirected_s>::allgraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMulti, GraphMutationTest, test_graphs<directed_s>::allgraphs);

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMulti, GraphMutationTest, test_graphs<bidirectional_s>::allgraphs);

}  // namespace
}  // namespace bagl
