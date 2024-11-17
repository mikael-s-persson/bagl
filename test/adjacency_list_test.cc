// Copyright 2024 Mikael Persson

#include "bagl/adjacency_list.h"

#include "gtest/gtest.h"
#include "test/graph_mutation_test_suite.h"
#include "test/graph_properties_test_suite.h"

namespace bagl {
namespace {

template <typename Kind, typename VertexProp, typename EdgeProp, typename GraphProp = no_property>
struct test_graphs {
  using LVL = adjacency_list<list_s, vec_s, Kind, VertexProp, EdgeProp, GraphProp>;
  using LLL = adjacency_list<list_s, list_s, Kind, VertexProp, EdgeProp, GraphProp>;
  using LPL = adjacency_list<list_s, pool_s, Kind, VertexProp, EdgeProp, GraphProp>;

  using allgraphs = ::testing::Types<LVL, LLL, LPL>;
  using allgraphs_prop =
      ::testing::Types<std::tuple<LVL, PropMapMaps>, std::tuple<LLL, PropMapMaps>, std::tuple<LPL, PropMapMaps>>;
  using allgraphs_bund =
      ::testing::Types<std::tuple<LVL, BundleMaps>, std::tuple<LLL, BundleMaps>, std::tuple<LPL, BundleMaps>>;
};

using UndirGTypes = test_graphs<undirected_s, test_vertex_property, test_edge_property>::allgraphs;
using UndirPTypes = test_graphs<undirected_s, VertexPropTest, EdgePropTest, GraphPropTest>::allgraphs_prop;
using UndirBTypes = test_graphs<undirected_s, VertexBundleTest, EdgeBundleTest, GraphBundleTest>::allgraphs_bund;

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMulti, GraphMutationTest, UndirGTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMultiTree, GraphTreeMutationTest, UndirGTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMultiProps, GraphPropertiesTest, UndirPTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListUndirMultiBundles, GraphPropertiesTest, UndirBTypes);

using DirGTypes = test_graphs<directed_s, test_vertex_property, test_edge_property>::allgraphs;
using DirPTypes = test_graphs<directed_s, VertexPropTest, EdgePropTest, GraphPropTest>::allgraphs_prop;
using DirBTypes = test_graphs<directed_s, VertexBundleTest, EdgeBundleTest, GraphBundleTest>::allgraphs_bund;

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMulti, GraphMutationTest, DirGTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMultiTree, GraphTreeMutationTest, DirGTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMultiProps, GraphPropertiesTest, DirPTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListDirMultiBundles, GraphPropertiesTest, DirBTypes);

using BidirGTypes = test_graphs<bidirectional_s, test_vertex_property, test_edge_property>::allgraphs;
using BidirPTypes = test_graphs<bidirectional_s, VertexPropTest, EdgePropTest, GraphPropTest>::allgraphs_prop;
using BidirBTypes = test_graphs<bidirectional_s, VertexBundleTest, EdgeBundleTest, GraphBundleTest>::allgraphs_bund;

INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMulti, GraphMutationTest, BidirGTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMultiTree, GraphTreeMutationTest, BidirGTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMultiProps, GraphPropertiesTest, BidirPTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBidirMultiBundles, GraphPropertiesTest, BidirBTypes);

}  // namespace
}  // namespace bagl
