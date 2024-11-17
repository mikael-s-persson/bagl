// Copyright 2024 Mikael Persson

#include "bagl/graph_concepts.h"
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/linked_tree.h"
#include "gtest/gtest.h"
#include "test/graph_properties_test_suite.h"
#include "test/tree_mutation_test_suite.h"

namespace bagl {
namespace {

template <typename Kind, typename VertexProp, typename EdgeProp, typename GraphProp = no_property>
struct test_graphs {
  using LVL = linked_tree<list_s, vec_s, Kind, VertexProp, EdgeProp, GraphProp>;
  using LLL = linked_tree<list_s, list_s, Kind, VertexProp, EdgeProp, GraphProp>;
  using LPL = linked_tree<list_s, pool_s, Kind, VertexProp, EdgeProp, GraphProp>;

  using allgraphs = ::testing::Types<LVL, LLL, LPL>;
  using allgraphs_prop =
      ::testing::Types<std::tuple<LVL, PropMapMaps>, std::tuple<LLL, PropMapMaps>, std::tuple<LPL, PropMapMaps>>;
  using allgraphs_bund =
      ::testing::Types<std::tuple<LVL, BundleMaps>, std::tuple<LLL, BundleMaps>, std::tuple<LPL, BundleMaps>>;
};

using BidirTypes = test_graphs<bidirectional_s, test_vertex_property, test_edge_property>::allgraphs;
using BidirPTypes = test_graphs<bidirectional_s, VertexPropTest, EdgePropTest, GraphPropTest>::allgraphs_prop;
using BidirBTypes = test_graphs<bidirectional_s, VertexBundleTest, EdgeBundleTest, GraphBundleTest>::allgraphs_bund;

INSTANTIATE_TYPED_TEST_SUITE_P(LinkedTreeBidirMulti, TreeMutationTest, BidirTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(LinkedTreeBidirMultiProps, GraphPropertiesTest, BidirPTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(LinkedTreeBidirMultiBundles, GraphPropertiesTest, BidirBTypes);

using DirTypes = test_graphs<directed_s, test_vertex_property, test_edge_property>::allgraphs;
using DirPTypes = test_graphs<directed_s, VertexPropTest, EdgePropTest, GraphPropTest>::allgraphs_prop;
using DirBTypes = test_graphs<directed_s, VertexBundleTest, EdgeBundleTest, GraphBundleTest>::allgraphs_bund;

INSTANTIATE_TYPED_TEST_SUITE_P(LinkedTreeDirMulti, TreeMutationTest, DirTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(LinkedTreeDirMultiProps, GraphPropertiesTest, DirPTypes);
INSTANTIATE_TYPED_TEST_SUITE_P(LinkedTreeDirMultiBundles, GraphPropertiesTest, DirBTypes);

}  // namespace
}  // namespace bagl
