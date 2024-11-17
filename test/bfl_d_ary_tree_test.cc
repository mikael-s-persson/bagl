// Copyright 2024 Mikael Persson

#include "bagl/bfl_d_ary_tree.h"

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "gtest/gtest.h"
#include "test/graph_properties_test_suite.h"
#include "test/tree_mutation_test_suite.h"

namespace bagl {
namespace {

using BFLType = bfl_d_ary_tree<4, test_vertex_property, test_edge_property>;
using BFLPType = std::tuple<bfl_d_ary_tree<4, VertexPropTest, EdgePropTest, GraphPropTest>, PropMapMaps>;
using BFLBType = std::tuple<bfl_d_ary_tree<4, VertexBundleTest, EdgeBundleTest, GraphBundleTest>, BundleMaps>;

INSTANTIATE_TYPED_TEST_SUITE_P(BFLTree, TreeMutationTest, BFLType);
INSTANTIATE_TYPED_TEST_SUITE_P(BFLTreeProps, GraphPropertiesTest, BFLPType);
INSTANTIATE_TYPED_TEST_SUITE_P(BFLTreeBundles, GraphPropertiesTest, BFLBType);

}  // namespace
}  // namespace bagl
