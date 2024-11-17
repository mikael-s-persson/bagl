// Copyright 2024 Mikael Persson

#include "bagl/vebl_d_ary_tree.h"

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "gtest/gtest.h"
#include "test/graph_properties_test_suite.h"
#include "test/tree_mutation_test_suite.h"

namespace bagl {
namespace {

using VEBLType = vebl_d_ary_tree<4, test_vertex_property, test_edge_property>;
using VEBLPType = std::tuple<vebl_d_ary_tree<4, VertexPropTest, EdgePropTest, GraphPropTest>, PropMapMaps>;
using VEBLBType = std::tuple<vebl_d_ary_tree<4, VertexBundleTest, EdgeBundleTest, GraphBundleTest>, BundleMaps>;

INSTANTIATE_TYPED_TEST_SUITE_P(VEBLTree, TreeMutationTest, VEBLType);
INSTANTIATE_TYPED_TEST_SUITE_P(VEBLTreeProps, GraphPropertiesTest, VEBLPType);
INSTANTIATE_TYPED_TEST_SUITE_P(VEBLTreeBundles, GraphPropertiesTest, VEBLBType);

}  // namespace
}  // namespace bagl
