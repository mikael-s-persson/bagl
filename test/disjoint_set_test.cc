// (C) Copyright Jeremy Siek 2002.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <array>

#include "bagl/disjoint_sets.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

template <typename T>
class DisjointSetsTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(DisjointSetsTest);

TYPED_TEST_P(DisjointSetsTest, Test) {
  using DisjointSet = TypeParam;
  // The following tests are pretty lame, just a basic sanity check.
  // Industrial strength tests still need to be written.
  std::array<std::size_t, 4> elts = {0, 1, 2, 3};

  DisjointSet ds(elts.size());

  ds.make_set(elts[0]);
  ds.make_set(elts[1]);
  ds.make_set(elts[2]);
  ds.make_set(elts[3]);

  EXPECT_NE(ds.find_set(0), ds.find_set(1));
  EXPECT_NE(ds.find_set(0), ds.find_set(2));
  EXPECT_NE(ds.find_set(0), ds.find_set(3));
  EXPECT_NE(ds.find_set(1), ds.find_set(2));
  EXPECT_NE(ds.find_set(1), ds.find_set(3));
  EXPECT_NE(ds.find_set(2), ds.find_set(3));

  ds.union_set(0, 1);
  ds.union_set(2, 3);
  EXPECT_NE(ds.find_set(0), ds.find_set(3));
  int a = ds.find_set(0);
  EXPECT_EQ(a, ds.find_set(1));
  int b = ds.find_set(2);
  EXPECT_EQ(b, ds.find_set(3));

  ds.link(a, b);
  EXPECT_EQ(ds.find_set(a), ds.find_set(b));
  EXPECT_EQ(ds.count_sets(elts), 1);

  ds.normalize_sets(elts);
  ds.compress_sets(elts);
  EXPECT_EQ(ds.count_sets(elts), 1);
}

REGISTER_TYPED_TEST_SUITE_P(DisjointSetsTest, Test);

using TestHalving = disjoint_sets_with_storage<identity_property_map, identity_property_map, find_with_path_halving>;
using TestFullComp =
    disjoint_sets_with_storage<identity_property_map, identity_property_map, find_with_full_path_compression>;
using TestSets = ::testing::Types<TestHalving, TestFullComp>;

INSTANTIATE_TYPED_TEST_SUITE_P(DisjointSets, DisjointSetsTest, TestSets);

}  // namespace
}  // namespace bagl
