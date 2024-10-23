// Copyright (C) Vladimir Prus 2003.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/copy.h"

#include "bagl/adjacency_list.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

class copier {
 public:
  template <class V1, class V2>
  void operator()(const V1&, const V2&) const {}
};

TEST(CopyTest, AdjList) {
  adjacency_list<vec_s, vec_s, directed_s, property<vertex_root_t, int> > g1;
  adjacency_list<vec_s, vec_s, directed_s, property<vertex_root_t, int> > g2;
  adjacency_list<multiset_s, vec_s, directed_s, property<vertex_index_t, int> > g3;

  copy_graph(g1, g2);
  copier c;
  copy_graph(g3, g1, c, default_edge_copier(g3, g1));
}

}  // namespace
}  // namespace bagl
