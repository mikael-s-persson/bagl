//  (C) Copyright Jeremy Siek 2004
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_DISJOINT_SETS_H_
#define BAGL_BAGL_DETAIL_DISJOINT_SETS_H_

#include <cassert>

#include "bagl/property_map.h"

namespace bagl::disjoint_sets_detail {

template <typename ParentPA>
auto find_representative_with_path_halving(ParentPA p, property_traits_value_t<ParentPA> v) {
  auto parent = get(p, v);
  auto grandparent = get(p, parent);
  while (parent != grandparent) {
    put(p, v, grandparent);
    v = grandparent;
    parent = get(p, v);
    grandparent = get(p, parent);
  }
  return parent;
}

template <typename ParentPA>
auto find_representative_with_full_compression(ParentPA parent, property_traits_value_t<ParentPA> v) {
  auto old = v;
  auto ancestor = get(parent, v);
  while (ancestor != v) {
    v = ancestor;
    ancestor = get(parent, v);
  }
  v = get(parent, old);
  while (ancestor != v) {
    put(parent, old, ancestor);
    old = v;
    v = get(parent, old);
  }
  return ancestor;
}

// the postcondition of link sets is:
// component_representative(i) == component_representative(j)
template <typename ParentPA, typename RankPA>
void link_sets(ParentPA p, RankPA rank, property_traits_value_t<ParentPA> i, property_traits_value_t<ParentPA> j) {
  assert(i == get(p, i));
  assert(j == get(p, j));
  if (i == j) {
    return;
  }
  if (get(rank, i) > get(rank, j)) {
    put(p, j, i);
  } else {
    put(p, i, j);
    if (get(rank, i) == get(rank, j)) {
      put(rank, j, get(rank, j) + 1);
    }
  }
}

// normalize components has the following postcondidition:
// i >= p[i]
// that is, the representative is the node with the smallest index in its
// class as its precondition it it assumes that the node container is
// compressed

template <typename ParentPA>
void normalize_node(ParentPA p, property_traits_value_t<ParentPA> i) {
  if (i > get(p, i) || get(p, get(p, i)) != get(p, i)) {
    put(p, i, get(p, get(p, i)));
  } else {
    put(p, get(p, i), i);
    put(p, i, i);
  }
}

}  // namespace bagl::disjoint_sets_detail

#endif  // BAGL_BAGL_DETAIL_DISJOINT_SETS_H_
