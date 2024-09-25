// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MUTABLE_HEAP_H_
#define BAGL_BAGL_MUTABLE_HEAP_H_

// There are a few things wrong with this set of functions.
//
// ExternalData should be removed, it is not part of the core
// algorithm. It can be handled inside the tree nodes.
//
// The min_element should be replaced by a fixed length loop
// (fixed at d for d-heaps).
//
// The member functions of TreeNode should be changed to global
// functions.
//
// These functions will be replaced by those in heap_tree.h

#include <algorithm>

namespace bagl {

template <class TreeNode, class Compare, class ExternalData>
TreeNode up_heap(TreeNode x, const Compare& comp, ExternalData& edata) {
  while (x.has_parent() && comp(x, x.parent())) {
    x.swap(x.parent(), edata);
  }
  return x;
}

template <class TreeNode, class Compare, class ExternalData>
TreeNode down_heap(TreeNode x, const Compare& comp, ExternalData& edata) {
  while (x.num_children() > 0) {
    auto children = x.children();
    auto child_it = std::min_element(children.begin(), children.end(), comp);
    if (!comp(*child_it, x)) {
      break;
    }
    x.swap(*child_it, edata);
  }
  return x;
}

template <class TreeNode, class Compare, class ExternalData>
void update_heap(TreeNode x, const Compare& comp, ExternalData& edata) {
  x = down_heap(x, comp, edata);
  (void)up_heap(x, comp, edata);
}

}  // namespace bagl

#endif  // BAGL_BAGL_MUTABLE_HEAP_H_
