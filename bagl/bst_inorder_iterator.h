// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_BST_INORDER_ITERATOR_H_
#define BAGL_BAGL_BST_INORDER_ITERATOR_H_

#include <iterator>

#include "bagl/graph_traits.h"
#include "bagl/tree_adaptor.h"
#include "bagl/tree_traits.h"

namespace bagl {

// This library implements a simple binary search tree in-order iterator for a BAGL tree.
// Mostly, this serves as an example usage of trees.

enum class bst_traversal_status { on_left_branch, on_middle_branch, on_right_branch };

namespace bst_detail {

template <typename TreeType, typename VertexType>
VertexType bst_go_down_left(const TreeType& tree, VertexType start) {
  // look down-left until we reach the leaf:
  while (out_degree(start, tree) > 0) {
    start = *children(start, tree).begin();
  }
  return start;
}

template <typename TreeType, typename VertexType>
VertexType bst_go_down_right(const TreeType& tree, VertexType start) {
  // look down-right until we reach the leaf:
  while (out_degree(start, tree) > 1) {
    start = *std::next(children(start, tree).begin());
  }
  return start;
}

template <typename TreeType, typename VertexType>
void bst_move_up_to_next(const TreeType& tree, VertexType& u, bst_traversal_status& status) {
  status = bst_traversal_status::on_right_branch;
  while (true) {
    if (in_degree(u, tree) == 0) {
      // at the root, go to the end (root, rightbranch)
      return;
    }
    VertexType v = parent(u, tree);
    VertexType v_first_child = *children(v, tree).begin();
    if (v_first_child == u) {  // u is the left child of v.
      u = v;
      status = bst_traversal_status::on_middle_branch;
      return;
    }
    // u must be the right child of v. keep going.
    u = v;
  }
}

template <typename TreeType, typename VertexType>
void bst_move_up_to_prev(const TreeType& tree, VertexType& w, bst_traversal_status& status) {
  status = bst_traversal_status::on_left_branch;
  VertexType u = w;
  while (true) {
    if (in_degree(u, tree) == 0) {
      // at the root, so, w must be the beginning node.
      return;
    }
    VertexType v = parent(u, tree);
    VertexType v_first_child = *children(v, tree).begin();
    if (v_first_child == u) {  // u is the left child of v.
      u = v;
      continue;
    }
    // u must be the right child of v. keep going.
    w = v;
    status = bst_traversal_status::on_middle_branch;
    return;
  }
}

}  // namespace bst_detail

template <typename CompleteBinaryTree, typename ValueType>
class bst_inorder_iterator {
 public:
  using self = bst_inorder_iterator<CompleteBinaryTree, ValueType>;
  using tree_type = CompleteBinaryTree;

  using difference_type = std::ptrdiff_t;
  using value_type = ValueType;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::bidirectional_iterator_tag;

 private:
  using vertex_type = typename graph_traits<tree_type>::vertex_descriptor;

  tree_type* tree_;
  vertex_type u_;
  bst_traversal_status status_;

  bst_inorder_iterator(tree_type* tree, vertex_type u, bst_traversal_status status)
      : tree_(tree), u_(u), status_(status) {}

 public:
  vertex_type base() const { return u_; }

  bst_inorder_iterator(tree_type* tree, vertex_type u)
      : tree_(tree), u_(u), status_(bst_traversal_status::on_right_branch) {
    // must figure out what the case is.
    if (tree_ == nullptr || u_ == graph_traits<tree_type>::null_vertex()) {
      return;
    }
    if (u_ == tree_root(*tree_)) {
      status_ = bst_traversal_status::on_middle_branch;
      return;
    }
    // first check if there are any children:
    if (out_degree(u_, *tree_)) {
      status_ = bst_traversal_status::on_middle_branch;  // not on leaf.
      return;
    }
    // then, check if u_ is a left or right child of its parent:
    vertex_type v = parent(u_, *tree_);
    const vertex_type v_first_child = *children(v, *tree_).begin();
    status_ = (v_first_child == u_ ? bst_traversal_status::on_left_branch : bst_traversal_status::on_right_branch);
  }

  static self begin(tree_type* tree) {
    if ((tree) && (num_vertices(*tree))) {
      return self(tree, bst_detail::bst_go_down_left(*tree, tree_root(*tree)), bst_traversal_status::on_left_branch);
    }
    return self(nullptr, graph_traits<tree_type>::null_vertex(), bst_traversal_status::on_right_branch);
  }

  static self end(tree_type* tree) {
    if ((tree) && (num_vertices(*tree))) {
      // Sentinel value.
      return self(tree, tree_root(*tree), bst_traversal_status::on_right_branch);
    }
    return self(nullptr, graph_traits<tree_type>::null_vertex(), bst_traversal_status::on_right_branch);
  }

  friend bool operator==(const self& lhs, const self& rhs) {
    return ((lhs.tree_ == rhs.tree_) && (lhs.u_ == rhs.u_) && (lhs.status_ == rhs.status_));
  }
  friend bool operator!=(const self& lhs, const self& rhs) {
    return ((lhs.tree_ != rhs.tree_) || (lhs.u_ != rhs.u_) || (lhs.status_ != rhs.status_));
  }

  self& operator++() {
    if (!tree_) {
      return *this;
    }
    switch (status_) {
      case bst_traversal_status::on_left_branch:
        // on a left-leaf, must move up to the parent as a middle-point
        {
          if (in_degree(u_, *tree_) == 0) {
            // at the root, go to the end (root, rightbranch).
            status_ = bst_traversal_status::on_right_branch;
            return *this;
          }
          u_ = parent(u_, *tree_);
          status_ = bst_traversal_status::on_middle_branch;
        }
        break;
      case bst_traversal_status::on_middle_branch:
        // on a middle-point, must move down to the right once, and then left to the bottom.
        if (out_degree(u_, *tree_) > 1) {
          // go to the right child.
          const vertex_type v_right_child = *std::next(children(u_, *tree_).begin());
          u_ = bst_detail::bst_go_down_left(*tree_, v_right_child);
          status_ =
              (u_ == v_right_child ? bst_traversal_status::on_right_branch : bst_traversal_status::on_left_branch);
          break;
        }
        [[fallthrough]];
        // this means that we must move up to the next value (no right child here).
      case bst_traversal_status::on_right_branch:
        bst_detail::bst_move_up_to_next(*tree_, u_, status_);
        break;
    }
    return *this;
  }

  self& operator--() {
    if (!tree_) {
      return *this;
    }
    switch (status_) {
      case bst_traversal_status::on_right_branch:
        if (u_ == tree_root(*tree_)) {
          // go to the left child, and down the right:
          u_ = bst_detail::bst_go_down_right(*tree_, u_);
          if (u_ == tree_root(*tree_)) {
            status_ = bst_traversal_status::on_middle_branch;
          } else {
            // On right leaf or middle
            status_ = (out_degree(u_, *tree_) == 0 ? bst_traversal_status::on_right_branch
                                                   : bst_traversal_status::on_middle_branch);
          }
          break;
        }
        [[fallthrough]];
        // this means that we are either on a right-leaf or on a mis-labeled middle-node,
        // in either case, try the middle-branch case:
      case bst_traversal_status::on_middle_branch:
        // on a middle-point or right-point, must move down to the left once (if exists), and then right to the bottom.
        if (out_degree(u_, *tree_) > 0) {
          // go to the left child, and down the right:
          const vertex_type u_first_child = *children(u_, *tree_).begin();
          u_ = bst_detail::bst_go_down_right(*tree_, u_first_child);
          const auto leaf_status =
              (u_ == u_first_child ? bst_traversal_status::on_left_branch : bst_traversal_status::on_right_branch);
          // On leaf or middle
          status_ = (out_degree(u_, *tree_) == 0 ? leaf_status : bst_traversal_status::on_middle_branch);
          break;
        }
        [[fallthrough]];
        // this means that we must move up to the previous value (no left child here).
      case bst_traversal_status::on_left_branch:
        bst_detail::bst_move_up_to_prev(*tree_, u_, status_);
        break;
    }
    return *this;
  }

  self operator++(int) {
    self result(*this);
    return ++result;
  }
  self operator--(int) {
    self result(*this);
    return --result;
  }

  reference operator*() const { return (*tree_)[u_]; };
  pointer operator->() const { return &(*tree_)[u_]; };
};

}  // namespace bagl

#endif // BAGL_BAGL_BST_INORDER_ITERATOR_H_
