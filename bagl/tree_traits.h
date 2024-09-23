// Copyright 1999 Jeremy Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TREE_TRAITS_H_
#define BAGL_BAGL_TREE_TRAITS_H_

#include <tuple>

#include "bagl/has_trait_member.h"

namespace bagl {

namespace tree_traits_detail {

BAGL_GRAPH_HAS_TRAIT_MEMBER(node_descriptor, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(children_range, void)

}  // namespace tree_traits_detail

template <typename T>
struct tree_traits {
  using node_descriptor = tree_traits_detail::get_node_descriptor_or_not<G>;
  using children_range = tree_traits_detail::get_children_range_or_not<G>;
};
template <typename T>
using tree_node_descriptor_t = typename tree_traits<T>::node_descriptor;
template <typename T>
using tree_children_range_t = typename tree_traits<T>::children_range;

template <typename Tree, typename TreeVisitor>
void traverse_tree(tree_node_descriptor_t<Tree> v, Tree& t, TreeVisitor visitor) {
  using StackElem =
      std::tuple<tree_node_descriptor_t<Tree>, bool, tree_children_range_t<Tree>, std::ranges::iterator_t<tree_children_range_t<Tree>>>;
  std::vector<StackElem> stack;
  auto vr = children(v, t);
  stack.emplace_back(v, false, vr, vr.begin());
  visitor.preorder(v, t);
  while (!stack.empty()) {
    auto& [u, ud, ur, ui] = stack.back();
    if (ui == ur.end()) {
      if (!ud) {
        visitor.inorder(u, t);
      }
      visitor.postorder(u, t);
      stack.pop();
      if (!stack.empty() && !std::get<1>(stack.back())) {
        visitor.inorder(std::get<0>(stack.back()), t);
        std::get<1>(stack.back()) = true;
      }
      continue;
    }
    auto w = *ui;
    ++ui;
    auto wr = children(w, t);
    stack.emplace_back(w, false, wr, wr.begin());
    visitor.preorder(w, t);
  }
}

struct null_tree_visitor {
  template <typename Node, typename Tree>
  void preorder(Node, Tree&) {}
  template <typename Node, typename Tree>
  void inorder(Node, Tree&) {}
  template <typename Node, typename Tree>
  void postorder(Node, Tree&) {}
};

}  // namespace bagl

#endif  // BAGL_BAGL_TREE_TRAITS_H_
