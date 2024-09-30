// Copyright 1999 Jeremy Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TREE_TRAITS_H_
#define BAGL_BAGL_TREE_TRAITS_H_

#include <tuple>
#include <vector>

#include "bagl/graph_traits.h"
#include "bagl/has_trait_member.h"
#include "bagl/partial_range.h"

namespace bagl {

namespace tree_traits_detail {

BAGL_GRAPH_HAS_TRAIT_MEMBER(node_descriptor, void)

}  // namespace tree_traits_detail

template <typename T>
struct tree_traits {
  using node_descriptor = tree_traits_detail::get_node_descriptor_or_not<T>;
};
template <typename T>
using tree_node_descriptor_t = typename tree_traits<T>::node_descriptor;

template <typename Tree, typename TreeVisitor>
void traverse_tree(tree_node_descriptor_t<Tree> v, Tree& t, TreeVisitor visitor) {
  using StackElem = std::tuple<tree_node_descriptor_t<Tree>, bool, partial_view<decltype(children(v, t))>>;
  std::vector<StackElem> stack;
  stack.emplace_back(v, false, partial_view(children(v, t)));
  visitor.preorder(v, t);
  while (!stack.empty()) {
    auto& [u, ud, ur] = stack.back();
    if (ur.empty()) {
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
    auto ui = ur.begin();
    auto w = *ui;
    ur.move_begin_to(++ui);
    stack.emplace_back(w, false, partial_view(children(w, t)));
    visitor.preorder(w, t);
  }
}

struct null_tree_visitor {
  template <typename Node, typename Tree>
  void preorder(Node /*unused*/, Tree& /*unused*/) {}
  template <typename Node, typename Tree>
  void inorder(Node /*unused*/, Tree& /*unused*/) {}
  template <typename Node, typename Tree>
  void postorder(Node /*unused*/, Tree& /*unused*/) {}
};

template <typename VertexDescriptor, typename EdgeDescriptor,
          typename StorageTag>
struct tree_storage {
  using type = typename StorageTag::template bind<VertexDescriptor,
                                                  EdgeDescriptor>::type;
};

template <typename StorageTag>
struct tree_storage_traits {
  static constexpr bool is_rand_access_v = false;
  static constexpr bool is_bidir_v = true;
  static constexpr bool is_directed_v = true;

  using directed_category = std::conditional_t<
      is_bidir_v, bidirectional_tag,
      std::conditional_t<is_directed_v, directed_tag, undirected_tag>>;

  using edge_parallel_category = disallow_parallel_edge_tag;

  using vertices_size_type = std::size_t;
  using vertex_descriptor = typename StorageTag::vertex_descriptor;
  using edges_size_type = std::size_t;
  using edge_descriptor = typename StorageTag::edge_descriptor;
};

}  // namespace bagl

#endif  // BAGL_BAGL_TREE_TRAITS_H_
