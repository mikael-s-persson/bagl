// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ARRAY_BINARY_TREE_H_
#define BAGL_BAGL_ARRAY_BINARY_TREE_H_

#include <concepts>
#include <functional>
#include <iterator>
#include <ranges>

namespace bagl {

// array_binary_tree is a very simple, completely balanced binary tree.
// For a more general and useful array trees, see bfl_d_ary_tree.h or vebl_d_ary_tree.h.
template <class RAIter, class ID>
class array_binary_tree_node {
 public:
  using self = array_binary_tree_node<RAIter, ID>;
  using rep_iterator = RAIter;
  using difference_type = std::iter_difference_t<RAIter>;
  using value_type = std::iter_value_t<RAIter>;
  using size_type = difference_type;

  array_binary_tree_node() = default;
  array_binary_tree_node(rep_iterator start, rep_iterator end, rep_iterator pos, const ID& _id)
      : r(start), i(pos - start), n(end - start), id(_id) {}
  array_binary_tree_node(rep_iterator rr, size_type ii, size_type nn, const ID& _id) : r(rr), i(ii), n(nn), id(_id) {}
  [[nodiscard]] value_type& value() { return *(r + i); }
  [[nodiscard]] const value_type& value() const { return *(r + i); }
  auto parent() const { return self{r, (i - 1) / 2, n, id}; }
  [[nodiscard]] bool has_parent() const { return i != 0; }
  [[nodiscard]] size_type num_children() const { 
    size_type c = 2 * i + 1;
    size_type s = 0;
    if (c + 1 < n) {
      s = 2;
    } else if (c < n) {
      s = 1;
    }
    return s;
  }
  [[nodiscard]] auto children() const {
    return std::views::iota(2 * i + 1, 2 * i + 1 + num_children()) | std::views::transform([this](size_type j) {
             return self{r, j, n, id};
           });
  }

  template <class ExternalData>
  void swap(self x, ExternalData& edata) {
    value_type tmp = std::move(x.value());

    /*swap external data*/
    edata[get(id, tmp)] = i;
    edata[get(id, value())] = x.i;

    x.value() = std::move(value());
    value() = std::move(tmp);
    i = x.i;
  }
  [[nodiscard]] size_type index() const { return i; }

  rep_iterator r;
  size_type i = 0;
  size_type n = 0;
  ID id;
};

template <class RAContainer, class Compare = std::less<typename RAContainer::value_type> >
struct compare_array_node {
  using value_type = typename RAContainer::value_type;
  explicit compare_array_node(const Compare& x) : comp(x) {}

  template <class node_type>
  [[nodiscard]] bool operator()(const node_type& x, const node_type& y) {
    return comp(x.value(), y.value());
  }

  template <class node_type>
  [[nodiscard]] bool operator()(const node_type& x, const node_type& y) const {
    return comp(x.value(), y.value());
  }
  Compare comp;
};

}  // namespace bagl

#endif // BAGL_BAGL_ARRAY_BINARY_TREE_H_
