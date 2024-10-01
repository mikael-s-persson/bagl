// Copyright 2002 Indiana University.
// Copyright 2009 Trustees of Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek, Michael Hansen
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_INCREMENTAL_COMPONENTS_H_
#define BAGL_BAGL_DETAIL_INCREMENTAL_COMPONENTS_H_

#include <iterator>

namespace bagl::inc_comp_detail {

// Iterator for a component index linked list.  The contents of
// each array element represent the next index in the list.  A
// special value (the maximum index + 1) is used to terminate a
// list.
template <typename IndexIter>
class component_index_iterator {
 private:
  using self = component_index_iterator<IndexIter>;

 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::iter_value_t<IndexIter>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::iter_difference_t<IndexIter>;

  // Constructor for "begin" iterator
  component_index_iterator(IndexIter index_iterator, value_type begin_index)
      : index_iterator_(index_iterator), current_index_(begin_index) {}

  // Constructor for "end" iterator (end_index should be the linked
  // list terminator).
  explicit component_index_iterator(value_type end_index = 0) : current_index_(end_index) {}

  reference operator*() const { return current_index_; }
  pointer operator->() const { return &current_index_; }

  self& operator++() {
    // Move to the next element in the linked list
    current_index_ = index_iterator_[current_index_];
    return *this;
  }

  self operator++(int) {
    self tmp{*this};
    // Move to the next element in the linked list
    current_index_ = index_iterator_[current_index_];
    return tmp;
  }

  bool operator==(const self& other_iterator) const { return current_index_ == other_iterator.current_index_; }
  bool operator!=(const self& other_iterator) const { return current_index_ != other_iterator.current_index_; }

 protected:
  IndexIter index_iterator_;
  value_type current_index_ = 0;
};

}  // namespace bagl::inc_comp_detail

#endif  // BAGL_BAGL_DETAIL_INCREMENTAL_COMPONENTS_H_
