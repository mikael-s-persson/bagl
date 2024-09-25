// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MUTABLE_QUEUE_H_
#define BAGL_BAGL_MUTABLE_QUEUE_H_

#include <algorithm>
#include <functional>
#include <iterator>
#include <ranges>
#include <vector>

#include "bagl/detail/array_binary_tree.h"
#include "bagl/is_heap.h"
#include "bagl/mutable_heap.h"
#include "bagl/property_map.h"

namespace bagl {

// The mutable queue whose elements are indexed
//
// This adaptor provides a special kind of priority queue that has
// and update operation. This allows the ordering of the items to
// change. After the ordering criteria for item x changes, one must
// call the Q.update(x)
//
// In order to efficiently find x in the queue, a functor must be
// provided to map value_type to a unique ID, which the
// mutable_queue will then use to map to the location of the
// item. The ID's generated must be between 0 and N, where N is the
// value passed to the constructor of mutable_queue

template <class IndexedType, class RAContainer = std::vector<IndexedType>,
          class Comp = std::less<typename RAContainer::value_type>, class ID = identity_property_map>
class mutable_queue {
 public:
  using value_type = IndexedType;
  using size_type = typename RAContainer::size_type;

 protected:
  using iterator = typename RAContainer::iterator;
  using Node = array_binary_tree_node<iterator, ID>;
  using Compare = compare_array_node<RAContainer, Comp>;
  using IndexArray = std::vector<size_type>;

 public:
  using value_compare = Compare;
  using id_generator = ID;

  mutable_queue(size_type n, const Comp& x, const ID& id) : index_array_(n), comp_(x), id_(id) { c_.reserve(n); }
  template <class ForwardIterator>
  mutable_queue(ForwardIterator first, ForwardIterator last, const Comp& x, const ID& id)
      : index_array_(std::distance(first, last)), comp_(x), id_(id) {
    while (first != last) {
      push(*first);
      ++first;
    }
  }

  [[nodiscard]] bool empty() const { return c_.empty(); }

  void pop() {
    using std::swap;
    swap(c_.front(), c_.back());

    size_type id_f = get(id_, c_.back());
    size_type id_b = get(id_, c_.front());
    swap(index_array_[id_f], index_array_[id_b]);

    c_.pop_back();
    Node node(c_.begin(), c_.end(), c_.begin(), id_);
    down_heap(node, comp_, index_array_);
  }
  void push(const IndexedType& x) {
    c_.push_back(x);
    /*set index-array*/
    index_array_[get(id_, x)] = c_.size() - 1;
    Node node(c_.begin(), c_.end(), c_.end() - 1, id_);
    up_heap(node, comp_, index_array_);
  }

  void update(const IndexedType& x) {
    size_type current_pos = index_array_[get(id_, x)];
    c_[current_pos] = x;

    Node node(c_.begin(), c_.end(), c_.begin() + current_pos, id_);
    update_heap(node, comp_, index_array_);
  }

  [[nodiscard]] value_type& front() { return c_.front(); }
  [[nodiscard]] value_type& top() { return c_.front(); }

  [[nodiscard]] const value_type& front() const { return c_.front(); }
  [[nodiscard]] const value_type& top() const { return c_.front(); }

  [[nodiscard]] size_type size() const { return c_.size(); }

  void clear() { c_.clear(); }

  [[nodiscard]] bool check_heap_property() const {
    return is_heap(c_.begin(), c_.end(), Comp());
  }

 protected:
  IndexArray index_array_;
  Compare comp_;
  RAContainer c_;
  ID id_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_MUTABLE_QUEUE_H_
