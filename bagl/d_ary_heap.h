// Copyright 2009 Trustees of Indiana University
// Authors: Jeremiah J. Willcock, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_D_ARY_HEAP_H_
#define BAGL_BAGL_D_ARY_HEAP_H_

#include <algorithm>
#include <cstddef>
#include <limits>
#include <span>
#include <utility>
#include <vector>

#include "bagl/property_map.h"

// WARNING: it is not safe to copy a d_ary_heap_indirect and then modify one of
// the copies.  The class is required to be copyable so it can be passed around
// (without move support from C++11), but it deep-copies the heap contents yet
// shallow-copies the index_in_heap_map.

namespace bagl {

// D-ary heap using an indirect compare operator (use identity_property_map
// as DistanceMap to get a direct compare operator).  This heap appears to be
// commonly used for Dijkstra's algorithm for its good practical performance
// on some platforms; asymptotically, it has an O(lg N) decrease-key
// operation while that can be done in constant time on a relaxed heap.  The
// implementation is mostly based on the binary heap page on Wikipedia and
// online sources that state that the operations are the same for d-ary
// heaps.  This code is not based on the old Boost d-ary heap code.
//
// - d_ary_heap_indirect is a model of UpdatableQueue as is needed for
//   dijkstra_shortest_paths.
//
// - Value must model Assignable.
// - Arity must be at least 2 (optimal value appears to be 4, both in my and
//   third-party experiments).
// - IndexInHeapMap must be a ReadWritePropertyMap from Value to
//   std::size_t (to store the index of each stored value within the
//   heap for decrease-key aka update).
// - DistanceMap must be a ReadablePropertyMap from Value to something
//   (typedef'ed as distance_type).
// - Compare must be a BinaryPredicate used as a less-than operator on
//   distance_type.
//
template <std::semiregular Value, std::size_t Arity, concepts::ReadWritePropertyMap<Value> IndexInHeapPropertyMap,
          concepts::ReadablePropertyMap<Value> DistanceMap,
          typename Compare = std::less<> >
class d_ary_heap_indirect {
  static_assert(Arity >= 2);

 public:
  using container_type = std::vector<Value>;
  using size_type = typename container_type::size_type;
  using value_type = Value;
  using key_type = property_traits_value_t<DistanceMap>;
  using key_map = DistanceMap;

  static size_type invalid_index() { return std::numeric_limits<size_type>::max(); }

  d_ary_heap_indirect(DistanceMap distance, IndexInHeapPropertyMap index_in_heap, Compare compare = Compare())
      : compare_(std::move(compare)), distance_(std::move(distance)), index_in_heap_(std::move(index_in_heap)) {}

  d_ary_heap_indirect(const d_ary_heap_indirect&) = delete;
  d_ary_heap_indirect& operator=(const d_ary_heap_indirect&) = delete;
  d_ary_heap_indirect(d_ary_heap_indirect&&) noexcept = default;
  d_ary_heap_indirect& operator=(d_ary_heap_indirect&&) noexcept = default;
  ~d_ary_heap_indirect() = default;

  [[nodiscard]] size_type size() const { return data_.size(); }

  [[nodiscard]] bool empty() const { return data_.empty(); }

  template <typename... Args>
  void emplace(Args&&... args) {
    size_type index = data_.size();
    data_.emplace_back(std::forward<Args>(args)...);
    put(index_in_heap_, data_.back(), index);
    preserve_heap_property_up(index);
  }

  void push(const Value& v) {
    emplace(v);
  }

  [[nodiscard]] Value& top() {
    assert(!this->empty());
    return data_[0];
  }

  [[nodiscard]] const Value& top() const {
    assert(!this->empty());
    return data_[0];
  }

  void pop() {
    assert(!this->empty());
    put(index_in_heap_, data_[0], invalid_index());
    if (data_.size() != 1) {
      data_[0] = data_.back();
      put(index_in_heap_, data_[0], 0);
      data_.pop_back();
      preserve_heap_property_down();
    } else {
      data_.pop_back();
    }
  }

  // This function assumes the key has been updated (using an external write
  // to the distance map or such)
  // See
  // http://coding.derkeiler.com/Archive/General/comp.theory/2007-05/msg00043.html
  void update(const Value& v) { /* decrease-key */
    size_type index = get(index_in_heap_, v);
    preserve_heap_property_up(index);
  }

  bool contains(const Value& v) const {
    size_type index = get(index_in_heap_, v);
    return (index != invalid_index());
  }

  template <typename... Args>
  void emplace_or_update(Args&&... args) { /* insert if not present, else update */
    Value v{std::forward<Args>(args)...};
    size_type index = get(index_in_heap_, v);
    if (index == invalid_index()) {
      index = data_.size();
      data_.emplace_back(std::move(v));
      put(index_in_heap_, data_.back(), index);
    }
    preserve_heap_property_up(index);
  }

  void push_or_update(const Value& v) {
    emplace_or_update(v);
  }

  [[nodiscard]] DistanceMap keys() const { return distance_; }

 private:
  Compare compare_;
  container_type data_;
  DistanceMap distance_;
  IndexInHeapPropertyMap index_in_heap_;

  // The distances being compared using compare and that are stored in the
  // distance map
  using distance_type = property_traits_value_t<DistanceMap>;

  // Get the parent of a given node in the heap
  static size_type parent(size_type index) { return (index - 1) / Arity; }

  // Get the first child of a given node
  static size_type first_child(size_type index) { return index * Arity + 1; }

  // Swap two elements in the heap by index, updating index_in_heap
  void swap_heap_elements(size_type index_a, size_type index_b) {
    using std::swap;
    Value value_a = data_[index_a];
    Value value_b = data_[index_b];
    data_[index_a] = value_b;
    data_[index_b] = value_a;
    put(index_in_heap_, value_a, index_b);
    put(index_in_heap_, value_b, index_a);
  }

  // Emulate the indirect_cmp that is now folded into this heap class
  bool compare_indirect(const Value& a, const Value& b) const { return compare_(get(distance_, a), get(distance_, b)); }

  // Starting at a node, move up the tree swapping elements to preserve the
  // heap property
  void preserve_heap_property_up(size_type index) {
    size_type orig_index = index;
    size_type num_levels_moved = 0;
    // The first loop just saves swaps that need to be done in order to
    // avoid aliasing issues in its search; there is a second loop that does
    // the necessary swap operations
    if (index == 0) {
      return;  // Do nothing on root
    }
    Value currently_being_moved = data_[index];
    distance_type currently_being_moved_dist = get(distance_, currently_being_moved);
    for (;;) {
      if (index == 0) {
        break;  // Stop at root
      }
      size_type parent_index = parent(index);
      Value parent_value = data_[parent_index];
      if (compare_(currently_being_moved_dist, get(distance_, parent_value))) {
        ++num_levels_moved;
        index = parent_index;
        continue;
      }
      break;  // Heap property satisfied
    }
    // Actually do the moves -- move num_levels_moved elements down in the
    // tree, then put currently_being_moved at the top
    index = orig_index;
    for (size_type i = 0; i < num_levels_moved; ++i) {
      size_type parent_index = parent(index);
      Value parent_value = data_[parent_index];
      put(index_in_heap_, parent_value, index);
      data_[index] = parent_value;
      index = parent_index;
    }
    data_[index] = currently_being_moved;
    put(index_in_heap_, currently_being_moved, index);
  }

  // From the root, swap elements (each one with its smallest child) if there
  // are any parent-child pairs that violate the heap property
  void preserve_heap_property_down() {
    if (data_.empty()) {
      return;
    }
    size_type index = 0;
    Value currently_being_moved = data_[0];
    distance_type currently_being_moved_dist = get(distance_, currently_being_moved);
    size_type heap_size = data_.size();
    auto data_span = std::span<Value>{data_};
    for (;;) {
      size_type first_child_index = first_child(index);
      if (first_child_index >= heap_size) {
        break; /* No children */
      }
      auto child_base_span = data_span.subspan(first_child_index);
      size_type smallest_child_index = 0;
      distance_type smallest_child_dist = get(distance_, child_base_span[smallest_child_index]);
      if (first_child_index + Arity <= heap_size) {
        // Special case for a statically known loop count (common case)
        for (size_t i = 1; i < Arity; ++i) {
          Value i_value = child_base_span[i];
          distance_type i_dist = get(distance_, i_value);
          if (compare_(i_dist, smallest_child_dist)) {
            smallest_child_index = i;
            smallest_child_dist = i_dist;
          }
        }
      } else {
        for (size_t i = 1; i < heap_size - first_child_index; ++i) {
          distance_type i_dist = get(distance_, child_base_span[i]);
          if (compare_(i_dist, smallest_child_dist)) {
            smallest_child_index = i;
            smallest_child_dist = i_dist;
          }
        }
      }
      if (compare_(smallest_child_dist, currently_being_moved_dist)) {
        swap_heap_elements(smallest_child_index + first_child_index, index);
        index = smallest_child_index + first_child_index;
        continue;
      }
      break;  // Heap property satisfied
    }
  }
};

template <std::semiregular Value, std::size_t Arity, concepts::ReadWritePropertyMap<Value> IndexInHeapPropertyMap,
          concepts::ReadablePropertyMap<Value> DistanceMap,
          typename Compare>
auto make_d_ary_heap_indirect(DistanceMap distance, IndexInHeapPropertyMap index_in_heap, Compare compare) {
    return d_ary_heap_indirect<Value, Arity, IndexInHeapPropertyMap, DistanceMap, Compare>(std::move(distance), std::move(index_in_heap), std::move(compare));
}

template <std::semiregular Value, std::size_t Arity, concepts::ReadWritePropertyMap<Value> IndexInHeapPropertyMap,
          concepts::ReadablePropertyMap<Value> DistanceMap>
auto make_d_ary_heap_indirect(DistanceMap distance, IndexInHeapPropertyMap index_in_heap) {
    return d_ary_heap_indirect<Value, Arity, IndexInHeapPropertyMap, DistanceMap>(std::move(distance), std::move(index_in_heap));
}

}  // namespace bagl

#endif  // BAGL_BAGL_D_ARY_HEAP_H_
