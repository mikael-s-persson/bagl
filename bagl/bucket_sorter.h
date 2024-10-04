// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

// Revision History:
//   13 June 2001: Changed some names for clarity. (Jeremy Siek)
//   01 April 2001: Modified to use new <boost/limits.hpp> header. (JMaddock)

#ifndef BAGL_BAGL_BUCKET_SORTER_H_
#define BAGL_BAGL_BUCKET_SORTER_H_

#include <cassert>
#include <limits>
#include <vector>

#include "bagl/property_map.h"

namespace bagl {

template <class BucketType, class ValueType, class Bucket, concepts::ReadablePropertyMap<ValueType> ValueIndexMap>
class bucket_sorter {
 public:
  using bucket_type = BucketType;
  using value_type = ValueType;
  using size_type = std::size_t;

  bucket_sorter(size_type length, bucket_type max_bucket, const Bucket& bucket = Bucket(),
                const ValueIndexMap& id = ValueIndexMap())
      : head_(max_bucket, invalid_value()),
        next_(length, invalid_value()),
        prev_(length, invalid_value()),
        id_to_value_(length),
        bucket_(bucket),
        id_(id) {}

  void remove(const value_type& x) {
    const size_type i = get(id_, x);
    const size_type& next_node = next_[i];
    const size_type& prev_node = prev_[i];

    // check if i is the end of the bucket list
    if (next_node != invalid_value()) {
      prev_[next_node] = prev_node;
    }
    // check if i is the begin of the bucket list
    if (prev_node != invalid_value()) {
      next_[prev_node] = next_node;
    } else {
      // need update head of current bucket list
      head_[bucket_[x]] = next_node;
    }
  }

  void push(const value_type& x) {
    id_to_value_[get(id_, x)] = x;
    (*this)[bucket_[x]].push(x);
  }

  void update(const value_type& x) {
    remove(x);
    (*this)[bucket_[x]].push(x);
  }
  //  private:
  //    with KCC, the nested stack class is having access problems
  //    despite the friend decl.
  static size_type invalid_value() { return std::numeric_limits<size_type>::max(); }

  using Iter = typename std::vector<size_type>::iterator;
  using IndexValueMap = typename std::vector<value_type>::iterator;

  class stack {
   public:
    stack(bucket_type bucket_id, Iter h, Iter n, Iter p, IndexValueMap v, const ValueIndexMap& id)
        : bucket_id_(bucket_id), head_(h), next_(n), prev_(p), value_(v), id_(id) {}

    // Avoid using default arg for ValueIndexMap so that the default
    // constructor of the ValueIndexMap is not required if not used.
    stack(bucket_type bucket_id, Iter h, Iter n, Iter p, IndexValueMap v)
        : bucket_id_(bucket_id), head_(h), next_(n), prev_(p), value_(v) {}

    void push(const value_type& x) {
      const size_type new_head = get(id_, x);
      const size_type current = head_[bucket_id_];
      if (current != invalid_value()) {
        prev_[current] = new_head;
      }
      prev_[new_head] = invalid_value();
      next_[new_head] = current;
      head_[bucket_id_] = new_head;
    }
    void pop() {
      size_type current = head_[bucket_id_];
      size_type next_node = next_[current];
      head_[bucket_id_] = next_node;
      if (next_node != invalid_value()) {
        prev_[next_node] = invalid_value();
      }
    }
    [[nodiscard]] value_type& top() { return value_[head_[bucket_id_]]; }
    [[nodiscard]] const value_type& top() const { return value_[head_[bucket_id_]]; }
    [[nodiscard]] bool empty() const { return head_[bucket_id_] == invalid_value(); }

   private:
    bucket_type bucket_id_;
    Iter head_;
    Iter next_;
    Iter prev_;
    IndexValueMap value_;
    ValueIndexMap id_;
  };

  stack operator[](const bucket_type& i) {
    assert(i < head_.size());
    return stack(i, head_.begin(), next_.begin(), prev_.begin(), id_to_value_.begin(), id_);
  }

 protected:
  std::vector<size_type> head_;
  std::vector<size_type> next_;
  std::vector<size_type> prev_;
  std::vector<value_type> id_to_value_;
  Bucket bucket_;
  ValueIndexMap id_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_BUCKET_SORTER_H_
