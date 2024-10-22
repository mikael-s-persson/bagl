// Copyright 2004, 2005 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ERDOS_RENYI_GENERATOR_H_
#define BAGL_BAGL_ERDOS_RENYI_GENERATOR_H_

#include <cassert>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <memory>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>

#include "bagl/graph_traits.h"

namespace bagl {

template <typename RandomGenerator, typename Graph>
class erdos_renyi_iterator {
  using directed_category = graph_directed_category_t<Graph>;
  using self = erdos_renyi_iterator<RandomGenerator, Graph>;

  static constexpr bool is_undirected = std::is_base_of_v<undirected_tag, directed_category>;

 public:
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::input_iterator_tag;

  erdos_renyi_iterator() = default;
  erdos_renyi_iterator(RandomGenerator gen, std::size_t n, double fraction = 0.0, bool allow_self_loops = false)
      : gen_(std::move(gen)),
        n_(n),
        edges_(static_cast<std::size_t>(fraction * n * n)),
        allow_self_loops_(allow_self_loops) {
    if constexpr (is_undirected) {
      edges_ = edges_ / 2;
    }
    next();
  }

  erdos_renyi_iterator(RandomGenerator gen, std::size_t n, std::size_t m, bool allow_self_loops = false)
      : gen_(std::move(gen)), n_(n), edges_(m), allow_self_loops_(allow_self_loops) {
    next();
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  self& operator++() {
    --edges_;
    next();
    return *this;
  }
  self operator++(int) {
    self t = *this;
    ++(*this);
    return t;
  }
  bool operator==(const self& x) const { return edges_ == x.edges_; }
  bool operator!=(const self& x) const { return edges_ != x.edges_; }

 private:
  void next() {
    std::uniform_int_distribution<std::size_t> rand_vertex(0, n_ - 1);
    current_.first = rand_vertex(gen_);
    do {
      current_.second = rand_vertex(gen_);
    } while (current_.first == current_.second && !allow_self_loops_);
  }

  RandomGenerator gen_;
  std::size_t n_ = 0;
  std::size_t edges_ = 0;
  bool allow_self_loops_ = false;
  std::pair<std::size_t, std::size_t> current_ = {0, 0};
};

template <typename RandomGenerator, typename Graph>
class sorted_erdos_renyi_iterator {
  using directed_category = graph_directed_category_t<Graph>;
  using self = sorted_erdos_renyi_iterator<RandomGenerator, Graph>;

  static constexpr bool is_undirected = std::is_base_of_v<undirected_tag, directed_category>;

 public:
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;
  using iterator_category = std::input_iterator_tag;

  sorted_erdos_renyi_iterator() = default;

  // NOTE: The default probability has been changed to be the same as that
  // used by the geometic distribution. It was previously 0.0, which would
  // cause an assertion.
  sorted_erdos_renyi_iterator(RandomGenerator gen, std::size_t n, double prob = 0.5, bool loops = false)
      : gen_(std::move(gen)), rand_vertex_(1. - prob), n_(n), allow_self_loops_(loops), prob_(prob) {
    if (prob != 0.0) {
      src_ = 0;
      return;
    }
    next();
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  bool operator==(const self& x) const { return src_ == x.src_ && tgt_index_ == x.tgt_index_; }
  bool operator!=(const self& x) const { return !(*this == x); }

  self& operator++() {
    next();
    return *this;
  }
  self operator++(int) {
    self tmp = *this;
    next();
    return tmp;
  }

 private:
  void next() {
    // In order to get the edges from the generator in sorted order, one
    // effective (but slow) procedure would be to use a
    // bernoulli_distribution for each legal (src, tgt_index) pair.  Because
    // of the O(|V|^2) cost of that, a geometric distribution is used.  The
    // geometric distribution tells how many times the
    // bernoulli_distribution would need to be run until it returns true.
    // Thus, this distribution can be used to step through the edges
    // which are actually present.
    assert(src_ != std::numeric_limits<std::size_t>::max());
    assert(src_ != n_);
    while (src_ != n_) {
      std::size_t increment = rand_vertex_(gen_);
      size_t tgt_index_limit = (is_undirected ? src_ + 1 : n_) + (allow_self_loops_ ? 0 : -1);
      if (tgt_index_ + increment >= tgt_index_limit) {
        // Overflowed this source; go to the next one and try again.
        ++src_;
        // This bias is because the geometric distribution always
        // returns values >=1, and we want to allow 0 as a valid target.
        tgt_index_ = std::numeric_limits<std::size_t>::max();
        continue;
      }

      tgt_index_ += increment;
      current_.first = src_;
      current_.second = tgt_index_ + (!allow_self_loops_ && !is_undirected && tgt_index_ >= src_ ? 1 : 0);
      break;
    }
    if (src_ == n_) {
      src_ = std::numeric_limits<std::size_t>::max();
    }
  }

  RandomGenerator gen_;
  std::geometric_distribution<std::size_t> rand_vertex_{0.5};
  std::size_t n_ = 0;
  bool allow_self_loops_ = false;
  std::size_t src_ = std::numeric_limits<std::size_t>::max();
  std::size_t tgt_index_ = std::numeric_limits<std::size_t>::max();
  std::pair<std::size_t, std::size_t> current_ = {0, 0};
  double prob_ = 0.5;
};

}  // namespace bagl

#endif  // BAGL_BAGL_ERDOS_RENYI_GENERATOR_H_
