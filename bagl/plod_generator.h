// Copyright 2004-2006 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLOD_GENERATOR_H_
#define BAGL_BAGL_PLOD_GENERATOR_H_

#include <cmath>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <type_traits>
#include <utility>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"

namespace bagl {
template <std::uniform_random_bit_generator RandomGenerator>
class out_directed_plod_iterator {
 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;

  out_directed_plod_iterator() : at_end_(true) {}

  out_directed_plod_iterator(RandomGenerator& gen, std::size_t n, double alpha, double beta, bool allow_self_loops)
      : gen_(&gen), n_(n), alpha_(alpha), beta_(beta), allow_self_loops_(allow_self_loops), current_(0, 0) {
    using std::pow;

    std::uniform_int_distribution<std::size_t> x(0, n - 1);
    std::size_t xv = x(gen);
    degree_ = (xv == 0 ? 0 : std::size_t{beta * pow(xv, -alpha)});
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  out_directed_plod_iterator& operator++() {
    using std::pow;

    std::uniform_int_distribution<std::size_t> x(0, n_ - 1);

    // Continue stepping through source nodes until the
    // (out)degree is > 0
    while (degree_ == 0) {
      // Step to the next source node. If we've gone past the
      // number of nodes we're responsible for, we're done.
      if (++current_.first >= n_) {
        at_end_ = true;
        return *this;
      }

      std::size_t xv = x(*gen_);
      degree_ = (xv == 0 ? 0 : std::size_t{beta_ * pow(xv, -alpha_)});
    }

    do {
      current_.second = x(*gen_);
    } while (current_.first == current_.second && !allow_self_loops_);
    --degree_;

    return *this;
  }

  out_directed_plod_iterator operator++(int) {
    out_directed_plod_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const out_directed_plod_iterator& other) const { return at_end_ == other.at_end_; }

  bool operator!=(const out_directed_plod_iterator& other) const { return !(*this == other); }

 private:
  RandomGenerator* gen_ = nullptr;
  std::size_t n_ = 0;
  double alpha_ = 0.0;
  double beta_ = 0.0;
  bool allow_self_loops_ = false;
  bool at_end_ = false;
  std::size_t degree_ = 0;
  value_type current_{};
};

template <std::uniform_random_bit_generator RandomGenerator>
class undirected_plod_iterator {
  using out_degrees_t = std::vector<std::pair<std::size_t, std::size_t>>;

 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;

  undirected_plod_iterator() = default;

  undirected_plod_iterator(RandomGenerator& gen, std::size_t n, double alpha, double beta,
                           bool allow_self_loops = false)
      : gen_(&gen), n_(n), out_degrees_(new out_degrees_t), allow_self_loops_(allow_self_loops) {
    using std::pow;

    std::uniform_int_distribution<std::size_t> x(0, n - 1);
    for (std::size_t i = 0; i != n; ++i) {
      std::size_t xv = x(gen);
      std::size_t degree = (xv == 0 ? 0 : std::size_t(beta * pow(xv, -alpha)));
      if (degree == 0) {
        degree = 1;
      } else if (degree >= n) {
        degree = n - 1;
      }
      out_degrees_->push_back(std::make_pair(i, degree));
      degrees_left_ += degree;
    }

    next();
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  undirected_plod_iterator& operator++() {
    next();
    return *this;
  }

  undirected_plod_iterator operator++(int) {
    undirected_plod_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const undirected_plod_iterator& other) const { return degrees_left_ == other.degrees_left_; }

  bool operator!=(const undirected_plod_iterator& other) const { return !(*this == other); }

 private:
  void next() {
    std::size_t source = 0;
    std::size_t target = 0;
    while (true) {
      /* We may get to the point where we can't actually find any
         new edges, so we just add some random edge and set the
         degrees left = 0 to signal termination. */
      if (out_degrees_->size() < 2) {
        std::uniform_int_distribution<std::size_t> x(0, n_ - 1);
        current_.first = x(*gen_);
        do {
          current_.second = x(*gen_);
        } while (current_.first == current_.second && !allow_self_loops_);
        degrees_left_ = 0;
        out_degrees_->clear();
        return;
      }

      std::uniform_int_distribution<std::size_t> x(0, out_degrees_->size() - 1);

      // Select source vertex
      source = x(*gen_);
      if ((*out_degrees_)[source].second == 0) {
        (*out_degrees_)[source] = out_degrees_->back();
        out_degrees_->pop_back();
        continue;
      }

      // Select target vertex
      target = x(*gen_);
      if ((*out_degrees_)[target].second == 0) {
        (*out_degrees_)[target] = out_degrees_->back();
        out_degrees_->pop_back();
        continue;
      }
      if (source != target || (allow_self_loops_ && (*out_degrees_)[source].second > 2)) {
        break;
      }
    }

    // Update degree counts
    --(*out_degrees_)[source].second;
    --degrees_left_;
    --(*out_degrees_)[target].second;
    --degrees_left_;
    current_.first = (*out_degrees_)[source].first;
    current_.second = (*out_degrees_)[target].first;
  }

  RandomGenerator* gen_ = nullptr;
  std::size_t n_ = 0;
  std::shared_ptr<out_degrees_t> out_degrees_;
  std::size_t degrees_left_ = 0;
  bool allow_self_loops_ = false;
  value_type current_{};
};

template <concepts::Graph G, std::uniform_random_bit_generator RandomGenerator>
auto plod_range(RandomGenerator& gen, std::size_t n, double alpha, double beta, bool allow_self_loops = false) {
  if constexpr (is_directed_graph_v<G>) {
    return std::ranges::subrange(out_directed_plod_iterator<RandomGenerator>(gen, n, alpha, beta, allow_self_loops),
                                 out_directed_plod_iterator<RandomGenerator>());
  } else {
    return std::ranges::subrange(undirected_plod_iterator<RandomGenerator>(gen, n, alpha, beta, allow_self_loops),
                                 undirected_plod_iterator<RandomGenerator>());
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_PLOD_GENERATOR_H_
