// Copyright 2004 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SMALL_WORLD_GENERATOR_H_
#define BAGL_BAGL_SMALL_WORLD_GENERATOR_H_

#include <iterator>
#include <random>
#include <ranges>
#include <utility>

namespace bagl {

// Assumes undirected
template <std::uniform_random_bit_generator RandomGenerator>
class small_world_iterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = int;

  small_world_iterator() = default;
  small_world_iterator(RandomGenerator& gen, std::size_t n, std::size_t k, double prob = 0.0,
                       bool allow_self_loops = false)
      : gen_(&gen),
        n_(n),
        k_(k),
        prob_(prob),
        target_(allow_self_loops ? 0 : 1),
        allow_self_loops_(allow_self_loops),
        current_(0, allow_self_loops ? 0 : 1) {}

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  small_world_iterator& operator++() {
    target_ = (target_ + 1) % n_;
    if (target_ == (source_ + k_ / 2 + 1) % n_) {
      ++source_;
      if (allow_self_loops_) {
        target_ = source_;
      } else {
        target_ = (source_ + 1) % n_;
      }
    }
    current_.first = source_;

    std::uniform_real_distribution<double> rand01{};
    std::uniform_int_distribution<std::size_t> rand_vertex_gen(0, n_ - 1);
    double x = rand01(*gen_);
    if (x < prob_) {
      std::size_t lower = (source_ + n_ - k_ / 2) % n_;
      std::size_t upper = (source_ + k_ / 2) % n_;
      do {
        current_.second = rand_vertex_gen(*gen_);
      } while ((current_.second >= lower && current_.second <= upper) ||
               (upper < lower && (current_.second >= lower || current_.second <= upper)));
    } else {
      current_.second = target_;
    }
    return *this;
  }

  small_world_iterator operator++(int) {
    small_world_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const small_world_iterator& other) const {
    if (!gen_ && other.gen_) {
      return other == *this;
    }
    if (gen_ && !other.gen_) {
      return source_ == n_;
    }
    if (!gen_ && !other.gen_) {
      return true;
    }
    return source_ == other.source_ && target_ == other.target_;
  }

  bool operator!=(const small_world_iterator& other) const { return !(*this == other); }

 private:
  void next() {
    std::uniform_int_distribution<std::size_t> rand_vertex(0, n_ - 1);
    current_.first = rand_vertex(*gen_);
    do {
      current_.second = rand_vertex(*gen_);
    } while (current_.first == current_.second && !allow_self_loops_);
  }

  RandomGenerator* gen_ = nullptr;
  std::size_t n_ = 0;
  std::size_t k_ = 0;
  double prob_ = 0.0;
  std::size_t source_ = 0;
  std::size_t target_ = 0;
  bool allow_self_loops_ = false;
  value_type current_{};
};

template <std::uniform_random_bit_generator RandomGenerator>
auto small_world_range(RandomGenerator& gen, std::size_t n, std::size_t k, double prob = 0.0,
                       bool allow_self_loops = false) {
  return std::ranges::subrange(small_world_iterator<RandomGenerator>(gen, n, k, prob, allow_self_loops),
                               small_world_iterator<RandomGenerator>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_SMALL_WORLD_GENERATOR_H_
