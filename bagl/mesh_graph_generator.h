// Copyright 2004, 2005 The Trustees of Indiana University.
//  Authors: Nick Edmonds
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MESH_GRAPH_GENERATOR_H_
#define BAGL_BAGL_MESH_GRAPH_GENERATOR_H_

#include <cassert>
#include <iterator>
#include <ranges>
#include <type_traits>
#include <utility>

#include "bagl/graph_traits.h"

namespace bagl {

template <typename Graph>
class mesh_iterator {
  using directed_category = graph_directed_category_t<Graph>;

  static constexpr bool is_undirected_v = is_undirected_graph_v<Graph>;

 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = void;

  mesh_iterator() : done_(true) {}

  // Vertices are numbered in row-major order
  // Assumes directed
  mesh_iterator(std::size_t x, std::size_t y, bool toroidal = true)
      : x_(x), y_(y), n_(x * y), target_(1), current_(0, 1), toroidal_(toroidal) {
    assert(x > 1 && y > 1);
  }

  [[nodiscard]] reference operator*() const { return current_; }
  [[nodiscard]] pointer operator->() const { return &current_; }

  mesh_iterator& operator++() {
    if constexpr (is_undirected_v) {
      if (!toroidal_) {
        if (target_ == source_ + 1) {
          if (source_ < x_ * (y_ - 1)) {
            target_ = source_ + x_;
          } else {
            source_++;
            target_ = (source_ % x_) < x_ - 1 ? source_ + 1 : source_ + x_;
            if (target_ > n_) {
              done_ = true;
            }
          }
        } else if (target_ == source_ + x_) {
          source_++;
          target_ = (source_ % x_) < x_ - 1 ? source_ + 1 : source_ + x_;
        }
      } else {
        if (target_ == source_ + 1 || target_ == source_ - (source_ % x_)) {
          target_ = (source_ + x_) % n_;
        } else if (target_ == (source_ + x_) % n_) {
          if (source_ == n_ - 1) {
            done_ = true;
          } else {
            source_++;
            target_ = (source_ % x_) < (x_ - 1) ? source_ + 1 : source_ - (source_ % x_);
          }
        }
      }
    } else {  // Directed
      if (!toroidal_) {
        if (target_ == source_ - x_) {
          target_ = source_ % x_ > 0 ? source_ - 1 : source_ + 1;
        } else if (target_ == source_ - 1) {
          if ((source_ % x_) < (x_ - 1)) {
            target_ = source_ + 1;
          } else if (source_ < x_ * (y_ - 1)) {
            target_ = source_ + x_;
          } else {
            done_ = true;
          }
        } else if (target_ == source_ + 1) {
          if (source_ < x_ * (y_ - 1)) {
            target_ = source_ + x_;
          } else {
            source_++;
            target_ = source_ - x_;
          }
        } else if (target_ == source_ + x_) {
          source_++;
          target_ = (source_ >= x_) ? source_ - x_ : source_ - 1;
        }
      } else {
        if (source_ == n_ - 1 && target_ == (source_ + x_) % n_) {
          done_ = true;
        } else if (target_ == source_ - 1 || target_ == source_ + x_ - 1) {
          target_ = (source_ + x_) % n_;
        } else if (target_ == source_ + 1 || target_ == source_ - (source_ % x_)) {
          target_ = (source_ - x_ + n_) % n_;
        } else if (target_ == (source_ - x_ + n_) % n_) {
          target_ = (source_ % x_ > 0) ? source_ - 1 : source_ + x_ - 1;
        } else if (target_ == (source_ + x_) % n_) {
          source_++;
          target_ = (source_ % x_) < (x_ - 1) ? source_ + 1 : source_ - (source_ % x_);
        }
      }
    }

    current_.first = source_;
    current_.second = target_;

    return *this;
  }

  mesh_iterator operator++(int) {
    mesh_iterator temp(*this);
    ++(*this);
    return temp;
  }

  [[nodiscard]] bool operator==(const mesh_iterator& other) const { return done_ == other.done_; }
  [[nodiscard]] bool operator!=(const mesh_iterator& other) const { return !(*this == other); }

 private:
  std::size_t x_ = 0;
  std::size_t y_ = 0;
  std::size_t n_ = 0;
  std::size_t source_ = 0;
  std::size_t target_ = 0;
  value_type current_{};
  bool toroidal_ = false;
  bool done_ = false;
};

template <typename Graph>
auto mesh_range(std::size_t x, std::size_t y, bool toroidal = true) {
  return std::ranges::subrange(mesh_iterator<Graph>(x, y, toroidal), mesh_iterator<Graph>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_MESH_GRAPH_GENERATOR_H_
