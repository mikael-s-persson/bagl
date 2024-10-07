// Copyright 2004, 2005 The Trustees of Indiana University.
//  Authors: Nick Edmonds
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SSCA_GRAPH_GENERATOR_H_
#define BAGL_BAGL_SSCA_GRAPH_GENERATOR_H_

#include <iterator>
#include <queue>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

namespace bagl {

enum class ssca_direction { forward = 1, backward = 2, both = forward | backward };

// This generator generates graphs according to the method specified
// in SSCA 1.1.  Current versions of SSCA use R-MAT graphs

template <std::uniform_random_bit_generator RandomGenerator>
class ssca_iterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = void;

  // No argument constructor, set to terminating condition
  ssca_iterator() = default;

  // Initialize for edge generation
  ssca_iterator(RandomGenerator& gen, std::size_t tot_vertices, std::size_t max_clique_size, double prob_unidirectional,
                int max_parallel_edges, double prob_interclique_edges)
      : gen_(&gen),
        tot_vertices_(tot_vertices),
        max_clique_size_(max_clique_size),
        prob_unidirectional_(prob_unidirectional),
        max_parallel_edges_(max_parallel_edges),
        prob_interclique_edges_(prob_interclique_edges),
        vertices_remaining_(tot_vertices_) {
    clique_num_ = std::vector<int>(tot_vertices_, -1);
    current_ = std::make_pair(0, 0);
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  ssca_iterator& operator++() {
    while (values_.empty() && vertices_remaining_ > 0) {  // If there are no values left, generate a new clique
      std::uniform_int_distribution<std::size_t> clique_size(1, max_clique_size_);
      std::uniform_int_distribution<std::size_t> rand_vertex(0, tot_vertices_ - 1);
      std::uniform_int_distribution<int> num_parallel_edges(1, max_parallel_edges_);
      std::uniform_int_distribution<std::int16_t> direction(0, 1);
      std::uniform_real_distribution<> prob{};
      std::vector<std::size_t> clique_vertices;

      clique_vertices.clear();
      std::size_t size = std::min(clique_size(*gen_), vertices_remaining_);
      while (clique_vertices.size() < size) {
        std::size_t v = rand_vertex(*gen_);
        if (clique_num_[v] == -1) {
          clique_num_[v] = current_clique_;
          clique_vertices.push_back(v);
          vertices_remaining_--;
        }
      }  // Nick: This is inefficient when only a few vertices remain...
         //       I should probably just select the remaining vertices
         //       in order when only a certain fraction remain.

      for (auto first = clique_vertices.begin(); first != clique_vertices.end(); ++first) {
        for (auto second = first + 1; second != clique_vertices.end(); ++second) {
          int edges = 0;

          ssca_direction d = prob(*gen_) < prob_unidirectional_
                                 ? (direction(*gen_) == 0 ? ssca_direction::forward : ssca_direction::backward)
                                 : ssca_direction::both;

          if ((d == ssca_direction::forward) || (d == ssca_direction::both)) {
            edges = num_parallel_edges(*gen_);
            for (int i = 0; i < edges; ++i) {
              values_.emplace(*first, *second);
            }
          }

          if ((d == ssca_direction::backward) || (d == ssca_direction::both)) {
            edges = num_parallel_edges(*gen_);
            for (int i = 0; i < edges; ++i) {
              values_.emplace(*second, *first);
            }
          }
        }
      }

      if (vertices_remaining_ == 0) {
        // Generate interclique edges
        for (std::size_t i = 0; i < tot_vertices_; ++i) {
          double p = prob_interclique_edges_;
          for (std::size_t d = 2; d < tot_vertices_ / 2; d *= 2, p /= 2) {
            std::size_t j = (i + d) % tot_vertices_;
            if (clique_num_[j] != clique_num_[i] && prob(*gen_) < p) {
              int edges = num_parallel_edges(*gen_);
              for (int k = 0; k < edges; ++k) {
                values_.emplace(k, j);
              }
            }
          }
        }
      }

      current_clique_++;
    }

    if (!values_.empty()) {  // If we're not done return a value
      current_ = values_.front();
      values_.pop();
    }

    return *this;
  }

  ssca_iterator operator++(int) {
    ssca_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const ssca_iterator& other) const {
    return vertices_remaining_ == other.vertices_remaining_ && values_.empty() && other.values_.empty();
  }

  bool operator!=(const ssca_iterator& other) const { return !(*this == other); }

 private:
  // Parameters
  RandomGenerator* gen_ = nullptr;
  std::size_t tot_vertices_ = 0;
  std::size_t max_clique_size_ = 0;
  double prob_unidirectional_ = 0.0;
  int max_parallel_edges_ = 0;
  double prob_interclique_edges_ = 0.0;

  // Internal data structures
  std::vector<int> clique_num_;
  std::queue<value_type> values_;
  int current_clique_ = 0;
  std::size_t vertices_remaining_ = 0;
  value_type current_{};
};

template <std::uniform_random_bit_generator RandomGenerator>
auto small_world_range(RandomGenerator& gen, std::size_t tot_vertices, std::size_t max_clique_size,
                       double prob_unidirectional, int max_parallel_edges, double prob_interclique_edges) {
  return std::ranges::subrange(ssca_iterator<RandomGenerator>(gen, tot_vertices, max_clique_size, prob_unidirectional,
                                                              max_parallel_edges, prob_interclique_edges),
                               ssca_iterator<RandomGenerator>());
}

}  // namespace bagl

#endif  // BAGL_BAGL_SSCA_GRAPH_GENERATOR_H_
