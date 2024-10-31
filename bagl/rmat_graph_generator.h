// Copyright 2004, 2005 The Trustees of Indiana University.
//  Authors: Nick Edmonds
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_RMAT_GRAPH_GENERATOR_H_
#define BAGL_BAGL_RMAT_GRAPH_GENERATOR_H_

#include <cassert>
#include <cmath>
#include <iterator>
#include <map>
#include <memory>
#include <queue>
#include <random>
#include <ranges>
#include <type_traits>
#include <utility>
#include <vector>

namespace bagl {

namespace rmat_detail {

// Returns floor(log_2(n)), and -1 when n is 0
template <typename IntegerType>
int int_log2(IntegerType n) {
  int l = 0;
  while (n > 0) {
    ++l;
    n >>= 1;
  }
  return l - 1;
}

template <std::uniform_random_bit_generator RandomGenerator, std::integral T>
void generate_permutation_vector(RandomGenerator& gen, std::vector<T>& vertexPermutation, T n) {
  vertexPermutation.resize(n);

  // Generate permutation map of vertex numbers
  std::uniform_int_distribution<T> rand_vertex(0, n - 1);
  for (T i = 0; i < n; ++i) {
    vertexPermutation[i] = i;
  }

  // Can't use std::random_shuffle unless we create another (synchronized)
  // PRNG
  for (T i = 0; i < n; ++i) {
    std::swap(vertexPermutation[i], vertexPermutation[rand_vertex(gen)]);
  }
}

template <std::uniform_random_bit_generator RandomGenerator, std::integral T>
std::pair<T, T> generate_edge(RandomGenerator& gen, std::uniform_real_distribution<> prob, T n, std::size_t scale,
                              double a, double b, double c, double d) {
  T u = 0;
  T v = 0;
  T step = n / 2;
  for (std::size_t j = 0; j < scale; ++j) {
    double p = prob(gen);

    if (p < a) {
      ;
    } else if (p >= a && p < a + b) {
      v += step;
    } else if (p >= a + b && p < a + b + c) {
      u += step;
    } else {  // p > a + b + c && p < a + b + c + d
      u += step;
      v += step;
    }

    step /= 2;

    // 0.2 and 0.9 are hardcoded in the reference SSCA implementation.
    // The maximum change in any given value should be less than 10%
    a *= 0.9 + 0.2 * prob(gen);
    b *= 0.9 + 0.2 * prob(gen);
    c *= 0.9 + 0.2 * prob(gen);
    d *= 0.9 + 0.2 * prob(gen);

    double s = a + b + c + d;

    a /= s;
    b /= s;
    c /= s;
    // d /= S;
    // Ensure all values add up to 1, regardless of floating point errors
    d = 1. - a - b - c;
  }

  return {u, v};
}

}  // namespace rmat_detail

struct keep_all_edges {
  template <typename T>
  bool operator()(const T& /*unused*/, const T& /*unused*/) {
    return true;
  }
};

template <typename Distribution, typename ProcessId>
struct keep_local_edges {
  keep_local_edges(const Distribution& distrib, const ProcessId& id) : distrib_(distrib), id_(id) {}

  template <typename T>
  bool operator()(const T& x, const T& y) {
    return distrib_(x) == id_ || distrib_(y) == id_;
  }

 private:
  const Distribution& distrib_;
  const ProcessId& id_;
};

/*
  Chakrabarti's R-MAT scale free generator.

  For all flavors of the R-MAT iterator a+b+c+d must equal 1 and for the
  unique_rmat_iterator 'm' << 'n^2'.  If 'm' is too close to 'n^2' the
  generator may be unable to generate sufficient unique edges

  To get a true scale free distribution {a, b, c, d : a > b, a > c, a > d}
*/

template <std::uniform_random_bit_generator RandomGenerator>
class rmat_iterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;  // Not used

  // No argument constructor, set to terminating condition
  rmat_iterator() = default;

  // Initialize for edge generation
  rmat_iterator(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c, double d,
                bool permute_vertices = true)
      : gen_(&gen),
        n_(n),
        a_(a),
        b_(b),
        c_(c),
        d_(d),
        edge_(m),
        permute_vertices_(permute_vertices),
        scale_(rmat_detail::int_log2(n))

  {
    if (permute_vertices) {
      rmat_detail::generate_permutation_vector(*gen_, vertex_permutation_, n);
    }

    // TODO: Generate the entire adjacency matrix then "Clip and flip" if
    // undirected graph

    // Generate the first edge
    auto [u, v] = rmat_detail::generate_edge(*gen_, prob_, n, scale_, a, b, c, d);

    if (permute_vertices) {
      current_ = {vertex_permutation_[u], vertex_permutation_[v]};
    } else {
      current_ = {u, v};
    }

    --edge_;
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  rmat_iterator& operator++() {
    auto [u, v] = rmat_detail::generate_edge(*gen_, prob_, n_, scale_, a_, b_, c_, d_);

    if (permute_vertices_) {
      current_ = {vertex_permutation_[u], vertex_permutation_[v]};
    } else {
      current_ = {u, v};
    }

    --edge_;

    return *this;
  }

  rmat_iterator operator++(int) {
    rmat_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const rmat_iterator& other) const { return edge_ == other.edge_; }

  bool operator!=(const rmat_iterator& other) const { return !(*this == other); }

 private:
  // Parameters
  RandomGenerator* gen_ = nullptr;
  std::uniform_real_distribution<> prob_;
  std::size_t n_ = 0;
  double a_ = 0.0;
  double b_ = 0.0;
  double c_ = 0.0;
  double d_ = 0.0;
  std::size_t edge_ = 0;
  bool permute_vertices_ = false;
  int scale_ = 0;

  // Internal data structures
  std::vector<std::size_t> vertex_permutation_;
  value_type current_{};
};

template <std::uniform_random_bit_generator RandomGenerator>
auto rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c, double d,
                bool permute_vertices = true) {
  return std::ranges::subrange(rmat_iterator<RandomGenerator>(gen, n, m, a, b, c, d, permute_vertices),
                               rmat_iterator<RandomGenerator>());
}

template <std::uniform_random_bit_generator RandomGenerator, typename EdgePredicate = keep_all_edges>
class sorted_rmat_iterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;  // Not used

  // No argument constructor, set to terminating condition
  sorted_rmat_iterator() : done_(true) {}

  // Initialize for edge generation
  sorted_rmat_iterator(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c, double d,
                       bool permute_vertices = true, EdgePredicate ep = keep_all_edges())
      : gen_(&gen),
        permute_vertices_(permute_vertices)

  {
    std::vector<std::size_t> vertex_permutation;
    if (permute_vertices) {
      rmat_detail::generate_permutation_vector(*gen_, vertex_permutation, n);
    }

    // TODO: "Clip and flip" if undirected graph
    int scale = rmat_detail::int_log2(n);

    for (std::size_t i = 0; i < m; ++i) {
      auto [u, v] = rmat_detail::generate_edge(*gen_, prob_, n, scale, a, b, c, d);

      if (permute_vertices) {
        if (ep(vertex_permutation[u], vertex_permutation[v])) {
          values_.emplace(vertex_permutation[u], vertex_permutation[v]);
        }
      } else {
        if (ep(u, v)) {
          values_.emplace(u, v);
        }
      }
    }

    current_ = values_.top();
    values_.pop();
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  sorted_rmat_iterator& operator++() {
    if (!values_.empty()) {
      current_ = values_.top();
      values_.pop();
    } else {
      done_ = true;
    }

    return *this;
  }

  sorted_rmat_iterator operator++(int) {
    sorted_rmat_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const sorted_rmat_iterator& other) const {
    return values_.empty() && other.values_.empty() && done_ && other.done_;
  }

  bool operator!=(const sorted_rmat_iterator& other) const { return !(*this == other); }

 private:
  // Parameters
  RandomGenerator* gen_ = nullptr;
  std::uniform_real_distribution<> prob_;
  bool permute_vertices_ = false;

  // Internal data structures
  std::priority_queue<value_type, std::vector<value_type>, std::greater<> > values_;
  value_type current_{};
  bool done_ = false;
};

template <std::uniform_random_bit_generator RandomGenerator>
auto sorted_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c, double d,
                       bool permute_vertices = true) {
  return std::ranges::subrange(sorted_rmat_iterator<RandomGenerator>(gen, n, m, a, b, c, d, permute_vertices),
                               sorted_rmat_iterator<RandomGenerator>());
}

template <std::uniform_random_bit_generator RandomGenerator, typename EdgePredicate>
auto sorted_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c, double d,
                       bool permute_vertices, EdgePredicate ep) {
  return std::ranges::subrange(
      sorted_rmat_iterator<RandomGenerator, EdgePredicate>(gen, n, m, a, b, c, d, permute_vertices, ep),
      sorted_rmat_iterator<RandomGenerator, EdgePredicate>());
}

// This version is slow but guarantees unique edges
template <std::uniform_random_bit_generator RandomGenerator, bool IsUndirected = false,
          typename EdgePredicate = keep_all_edges>
class unique_rmat_iterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;  // Not used

  // No argument constructor, set to terminating condition
  unique_rmat_iterator() : done_(true) {}

  // Initialize for edge generation
  unique_rmat_iterator(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c, double d,
                       bool permute_vertices = true, EdgePredicate ep = {})
      : gen_(&gen)

  {
    std::vector<std::size_t> vertex_permutation;
    if (permute_vertices) {
      rmat_detail::generate_permutation_vector(*gen_, vertex_permutation, n);
    }

    int scale = rmat_detail::int_log2(n);

    std::map<value_type, bool> edge_map;

    std::size_t edges = 0;
    do {
      auto [u, v] = rmat_detail::generate_edge(*gen_, prob_, n, scale, a, b, c, d);

      // Lowest vertex number always comes first
      // (this means we don't have to worry about i->j and j->i being in
      // the edge list)
      if (u > v && IsUndirected) {
        std::swap(u, v);
      }

      if (edge_map.find({u, v}) == edge_map.end()) {
        edge_map[{u, v}] = true;

        if (permute_vertices) {
          if (ep(vertex_permutation[u], vertex_permutation[v])) {
            values_.emplace_back(vertex_permutation[u], vertex_permutation[v]);
          }
        } else {
          if (ep(u, v)) {
            values_.emplace_back(u, v);
          }
        }

        edges++;
      }
    } while (edges < m);
    // NGE - Asking for more than n^2 edges will result in an infinite loop
    // here
    //       Asking for a value too close to n^2 edges may as well

    current_ = values_.back();
    values_.pop_back();
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  unique_rmat_iterator& operator++() {
    if (!values_.empty()) {
      current_ = values_.back();
      values_.pop_back();
    } else {
      done_ = true;
    }

    return *this;
  }

  unique_rmat_iterator operator++(int) {
    unique_rmat_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const unique_rmat_iterator& other) const {
    return values_.empty() && other.values_.empty() && done_ && other.done_;
  }

  bool operator!=(const unique_rmat_iterator& other) const { return !(*this == other); }

 private:
  // Parameters
  RandomGenerator* gen_;
  std::uniform_real_distribution<> prob_;

  // Internal data structures
  std::vector<value_type> values_;
  value_type current_{};
  bool done_ = false;
};

template <std::uniform_random_bit_generator RandomGenerator>
auto unique_undirected_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                                  double d, bool permute_vertices = true) {
  return std::ranges::subrange(unique_rmat_iterator<RandomGenerator, true>(gen, n, m, a, b, c, d, permute_vertices),
                               unique_rmat_iterator<RandomGenerator, true>());
}

template <std::uniform_random_bit_generator RandomGenerator, typename EdgePredicate>
auto unique_undirected_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                                  double d, bool permute_vertices, EdgePredicate ep) {
  return std::ranges::subrange(
      unique_rmat_iterator<RandomGenerator, true, EdgePredicate>(gen, n, m, a, b, c, d, permute_vertices, ep),
      unique_rmat_iterator<RandomGenerator, true, EdgePredicate>());
}

template <std::uniform_random_bit_generator RandomGenerator>
auto unique_directed_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                                double d, bool permute_vertices = true) {
  return std::ranges::subrange(unique_rmat_iterator<RandomGenerator, false>(gen, n, m, a, b, c, d, permute_vertices),
                               unique_rmat_iterator<RandomGenerator, false>());
}

template <std::uniform_random_bit_generator RandomGenerator, typename EdgePredicate>
auto unique_directed_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                                double d, bool permute_vertices, EdgePredicate ep) {
  return std::ranges::subrange(
      unique_rmat_iterator<RandomGenerator, false, EdgePredicate>(gen, n, m, a, b, c, d, permute_vertices, ep),
      unique_rmat_iterator<RandomGenerator, false, EdgePredicate>());
}

template <bool IsUndirected, typename... Args>
auto unique_rmat_range(Args&&... args) {
  if constexpr (IsUndirected) {
    return unique_undirected_rmat_range(std::forward<Args>(args)...);
  } else {
    return unique_directed_rmat_range(std::forward<Args>(args)...);
  }
}

// This version is slow but guarantees unique edges
template <std::uniform_random_bit_generator RandomGenerator, bool IsUndirected = false,
          typename EdgePredicate = keep_all_edges>
class sorted_unique_rmat_iterator {
 public:
  using iterator_category = std::input_iterator_tag;
  using value_type = std::pair<std::size_t, std::size_t>;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;  // Not used

  // No argument constructor, set to terminating condition
  sorted_unique_rmat_iterator() : done_(true) {}

  // Initialize for edge generation
  sorted_unique_rmat_iterator(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                              double d, bool bidirectional = false, bool permute_vertices = true, EdgePredicate ep = {})
      : gen_(&gen),
        bidirectional_(bidirectional)

  {
    std::vector<std::size_t> vertex_permutation;
    if (permute_vertices) {
      rmat_detail::generate_permutation_vector(*gen_, vertex_permutation, n);
    }

    int scale = rmat_detail::int_log2(n);

    std::map<value_type, bool> edge_map;

    std::size_t edges = 0;
    do {
      auto [u, v] = rmat_detail::generate_edge(*gen_, prob_, n, scale, a, b, c, d);

      if (bidirectional) {
        if (edge_map.find({u, v}) == edge_map.end()) {
          edge_map[{u, v}] = true;
          edge_map[{v, u}] = true;

          if (ep(u, v)) {
            if (permute_vertices) {
              values_.emplace(vertex_permutation[u], vertex_permutation[v]);
              values_.emplace(vertex_permutation[v], vertex_permutation[u]);
            } else {
              values_.emplace(u, v);
              values_.emplace(v, u);
            }
          }

          ++edges;
        }
      } else {
        // Lowest vertex number always comes first
        // (this means we don't have to worry about i->j and j->i being
        // in the edge list)
        if (u > v && IsUndirected) {
          std::swap(u, v);
        }

        if (edge_map.find({u, v}) == edge_map.end()) {
          edge_map[{u, v}] = true;

          if (permute_vertices) {
            if (ep(vertex_permutation[u], vertex_permutation[v])) {
              values_.emplace(vertex_permutation[u], vertex_permutation[v]);
            }
          } else {
            if (ep(u, v)) {
              values_.emplace(u, v);
            }
          }

          ++edges;
        }
      }

    } while (edges < m);
    // NGE - Asking for more than n^2 edges will result in an infinite loop
    // here
    //       Asking for a value too close to n^2 edges may as well

    current_ = values_.top();
    values_.pop();
  }

  reference operator*() const { return current_; }
  pointer operator->() const { return &current_; }

  sorted_unique_rmat_iterator& operator++() {
    if (!values_.empty()) {
      current_ = values_.top();
      values_.pop();
    } else {
      done_ = true;
    }

    return *this;
  }

  sorted_unique_rmat_iterator operator++(int) {
    sorted_unique_rmat_iterator temp(*this);
    ++(*this);
    return temp;
  }

  bool operator==(const sorted_unique_rmat_iterator& other) const {
    return values_.empty() && other.values_.empty() && done_ && other.done_;
  }

  bool operator!=(const sorted_unique_rmat_iterator& other) const { return !(*this == other); }

 private:
  // Parameters
  RandomGenerator* gen_ = nullptr;
  std::uniform_real_distribution<> prob_;
  bool bidirectional_ = false;

  // Internal data structures
  std::priority_queue<value_type, std::vector<value_type>, std::greater<> > values_;
  value_type current_{};
  bool done_ = false;
};

template <std::uniform_random_bit_generator RandomGenerator>
auto sorted_unique_undirected_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b,
                                         double c, double d, bool bidirectional = false, bool permute_vertices = true) {
  return std::ranges::subrange(
      sorted_unique_rmat_iterator<RandomGenerator, true>(gen, n, m, a, b, c, d, bidirectional, permute_vertices),
      sorted_unique_rmat_iterator<RandomGenerator, true>());
}

template <std::uniform_random_bit_generator RandomGenerator, typename EdgePredicate>
auto sorted_unique_undirected_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b,
                                         double c, double d, bool bidirectional, bool permute_vertices,
                                         EdgePredicate ep) {
  return std::ranges::subrange(sorted_unique_rmat_iterator<RandomGenerator, true, EdgePredicate>(
                                   gen, n, m, a, b, c, d, bidirectional, permute_vertices, ep),
                               sorted_unique_rmat_iterator<RandomGenerator, true, EdgePredicate>());
}

template <std::uniform_random_bit_generator RandomGenerator>
auto sorted_unique_directed_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                                       double d, bool bidirectional = false, bool permute_vertices = true) {
  return std::ranges::subrange(
      sorted_unique_rmat_iterator<RandomGenerator, false>(gen, n, m, a, b, c, d, bidirectional, permute_vertices),
      sorted_unique_rmat_iterator<RandomGenerator, false>());
}

template <std::uniform_random_bit_generator RandomGenerator, typename EdgePredicate>
auto sorted_unique_directed_rmat_range(RandomGenerator& gen, std::size_t n, std::size_t m, double a, double b, double c,
                                       double d, bool bidirectional, bool permute_vertices, EdgePredicate ep) {
  return std::ranges::subrange(sorted_unique_rmat_iterator<RandomGenerator, false, EdgePredicate>(
                                   gen, n, m, a, b, c, d, bidirectional, permute_vertices, ep),
                               sorted_unique_rmat_iterator<RandomGenerator, false, EdgePredicate>());
}

template <bool IsUndirected, typename... Args>
auto sorted_unique_rmat_range(Args&&... args) {
  if constexpr (IsUndirected) {
    return sorted_unique_undirected_rmat_range(std::forward<Args>(args)...);
  } else {
    return sorted_unique_directed_rmat_range(std::forward<Args>(args)...);
  }
}

}  // namespace bagl

#endif  // BAGL_BAGL_RMAT_GRAPH_GENERATOR_H_
