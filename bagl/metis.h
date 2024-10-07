// Copyright 2005 The Trustees of Indiana University.
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_METIS_H_
#define BAGL_BAGL_METIS_H_

#include <math.h>

#include <algorithm>
#include <cstddef>
#include <exception>
#include <iostream>
#include <iterator>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace bagl {

class metis_exception : public std::exception {};
class metis_input_exception : public metis_exception {};

class metis_reader {
 public:
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  using vertex_weight_type = double;
  using edge_weight_type = double;

  class edge_iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = std::pair<vertices_size_type, vertices_size_type>;
    using reference = const value_type&;
    using pointer = const value_type*;
    using difference_type = std::ptrdiff_t;

    edge_iterator& operator++() {
      advance(false);
      return *this;
    }

    std::optional<value_type> operator++(int) {
      std::optional<value_type> result = **this;
      ++(*this);
      return result;
    }

    [[nodiscard]] reference operator*() const { return self_->edge_; }
    [[nodiscard]] pointer operator->() const { return &self_->edge_; }

    [[nodiscard]] bool operator==(const edge_iterator& y) const {
      return (self_ == y.self_ || (self_ != nullptr && self_->edge_.first == self_->num_vertices()) ||
              (y.self_ != nullptr && y.self_->edge_.first == y.self_->num_vertices()));
    }

    [[nodiscard]] bool operator!=(const edge_iterator& y) const { return !(*this == y); }

    // Default copy constructor and assignment operator are okay

    // Should be private, but iterators have to be default constructible.
    edge_iterator() : self_(nullptr) {}

   private:
    explicit edge_iterator(metis_reader* self) : self_(self) {
      if (self_ != nullptr) {
        advance(true);
      }
    }
    void advance(bool skip_initial_read);

    metis_reader* self_;

    friend class metis_reader;
  };
  friend class edge_iterator;

  class edge_weight_iterator {
   public:
    using iterator_category = std::input_iterator_tag;
    using value_type = edge_weight_type;
    using reference = const value_type&;
    using pointer = const value_type*;
    using difference_type = std::ptrdiff_t;

    // Default copy constructor and assignment operator are okay

    [[nodiscard]] reference operator*() const { return self_->edge_weight_; }
    [[nodiscard]] pointer operator->() const { return &self_->edge_weight_; }

    edge_weight_iterator& operator++() { return *this; }
    edge_weight_iterator operator++(int) { return *this; }

    [[nodiscard]] bool operator==(const edge_weight_iterator& y) const { return (self_ == y.self_); }
    [[nodiscard]] bool operator!=(const edge_weight_iterator& y) const { return !(*this == y); }

    // Should be private, but iterators have to be default constructible.
    edge_weight_iterator() : self_(nullptr) {}

   private:
    explicit edge_weight_iterator(metis_reader* self) : self_(self) {}

    metis_reader* self_;

    friend class metis_reader;
  };

  explicit metis_reader(std::istream& in) : in_(in), edge_weight_(1.0) { start(); }

  [[nodiscard]] auto edges() {
    if (edge_.first != 0) {
      start();
    }
    return std::ranges::subrange{edge_iterator(this), edge_iterator(nullptr)};
  }
  [[nodiscard]] auto weights() {
    return std::ranges::subrange{edge_weight_iterator(this), edge_weight_iterator(nullptr)};
  }

  [[nodiscard]] vertices_size_type num_vertices() const { return n_vertices_; }
  [[nodiscard]] edges_size_type num_edges() const { return n_edges_; }

  [[nodiscard]] std::size_t num_vertex_weights() const { return n_vertex_weights_; }

  [[nodiscard]] vertex_weight_type vertex_weight(vertices_size_type v, std::size_t n) {
    return vertex_weights_[v * num_vertex_weights() + n];
  }

  [[nodiscard]] bool has_edge_weights() const { return edge_weights_; }

 private:
  void start();

  // Configuration information
  std::istream& in_;

  // Information about the current METIS file
  vertices_size_type n_vertices_ = 0;
  edges_size_type n_edges_ = 0;
  std::size_t n_vertex_weights_ = 0;
  bool edge_weights_ = false;

  // Information about the current edge/vertex
  std::istringstream line_in_;
  std::pair<vertices_size_type, vertices_size_type> edge_{};
  std::vector<vertex_weight_type> vertex_weights_;
  edge_weight_type edge_weight_ = 0.0;
};

class metis_distribution {
 public:
  using process_id_type = int;
  using size_type = std::size_t;
  using iterator = std::vector<process_id_type>::iterator;

  metis_distribution(std::istream& in, process_id_type my_id)
      : my_id_(my_id),
        vertices_(std::istream_iterator<process_id_type>(in), std::istream_iterator<process_id_type>()) {}

  [[nodiscard]] size_type block_size(process_id_type id, size_type) const {
    return std::count(vertices_.begin(), vertices_.end(), id);
  }
  [[nodiscard]] process_id_type operator()(size_type n) const { return vertices_[n]; }
  [[nodiscard]] size_type local(size_type n) const {
    return std::count(vertices_.begin(), vertices_.begin() + static_cast<std::ptrdiff_t>(n), vertices_[n]);
  }
  [[nodiscard]] size_type global(size_type n) const { return global(my_id_, n); }
  [[nodiscard]] size_type global(process_id_type id, size_type n) const {
    auto i = vertices_.begin();
    while (*i != id) {
      ++i;
    }

    while (n > 0) {
      do {
        ++i;
      } while (*i != id);
      --n;
    }

    return i - vertices_.begin();
  }

  [[nodiscard]] iterator begin() { return vertices_.begin(); }
  [[nodiscard]] iterator end() { return vertices_.end(); }

 private:
  process_id_type my_id_;
  std::vector<process_id_type> vertices_;
};

inline void metis_reader::edge_iterator::advance(bool skip_initial_read) {
  do {
    if (!skip_initial_read) {
      // Try to read the next edge
      if (self_->line_in_ >> std::ws >> self_->edge_.second) {
        --self_->edge_.second;
        if (self_->has_edge_weights() && !(self_->line_in_ >> self_->edge_weight_)) {
          throw metis_input_exception();
        }
        return;
      }

      // Check if we're done
      ++self_->edge_.first;
      if (self_->edge_.first == self_->num_vertices()) {
        return;
      }
    }

    // Find the next line
    std::string line;
    while (getline(self_->in_, line) && !line.empty() && line[0] == '%') {
      /* Keep reading lines in the loop header... */
    }
    if (!self_->in_) {
      throw metis_input_exception();
    }
    self_->line_in_.str(line);
    self_->line_in_.clear();

    // Read the next line
    std::size_t weights_left = self_->n_vertex_weights_;
    vertex_weight_type weight = 0.0;
    while (weights_left > 0) {
      if (self_->line_in_ >> weight) {
        self_->vertex_weights_.push_back(weight);
      } else {
        throw metis_input_exception();
      }
      --weights_left;
    }

    // Successive iterations will pick up edges for this vertex.
    skip_initial_read = false;
  } while (true);
}

inline void metis_reader::start() {
  in_.seekg(0, std::ios::beg);
  std::string line;
  while (getline(in_, line) && !line.empty() && line[0] == '%') {
    /* Keep getting lines in loop header. */
  }

  if (!in_ || line.empty()) {
    throw metis_input_exception();
  }

  // Determine number of vertices and edges in the graph
  line_in_.str(line);
  if (!(line_in_ >> n_vertices_ >> n_edges_)) {
    throw metis_input_exception();
  }

  // Determine whether vertex or edge weights are included in the graph
  int fmt = 0;
  line_in_ >> fmt;
  n_vertex_weights_ = fmt / 10;
  edge_weights_ = (fmt % 10 == 1);

  // Determine how many (if any!) vertex weights are included
  if (n_vertex_weights_ != 0) {
    line_in_ >> n_vertex_weights_;
  }

  // Setup the iteration data structures
  edge_weight_ = 1.0;
  edge_.first = 0;
  edge_.second = 0;
  vertex_weights_.reserve(n_vertex_weights_ * num_vertices());
}

}  // namespace bagl

#endif  // BAGL_BAGL_METIS_H_
