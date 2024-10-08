// Copyright 2005 The Trustees of Indiana University.
//  Authors: Alex Breuer
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DIMACS_H_
#define BAGL_BAGL_DIMACS_H_

#include <cassert>
#include <charconv>
#include <exception>
#include <iostream>
#include <iterator>
#include <queue>
#include <ranges>
#include <sstream>
#include <string>
#include <vector>

namespace bagl {

class dimacs_exception : public std::exception {};

class dimacs_basic_reader {
 public:
  using vertices_size_type = std::size_t;
  using edges_size_type = std::size_t;
  using vertex_weight_type = double;
  using edge_weight_type = double;
  using edge_type = std::pair<vertices_size_type, vertices_size_type>;
  enum class incr_mode { edge, edge_weight };

  explicit dimacs_basic_reader(std::istream& in, bool want_weights = true) : inpt_(in), want_weights_(want_weights) {
    while (getline(inpt_, buf_) && !buf_.empty() && buf_[0] == 'c') {
    }

    if (buf_[0] != 'p') {
      throw dimacs_exception();
    }

    std::stringstream instr(buf_);
    std::string junk;

    instr >> junk >> junk >> num_vertices_ >> num_edges_;
    read_edge_weights_.push(-1);
    incr(incr_mode::edge_weight);
  }

  [[nodiscard]] const edge_type& edge_deref() const {
    assert(!read_edges_.empty());
    return read_edges_.front();
  }

  [[nodiscard]] const edge_type* edge_ref() const {
    assert(!read_edges_.empty());
    return &read_edges_.front();
  }

  [[nodiscard]] const edge_weight_type& edge_weight_deref() const {
    assert(!read_edge_weights_.empty());
    return read_edge_weights_.front();
  }

  [[nodiscard]] const edge_weight_type* edge_weight_ref() const {
    assert(!read_edge_weights_.empty());
    return &read_edge_weights_.front();
  }

  dimacs_basic_reader incr(incr_mode mode) {
    if (mode == incr_mode::edge) {
      assert(!read_edges_.empty());
      read_edges_.pop();
    } else if (mode == incr_mode::edge_weight) {
      assert(!read_edge_weights_.empty());
      read_edge_weights_.pop();
    }

    if ((mode == incr_mode::edge && read_edges_.empty()) ||
        (mode == incr_mode::edge_weight && read_edge_weights_.empty())) {
      if (seen_edges_ > num_edges_) {
        throw dimacs_exception();
      }

      while (getline(inpt_, buf_) && !buf_.empty() && buf_[0] == 'c') {
      }

      if (!inpt_.eof()) {
        int source = 0;
        int dest = 0;
        int weight = 0;
        read_edge_line(buf_, source, dest, weight);

        seen_edges_++;
        source--;
        dest--;

        read_edges_.push(edge_type(source, dest));
        if (want_weights_) {
          read_edge_weights_.push(weight);
        }
      }
      assert(read_edges_.size() < 100);
      assert(read_edge_weights_.size() < 100);
    }

    // the 1000000 just happens to be about how many edges can be read
    // in 10s
    //     if( !(seen_edges % 1000000) && !process_id( pg ) && mode ==
    //     edge ) {
    //       std::cout << "read " << seen_edges << " edges" <<
    //       std::endl;
    //     }
    return *this;
  }

  [[nodiscard]] bool done_edges() const { return inpt_.eof() && read_edges_.empty(); }

  [[nodiscard]] bool done_edge_weights() const { return inpt_.eof() && read_edge_weights_.empty(); }

  [[nodiscard]] vertices_size_type n_vertices() const { return num_vertices_; }

  [[nodiscard]] vertices_size_type processed_edges() const { return seen_edges_ - read_edges_.size(); }

  [[nodiscard]] vertices_size_type processed_edge_weights() const { return seen_edges_ - read_edge_weights_.size(); }

  [[nodiscard]] vertices_size_type n_edges() const { return num_edges_; }

 protected:
  static bool read_edge_line(std::string_view line, int& from, int& to, int& weight) {
    auto next_str_element = [&line](int& base_i) {
      int i = base_i;
      for (; i < line.size() && line[i] != '\n' && line[i] != '\0'; ++i) {
        if (line[i] == ' ') {
          break;
        }
      }
      auto res = line.substr(base_i, i - base_i);
      base_i = i + 1;
      return res;
    };

    switch (line[0]) {
      case 'e': {
        int base_i = 2;
        std::string_view from_str = next_str_element(base_i);
        std::string_view to_str = next_str_element(base_i);
        if (from_str.empty() || to_str.empty() ||
            std::from_chars(from_str.data(), from_str.data() + from_str.size(), from).ec != std::errc{} ||
            std::from_chars(to_str.data(), to_str.data() + to_str.size(), to).ec != std::errc{}) {
          return false;
        }
        weight = 0;
        return true;
      }
      case 'a': {
        int base_i = 2;
        std::string_view from_str = next_str_element(base_i);
        std::string_view to_str = next_str_element(base_i);
        std::string_view weight_str = next_str_element(base_i);
        return !(from_str.empty() || to_str.empty() || weight_str.empty() ||
                 std::from_chars(from_str.data(), from_str.data() + from_str.size(), from).ec != std::errc{} ||
                 std::from_chars(to_str.data(), to_str.data() + to_str.size(), to).ec != std::errc{} ||
                 std::from_chars(weight_str.data(), weight_str.data() + weight_str.size(), weight).ec != std::errc{});
      }
      default:
        return false;
    }
  }

  std::queue<edge_type> read_edges_;
  std::queue<edge_weight_type> read_edge_weights_;

  std::istream& inpt_;
  std::string buf_;
  vertices_size_type num_vertices_ = 0;
  vertices_size_type num_edges_ = 0;
  vertices_size_type seen_edges_ = 0;
  bool want_weights_ = true;
};

template <typename T>
class dimacs_edge_iterator {
 public:
  using edge_type = dimacs_basic_reader::edge_type;
  using incr_mode = dimacs_basic_reader::incr_mode;

  using iterator_category = std::input_iterator_tag;
  using value_type = edge_type;
  using reference = const value_type&;
  using pointer = const value_type*;
  using difference_type = std::ptrdiff_t;

  dimacs_edge_iterator() = default;
  explicit dimacs_edge_iterator(T& reader) : reader_(&reader) {}

  dimacs_edge_iterator& operator++() {
    reader_->incr(incr_mode::edge);
    return *this;
  }
  std::optional<value_type> operator++(int) {
    // Using optional as a proxy type (with dereference operators).
    std::optional<value_type> result = reader_->edge_deref();
    ++(*this);
    return result;
  }

  [[nodiscard]] reference operator*() const { return reader_->edge_deref(); }

  [[nodiscard]] pointer operator->() const { return reader_->edge_ref(); }

  // don't expect this to do the right thing if you're not comparing
  // against a general past-the-end-iterator (default constructed)
  [[nodiscard]] bool operator==(const dimacs_edge_iterator& arg) const {
    // If both are at the end (done or null), then they are equal, otherwise not.
    return ((reader_ == nullptr || reader_->done_edges()) && (arg.reader_ == nullptr || arg.reader_->done_edges()));
  }

  [[nodiscard]] bool operator!=(const dimacs_edge_iterator& arg) const { return !(*this == arg); }

 private:
  T* reader_ = nullptr;
};

inline auto dimacs_edge_range(dimacs_basic_reader& reader) {
  return std::ranges::subrange(dimacs_edge_iterator<dimacs_basic_reader>(reader),
                               dimacs_edge_iterator<dimacs_basic_reader>());
}

template <typename T>
class dimacs_edge_weight_iterator {
 public:
  using edge_weight_type = dimacs_basic_reader::edge_weight_type;
  using incr_mode = dimacs_basic_reader::incr_mode;

  using iterator_category = std::input_iterator_tag;
  using value_type = edge_weight_type;
  using reference = value_type;
  using pointer = const edge_weight_type*;
  using difference_type = std::ptrdiff_t;

  dimacs_edge_weight_iterator() = default;
  explicit dimacs_edge_weight_iterator(T& reader) : reader_(&reader) {}

  dimacs_edge_weight_iterator& operator++() {
    reader_->incr(incr_mode::edge_weight);
    return *this;
  }
  std::optional<value_type> operator++(int) {
    // Using optional as a proxy type (with dereference operators).
    std::optional<value_type> result = reader_->edge_weight_deref();
    ++(*this);
    return result;
  }

  [[nodiscard]] edge_weight_type operator*() const { return reader_->edge_weight_deref(); }
  [[nodiscard]] const edge_weight_type* operator->() const { return reader_->edge_weight_ref(); }

  // don't expect this to do the right thing if you're not comparing
  // against a general past-the-end-iterator made with the default
  // constructor for dimacs_basic_reader
  [[nodiscard]] bool operator==(const dimacs_edge_weight_iterator& arg) const {
    // If both are at the end (done or null), then they are equal, otherwise not.
    return ((reader_ == nullptr || reader_->done_edge_weights()) &&
            (arg.reader_ == nullptr || arg.reader_->done_edge_weights()));
  }

  [[nodiscard]] bool operator!=(const dimacs_edge_weight_iterator& arg) const { return !(*this == arg); }

 private:
  T* reader_ = nullptr;
};

inline auto dimacs_edge_weight_range(dimacs_basic_reader& reader) {
  return std::ranges::subrange(dimacs_edge_weight_iterator<dimacs_basic_reader>(reader),
                               dimacs_edge_weight_iterator<dimacs_basic_reader>());
}

}  // namespace bagl
#endif  // BAGL_BAGL_DIMACS_H_
