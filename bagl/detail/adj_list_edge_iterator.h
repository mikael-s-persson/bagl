// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_ADJ_LIST_EDGE_ITERATOR_H_
#define BAGL_BAGL_DETAIL_ADJ_LIST_EDGE_ITERATOR_H_

#include <iterator>
#include <optional>
#include <utility>

namespace bagl::detail {

template <class VertexIterator, class OutEdgeIterator, class Graph>
class adj_list_edge_iterator {
  using self = adj_list_edge_iterator;

 public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = typename std::iterator_traits<OutEdgeIterator>::value_type;
  using reference = typename std::iterator_traits<OutEdgeIterator>::reference;
  using pointer = typename std::iterator_traits<OutEdgeIterator>::pointer;
  using difference_type = typename std::iterator_traits<OutEdgeIterator>::difference_type;
  using distance_type = difference_type;

  adj_list_edge_iterator() = default;

  template <class G>
  adj_list_edge_iterator(VertexIterator c, VertexIterator e, const G& g)
      : curr_(c), end_(e), g_(&g) {
    if (curr_ != end_) {
      while (curr_ != end_ && out_degree(*curr_, *g_) == 0) {
        ++curr_;
      }
      if (curr_ != end_) {
        edges_ = out_edges(*curr_, *g_);
      }
    }
  }

  // Note:
  //  In the directed graph cases, it is fine.
  //  For undirected graphs, one edge go through twice.
  self& operator++() {
    ++edges_->first;
    if (edges_->first == edges_->second) {
      ++curr_;
      while (curr_ != end_ && out_degree(*curr_, *g_) == 0) {
        ++curr_;
      }
      if (curr_ != end_) {
        edges_ = out_edges(*curr_, *g_);
      }
    }
    return *this;
  }
  self operator++(int) {
    self tmp = *this;
    ++(*this);
    return tmp;
  }
  value_type operator*() const { return *edges_->first; }
  bool operator==(const self& x) const {
    return curr_ == x.curr_ && (curr_ == end_ || edges_->first == x.edges_->first);
  }
  bool operator!=(const self& x) const {
    return curr_ != x.curr_ || (curr_ != end_ && edges_->first != x.edges_->first);
  }

 protected:
  VertexIterator curr_;
  VertexIterator end_;
  std::optional<std::pair<OutEdgeIterator, OutEdgeIterator> > edges_;
  const Graph* g_;
};

}  // namespace bagl::detail

#endif  // BAGL_BAGL_DETAIL_ADJ_LIST_EDGE_ITERATOR_H_
