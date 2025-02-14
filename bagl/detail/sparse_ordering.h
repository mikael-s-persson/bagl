// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Copyright 2004, 2005 Trustees of Indiana University
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek,
//          Doug Gregor, D. Kevin McGrath
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_SPARSE_ORDERING_H_
#define BAGL_BAGL_DETAIL_SPARSE_ORDERING_H_

#include <queue>
#include <vector>

#include "bagl/breadth_first_search.h"
#include "bagl/depth_first_search.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

namespace sparse_detail {

// rcm_queue
//
// This is a custom queue type used in the
// *_ordering algorithms.
// In addition to the normal queue operations, the
// rcm_queue provides:
//
//   int eccentricity() const;
//   value_type spouse() const;
//

// yes, it's a bad name...but it works, so use it
template <typename Vertex, typename DegreeMap, typename Container = std::deque<Vertex> >
class rcm_queue : public std::queue<Vertex, Container> {
  using base = std::queue<Vertex, Container>;

 public:
  using value_type = typename base::value_type;
  using size_type = typename base::size_type;

  // SGI queue does not have a contructor queue(const Container&)
  explicit rcm_queue(DegreeMap deg) : degree_(deg) {}

  void pop() {
    if (!size_) {
      q_size_ = base::size();
    }

    base::pop();
    if (size_ == q_size_ - 1) {
      size_ = 0;
      ++eccen_;
    } else {
      ++size_;
    }
  }

  [[nodiscard]] value_type& front() {
    value_type& u = base::front();
    if (size_ == 0) {
      w_ = u;
    } else if (get(degree_, u) < get(degree_, w_)) {
      w_ = u;
    }
    return u;
  }

  [[nodiscard]] const value_type& front() const {
    const value_type& u = base::front();
    if (size_ == 0) {
      w_ = u;
    } else if (get(degree_, u) < get(degree_, w_)) {
      w_ = u;
    }
    return u;
  }

  [[nodiscard]] value_type& top() { return front(); }
  [[nodiscard]] const value_type& top() const { return front(); }

  [[nodiscard]] size_type size() const { return base::size(); }

  [[nodiscard]] size_type eccentricity() const { return eccen_; }
  [[nodiscard]] value_type spouse() const { return w_; }

 protected:
  size_type size_ = 0;
  size_type q_size_ = 1;
  int eccen_ = -1;
  mutable value_type w_;
  DegreeMap degree_;
};

template <typename Tp, typename Container = std::deque<Tp> >
class sparse_ordering_queue : public std::queue<Tp, Container> {
 public:
  decltype(auto) top() const { return this->front(); }
  decltype(auto) top() { return this->front(); }
  auto begin() { return this->c.begin(); }
  auto rbegin() { return this->c.rbegin(); }
  auto end() { return this->c.end(); }
  auto rend() { return this->c.rend(); }
  auto& operator[](int n) { return this->c[n]; }
  auto size() { return this->c.size(); }
};

}  // namespace sparse_detail

// Compute Pseudo peripheral
//
// To compute an approximated peripheral for a given vertex.
// Used in <tt>king_ordering</tt> algorithm.
//
template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          concepts::ReadableVertexPropertyMap<G> DegreeMap>
auto pseudo_peripheral_pair(const G& g, graph_vertex_descriptor_t<G> u, int& ecc, ColorMap color, DegreeMap degree) {
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;

  sparse_detail::rcm_queue<graph_vertex_descriptor_t<G>, DegreeMap> q(degree);
  for (auto v : vertices(g)) {
    if (get(color, v) != Color::red()) {
      put(color, v, Color::white());
    }
  }
  breadth_first_visit(g, u, q, default_bfs_visitor{}, color);

  ecc = q.eccentricity();
  return q.spouse();
}

// Find a good starting node
//
// This is to find a good starting node for the
// king_ordering algorithm. "good" is in the sense
// of the ordering generated by RCM.
//
template <concepts::VertexListGraph G, concepts::ReadWriteVertexPropertyMap<G> ColorMap,
          concepts::ReadableVertexPropertyMap<G> DegreeMap>
auto find_starting_node(const G& g, graph_vertex_descriptor_t<G> r, ColorMap color, DegreeMap degree) {
  int eccen_r = -1;
  int eccen_x = -1;
  auto x = pseudo_peripheral_pair(g, r, eccen_r, color, degree);
  auto y = pseudo_peripheral_pair(g, x, eccen_x, color, degree);

  while (eccen_x > eccen_r) {
    r = x;
    eccen_r = eccen_x;
    x = y;
    y = pseudo_peripheral_pair(g, x, eccen_x, color, degree);
  }
  return x;
}

}  // namespace bagl

#endif  // BAGL_BAGL_DETAIL_SPARSE_ORDERING_H_
