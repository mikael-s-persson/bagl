// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_EDGE_LIST_H_
#define BAGL_BAGL_EDGE_LIST_H_

#include <iterator>
#include <ranges>
#include <type_traits>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"

namespace bagl {

//
// The edge_list class is an EdgeListGraph module that is constructed
// from a pair of iterators whose value type is a pair of vertex
// descriptors.
//
// For example:
//
//  typedef std::pair<int,int> E;
//  list<E> elist;
//  ...
//  typedef edge_list<list<E>::iterator> Graph;
//  Graph g(elist.begin(), elist.end());
//
// If the iterators are random access, then Graph::edge_descriptor
// is of Integral type, otherwise it is a struct, though it is
// convertible to an Integral type.
//

struct edge_list_tag {};

// The implementation class for edge_list.
template <class G, class EdgeIter, class T, class D>
class edge_list_impl {
 public:
  using edge_id = D;
  using Vpair = T;
  using V = typename Vpair::first_type;
  using vertex_descriptor = V;
  using graph_tag = edge_list_tag;
  using edge_property_type = void;

  struct edge_descriptor {
    edge_descriptor() = default;
    edge_descriptor(EdgeIter p, edge_id id) : _ptr(p), _id(id) {}
    explicit operator edge_id() { return _id; }
    EdgeIter _ptr;
    edge_id _id;
  };
  using E = edge_descriptor;

  struct edge_iterator {
    using self = edge_iterator;
    using value_type = E;
    using reference = E&;
    using pointer = E*;
    using difference_type = std::ptrdiff_t;
    using iterator_category = std::input_iterator_tag;
    edge_iterator() = default;
    explicit edge_iterator(EdgeIter iter) : _iter(iter), _i(0) {}
    E operator*() { return E(_iter, _i); }
    self& operator++() {
      ++_iter;
      ++_i;
      return *this;
    }
    self operator++(int) {
      self t = *this;
      ++(*this);
      return t;
    }
    bool operator==(const self& x) { return _iter == x._iter; }
    bool operator!=(const self& x) { return _iter != x._iter; }
    EdgeIter _iter;
    edge_id _i;
  };
};

template <class G, class EI, class T, class D>
auto edges(const edge_list_impl<G, EI, T, D>& g_) {
  const G& g = static_cast<const G&>(g_);
  using edge_iterator = typename edge_list_impl<G, EI, T, D>::edge_iterator;
  return std::ranges::subrange(edge_iterator(g.first_), edge_iterator(g.last_));
}
template <class G, class EI, class T, class D>
auto source(typename edge_list_impl<G, EI, T, D>::edge_descriptor e, const edge_list_impl<G, EI, T, D>& /*unused*/) {
  return (*e._ptr).first;
}
template <class G, class EI, class T, class D>
auto target(typename edge_list_impl<G, EI, T, D>::edge_descriptor e, const edge_list_impl<G, EI, T, D>& /*unused*/) {
  return (*e._ptr).second;
}

template <class D, class E>
class el_edge_property_map : public put_get_helper<el_edge_property_map<D, E> > {
 public:
  using key_type = E;
  using value_type = D;
  using reference = D;

  value_type operator[](key_type e) const { return e._i; }
};
struct edge_list_edge_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using type = el_edge_property_map<typename Graph::edge_id, typename Graph::edge_descriptor>;
    using const_type = type;
  };
};
template <>
struct edge_property_selector<edge_list_tag> {
  using type = edge_list_edge_property_selector;
};

template <class G, class EI, class T, class D>
auto get(edge_index_t /*unused*/, const edge_list_impl<G, EI, T, D>& /*unused*/) {
  using EdgeIndexMap = property_map_t<edge_list_impl<G, EI, T, D>, edge_index_t>;
  return EdgeIndexMap();
}

template <class G, class EI, class T, class D>
D get(edge_index_t /*unused*/, const edge_list_impl<G, EI, T, D>& /*unused*/,
      typename edge_list_impl<G, EI, T, D>::edge_descriptor e) {
  return e._i;
}

// A specialized implementation for when the iterators are random access.

struct edge_list_ra_tag {};

template <class G, class EdgeIter, class T, class D>
class edge_list_impl_ra {
 public:
  using edge_id = D;
  using Vpair = T;
  using V = typename Vpair::first_type;
  using graph_tag = edge_list_ra_tag;
  using edge_property_type = void;

  using edge_descriptor = edge_id;
  using vertex_descriptor = V;
};

template <class G, class EI, class T, class D>
auto edges(const edge_list_impl_ra<G, EI, T, D>& g_) {
  const G& g = static_cast<const G&>(g_);
  return std::ranges::iota_view(D{0}, D{g.last_ - g.first_});
}
template <class G, class EI, class T, class D>
auto source(typename edge_list_impl_ra<G, EI, T, D>::edge_descriptor e, const edge_list_impl_ra<G, EI, T, D>& g_) {
  const G& g = static_cast<const G&>(g_);
  return g.first_[e].first;
}
template <class G, class EI, class T, class D>
auto target(typename edge_list_impl_ra<G, EI, T, D>::edge_descriptor e, const edge_list_impl_ra<G, EI, T, D>& g_) {
  const G& g = static_cast<const G&>(g_);
  return g.first_[e].second;
}
template <class E>
class el_ra_edge_property_map : public put_get_helper<el_ra_edge_property_map<E> > {
 public:
  using key_type = E;
  using value_type = E;
  using reference = E;

  value_type operator[](key_type e) const { return e; }
};
struct edge_list_ra_edge_property_selector {
  template <class Graph, class Property, class Tag>
  struct bind_ {
    using type = el_ra_edge_property_map<typename Graph::edge_descriptor>;
    using const_type = type;
  };
};
template <>
struct edge_property_selector<edge_list_ra_tag> {
  using type = edge_list_ra_edge_property_selector;
};
template <class G, class EI, class T, class D>
auto get(edge_index_t /*unused*/, const edge_list_impl_ra<G, EI, T, D>& /*unused*/) {
  using EdgeIndexMap = property_map_t<edge_list_impl_ra<G, EI, T, D>, edge_index_t>;
  return EdgeIndexMap();
}

template <class G, class EI, class T, class D>
D get(edge_index_t /*unused*/, const edge_list_impl_ra<G, EI, T, D>& /*unused*/,
      typename edge_list_impl_ra<G, EI, T, D>::edge_descriptor e) {
  return e;
}

// The edge_list class conditionally inherits from one of the
// above two classes.

template <class EdgeIter, class T = typename std::iterator_traits<EdgeIter>::value_type,
          class D = typename std::iterator_traits<EdgeIter>::difference_type,
          class Cat = typename std::iterator_traits<EdgeIter>::iterator_category>
class edge_list : public std::conditional_t<std::is_convertible_v<Cat, std::random_access_iterator_tag>,
                                            edge_list_impl_ra<edge_list<EdgeIter, T, D, Cat>, EdgeIter, T, D>,
                                            edge_list_impl<edge_list<EdgeIter, T, D, Cat>, EdgeIter, T, D> > {
 public:
  using directed_category = directed_tag;
  using edge_parallel_category = allow_parallel_edge_tag;
  using traversal_category = edge_list_graph_tag;
  edge_list(EdgeIter first, EdgeIter last) : first_(first), last_(last), num_edges_(std::distance(first, last)) {}
  edge_list(EdgeIter first, EdgeIter last, std::size_t E) : first_(first), last_(last), num_edges_(E) {}

  EdgeIter first_;
  EdgeIter last_;
  std::size_t num_edges_;
};

template <class EdgeIter, class T, class D, class Cat>
std::size_t num_edges(const edge_list<EdgeIter, T, D, Cat>& el) {
  return el.num_edges_;
}

template <class EdgeIter>
edge_list<EdgeIter> make_edge_list(EdgeIter first, EdgeIter last) {
  return edge_list<EdgeIter>(first, last);
}

}  // namespace bagl

#endif  // BAGL_BAGL_EDGE_LIST_H_
