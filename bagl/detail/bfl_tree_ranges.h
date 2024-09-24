// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_DETAIL_BFL_TREE_RANGES_H_
#define BAGL_BAGL_DETAIL_BFL_TREE_RANGES_H_

#include <iterator>
#include <limits>
#include <ranges>
#include <type_traits>
#include <variant>

#include "bagl/detail/container_generators.h"

namespace bagl::bfl_detail {

// This library provides detail-implementations of BFL-tree ranges.

template <typename VertexProperties, bool isEmpty>
struct bfltree_vertex_value_base : VertexProperties {
  const VertexProperties& vertex() const { return *this; }
  VertexProperties& vertex() { return *this; }
};

template <typename VertexProperties>
struct bfltree_vertex_value_base<VertexProperties, false> {
  VertexProperties v_data;
  const VertexProperties& vertex() const { return this->v_data; }
  VertexProperties& vertex() { return this->v_data; }
};

template <typename EdgeProperties, bool isEmpty>
struct bfltree_edge_value_base : EdgeProperties {
  const EdgeProperties& edge() const { return *this; }
  EdgeProperties& edge() { return *this; }
};

template <typename EdgeProperties>
struct bfltree_edge_value_base<EdgeProperties, false> {
  EdgeProperties e_data;
  const EdgeProperties& edge() const { return this->e_data; }
  EdgeProperties& edge() { return this->e_data; }
};

// uses the empty base-class optimization, if possible:
template <typename VertexProperties, typename EdgeProperties>
struct bfltree_value_type
    : bfltree_vertex_value_base<VertexProperties,
                                std::is_empty_v<VertexProperties>>,
      bfltree_edge_value_base<EdgeProperties, std::is_empty_v<EdgeProperties>> {
  std::size_t out_degree;

  bfltree_value_type()
      : bfltree_vertex_value_base<VertexProperties,
                                  std::is_empty_v<VertexProperties>>(),
        bfltree_edge_value_base<EdgeProperties,
                                std::is_empty_v<EdgeProperties>>(),
        out_degree((std::numeric_limits<std::size_t>::max)()) {}
};

struct bfltree_edge_desc {
  std::size_t target_vertex;
  explicit bfltree_edge_desc(
      std::size_t aTarget = (std::numeric_limits<std::size_t>::max)())
      : target_vertex(aTarget) {}
};
inline bool operator==(const bfltree_edge_desc& lhs,
                       const bfltree_edge_desc& rhs) {
  return (lhs.target_vertex == rhs.target_vertex);
}
inline bool operator!=(const bfltree_edge_desc& lhs,
                       const bfltree_edge_desc& rhs) {
  return (lhs.target_vertex != rhs.target_vertex);
}
inline bool operator<(const bfltree_edge_desc& lhs,
                      const bfltree_edge_desc& rhs) {
  return (lhs.target_vertex < rhs.target_vertex);
}
inline bool operator<=(const bfltree_edge_desc& lhs,
                       const bfltree_edge_desc& rhs) {
  return (lhs.target_vertex <= rhs.target_vertex);
}
inline bool operator>(const bfltree_edge_desc& lhs,
                      const bfltree_edge_desc& rhs) {
  return (lhs.target_vertex > rhs.target_vertex);
}
inline bool operator>=(const bfltree_edge_desc& lhs,
                       const bfltree_edge_desc& rhs) {
  return (lhs.target_vertex >= rhs.target_vertex);
}
struct bfltree_edge_desc_from_target {
  bfltree_edge_desc operator()(std::size_t v) const { return bfltree_edge_desc{v}; }
};

template <typename VProp>
bool bfltree_is_vertex_valid(const VProp& vp) {
  return (vp.out_degree != (std::numeric_limits<std::size_t>::max)());
}

template <typename Container>
struct bfltree_vertex_validity {
  const Container* p_cont;
  explicit bfltree_vertex_validity(const Container* aPCont = nullptr)
      : p_cont(aPCont) {}
  bool operator()(std::size_t d) {
    return ((d < p_cont->size()) && (bfltree_is_vertex_valid((*p_cont)[d])));
  }
};

template <typename Container>
struct bfltree_edge_validity {
  const Container* p_cont;
  explicit bfltree_edge_validity(const Container* aPCont = nullptr)
      : p_cont(aPCont) {}
  bool operator()(bfltree_edge_desc d) {
    return ((d.target_vertex < p_cont->size()) &&
            (bfltree_is_vertex_valid((*p_cont)[d.target_vertex])));
  }
};
#if 0
struct bfltree_eiter
    : public iterator_facade<bfltree_eiter, const bfltree_edge_desc,
                             std::random_access_iterator_tag> {
 public:
  using self = bfltree_eiter;

  explicit bfltree_eiter(bfltree_edge_desc aE = bfltree_edge_desc()) : e(aE) {}

  // private:
  friend class iterator_core_access;

  void increment() { ++e.target_vertex; }
  void decrement() { --e.target_vertex; }
  bool equal(const self& rhs) const { return (this->e == rhs.e); }
  const bfltree_edge_desc& dereference() const { return e; }

  void advance(std::ptrdiff_t i) { e.target_vertex += i; }
  std::ptrdiff_t distance_to(const self& rhs) const {
    return rhs.e.target_vertex - this->e.target_vertex;
  }

  bfltree_edge_desc e;
};

struct bfltree_viter : public iterator_facade<bfltree_viter, const std::size_t,
                                              std::random_access_iterator_tag> {
 public:
  using self = bfltree_viter;

  explicit bfltree_viter(std::size_t aVIt = 0) : v_it(aVIt) {}

  // private:
  friend class iterator_core_access;

  void increment() { ++v_it; }
  void decrement() { --v_it; }
  bool equal(const self& rhs) const { return (this->v_it == rhs.v_it); }
  const std::size_t& dereference() const { return v_it; }

  void advance(std::ptrdiff_t i) { v_it += i; }
  std::ptrdiff_t distance_to(const self& rhs) const {
    return rhs.v_it - this->v_it;
  }

  std::size_t v_it;
};
#endif
}  // namespace bagl::bfl_detail

#endif // BAGL_BAGL_DETAIL_BFL_TREE_RANGES_H_
