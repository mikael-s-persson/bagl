// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_DETAIL_BFL_TREE_RANGES_H_
#define BAGL_BAGL_DETAIL_BFL_TREE_RANGES_H_

#include <compare>
#include <limits>
#include <type_traits>

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
struct bfltree_value_type : bfltree_vertex_value_base<VertexProperties, std::is_empty_v<VertexProperties>>,
                            bfltree_edge_value_base<EdgeProperties, std::is_empty_v<EdgeProperties>> {
  std::size_t out_degree;

  bfltree_value_type()
      : bfltree_vertex_value_base<VertexProperties, std::is_empty_v<VertexProperties>>(),
        bfltree_edge_value_base<EdgeProperties, std::is_empty_v<EdgeProperties>>(),
        out_degree((std::numeric_limits<std::size_t>::max)()) {}
};

struct bfltree_edge_desc {
  std::size_t target_vertex;
  explicit bfltree_edge_desc(std::size_t aTarget = (std::numeric_limits<std::size_t>::max)())
      : target_vertex(aTarget) {}
  [[nodiscard]] constexpr auto operator<=>(const bfltree_edge_desc& rhs) const {
    return target_vertex <=> rhs.target_vertex;
  }
};

template <typename VProp>
bool bfltree_is_vertex_valid(const VProp& vp) {
  return (vp.out_degree != (std::numeric_limits<std::size_t>::max)());
}

template <typename Container>
struct bfltree_vertex_validity {
  const Container* p_cont;
  explicit bfltree_vertex_validity(const Container* aPCont = nullptr) : p_cont(aPCont) {}
  bool operator()(std::size_t d) { return ((d < p_cont->size()) && (bfltree_is_vertex_valid((*p_cont)[d]))); }
};

}  // namespace bagl::bfl_detail

#endif  // BAGL_BAGL_DETAIL_BFL_TREE_RANGES_H_
