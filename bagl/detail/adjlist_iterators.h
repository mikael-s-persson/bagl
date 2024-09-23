// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_ADJLIST_ITERATORS_H_
#define BAGL_BAGL_ADJLIST_ITERATORS_H_

#include <iterator>
#include <ranges>
#include <utility>
#include <variant>

#include "bagl/detail/container_generators.h"

namespace bagl::adjlist_detail {

template <typename DirectedS>
struct adjlist_traversal_tag : public virtual vertex_list_graph_tag,
                               public virtual incidence_graph_tag,
                               public virtual adjacency_graph_tag,
                               public virtual edge_list_graph_tag {};

template <>
struct adjlist_traversal_tag<undirected_s> : public virtual vertex_list_graph_tag,
                                             public virtual incidence_graph_tag,
                                             public virtual adjacency_graph_tag,
                                             public virtual edge_list_graph_tag,
                                             public virtual bidirectional_graph_tag {};

template <>
struct adjlist_traversal_tag<bidirectional_s> : public virtual vertex_list_graph_tag,
                                                public virtual incidence_graph_tag,
                                                public virtual adjacency_graph_tag,
                                                public virtual edge_list_graph_tag,
                                                public virtual bidirectional_graph_tag {};

template <typename ListS, typename Container>
struct adjlist_select_vertex_range {
  static auto create_range(Container& cont) { return std::views::iota(cont.begin(), cont.end()); }
};

template <typename Container>
struct adjlist_select_vertex_range<vec_s, Container> {
  static auto create_range(Container& cont) { return std::views::iota(0, cont.size()); }
};

template <typename Container>
struct adjlist_select_vertex_range<pool_s, Container> {
  static auto create_range(Container& cont) {
    return std::views::iota(0, cont.m_data.size()) | std::views::filter(container_detail::is_not_hole(&cont.m_data));
  }
};

// This is a concatenation of in-edges and out-edges ranges to make a complete
// range of incident edges for undirected graphs.
template <typename EDesc, typename IERange, typename OERange>
class adjlist_undir_ioerange : public std::ranges::view_interface<adjlist_undir_ioerange<EDesc, IERange, OERange>> {
 public:
  using self_range = adjlist_undir_ioerange<EDesc, IERange, OERange>;
  using ie_iter = decltype(std::declval<IERange>().begin());
  using oe_iter = decltype(std::declval<OERange>().begin());

  struct iterator {
    using value_type = EDesc;
    using reference = const EDesc&;
    using pointer = const EDesc*;
    using difference_type = int;
    using iterator_category = std::bidirectional_iterator_tag;

    iterator& operator++() {
      if (cur_ie != parent->ie_range.end()) {
        ++cur_ie;
      } else if (cur_oe != parent->oe_range.end()) {
        ++cur_oe;
      }
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator& operator--() {
      if (cur_oe == parent->oe_range.begin()) {
        --cur_ie;
      } else {
        --cur_oe;
      }
      return *this;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      --(*this);
      return tmp;
    }
    [[nodiscard]] bool operator==(const iterator& rhs) const {
      return (this->cur_ie == rhs.cur_ie) && (this->cur_oe == rhs.cur_oe);
    }
    [[nodiscard]] bool operator!=(const iterator& rhs) const { return !(*this == rhs); }
    EDesc operator*() const {
      if (cur_ie != parent->ie_range.end()) {
        return EDesc(*cur_ie, parent->as_out_edges);
      }
      return EDesc(*cur_oe, !parent->as_out_edges);
    }
    const EDesc* operator->() const { return &*(*this); }

    ie_iter cur_ie;
    oe_iter cur_oe;
    const self_range* parent;
  };

  iterator begin() const { return {ie_range.begin(), oe_range.begin(), this}; }

  iterator end() const { return {ie_range.end(), oe_range.end(), this}; }

  bool as_out_edges = false;
  IERange ie_range;
  OERange oe_range;
};

template <typename ListS, typename Container, typename EDesc>
struct adjlist_select_out_edge_range {
  template <typename Vertex>
  static auto create_range(Vertex u, Container& cont) {
    return std::views::iota(cont.begin(), cont.end()) |
           std::views::transform([u](const auto& e_it) { return EDesc(u, e_it); });
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container* cont) {
    return create_range(u, *cont);
  }
};

template <typename Container, typename EDesc>
struct adjlist_select_out_edge_range<vec_s, Container, EDesc> {
  template <typename Vertex>
  static auto create_range(Vertex u, Container& cont) {
    return std::views::iota(0, cont.size()) | std::views::transform([u](const auto& e_id) { return EDesc(u, e_id); });
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container* cont) {
    return create_range(u, *cont);
  }
};

template <typename Container, typename EDesc>
struct adjlist_select_out_edge_range<pool_s, Container, EDesc> {
  template <typename Vertex>
  static auto create_range(Vertex u, Container& cont) {
    return std::views::iota(0, cont.m_data.size()) |
           std::views::transform([u](const auto& e_id) { return EDesc(u, e_id); }) |
           std::views::filter(container_detail::edge_is_not_hole(&cont.m_data));
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container* cont) {
    return create_range(u, *cont);
  }
};

// Simply cast bidirectional edge descriptors to undirected ones for edges iteration.
template <typename BidirEdgeRange>
auto adjlist_undir_eiter_range(BidirEdgeRange&& r) {
  return std::views::transform(r, [](const auto& e) { return container_detail::undir_edge_desc(e); });
}

}  // namespace bagl::adjlist_detail

#endif  // BAGL_BAGL_ADJLIST_ITERATORS_H_
