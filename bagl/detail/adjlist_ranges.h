// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_ADJLIST_RANGES_H_
#define BAGL_BAGL_ADJLIST_RANGES_H_

#include <iterator>
#include <ranges>
#include <utility>
#include <variant>

#include "bagl/detail/container_generators.h"
#include "bagl/partial_range.h"

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

template <typename Container>
auto iota_view_for_list_container(Container* cont) {
  using Iter = decltype(cont->begin());
  Iter beg = (cont != nullptr ? cont->begin() : Iter{});
  Iter end = (cont != nullptr ? cont->end() : Iter{});
  return std::views::iota(beg, end);
}

template <typename Container>
auto iota_view_for_container(Container* cont) {
  using SizeT = decltype(cont->size());
  SizeT sz = (cont != nullptr ? cont->size() : SizeT{0});
  return std::views::iota(SizeT{0}, sz);
}

template <typename ListS, typename Container>
struct adjlist_select_vertex_range {
  static auto create_range(Container& cont) { return iota_view_for_list_container(&cont); }
};

template <typename Container>
struct adjlist_select_vertex_range<vec_s, Container> {
  static auto create_range(Container* cont) { return iota_view_for_container(cont); }
  static auto create_range(Container& cont) { return create_range(&cont); }
  static auto create_empty_range() { return create_range(nullptr); }
};

template <typename Container>
struct adjlist_select_vertex_range<pool_s, Container> {
  static auto create_range_impl(decltype(Container::m_data)* cont_data) {
    return iota_view_for_container(cont_data) | std::views::filter(container_detail::is_not_hole(cont_data));
  }
  static auto create_range(Container& cont) { return create_range_impl(&cont.m_data); }
  static auto create_range(Container* cont) { return create_range_impl(cont != nullptr ? &cont->m_data : nullptr); }
  static auto create_empty_range() { return create_range_impl(nullptr); }
};

// This is a concatenation of in-edges and out-edges ranges to make a complete
// range of incident edges for undirected graphs.
template <typename EDesc, typename IERange, typename OERange>
class adjlist_undir_ioerange : public std::ranges::view_interface<adjlist_undir_ioerange<EDesc, IERange, OERange>> {
 public:
  using self_range = adjlist_undir_ioerange<EDesc, IERange, OERange>;
  using ie_iter = std::decay_t<decltype(std::declval<IERange>().begin())>;
  using oe_iter = std::decay_t<decltype(std::declval<OERange>().begin())>;

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

    iterator(ie_iter&& a_cur_ie, oe_iter&& a_cur_oe, const self_range* a_parent)
        : cur_ie(std::move(a_cur_ie)), cur_oe(std::move(a_cur_oe)), parent(a_parent) {}
    iterator() = default;

    ie_iter cur_ie;
    oe_iter cur_oe;
    const self_range* parent;
  };

  iterator begin() const { return iterator{ie_range.begin(), oe_range.begin(), this}; }
  iterator end() const { return iterator{ie_range.end(), oe_range.end(), this}; }

  adjlist_undir_ioerange() = default;
  adjlist_undir_ioerange(bool a_as_out_edges, IERange a_ie_range, OERange a_oe_range)
      : as_out_edges(a_as_out_edges), ie_range(a_ie_range), oe_range(a_oe_range) {}

  bool as_out_edges = false;
  // This has to be mutable because standard ranges and views are allowed to mutate.
  mutable IERange ie_range;
  mutable OERange oe_range;
};

template <typename EDesc, typename IERange, typename OERange>
auto make_adjlist_undir_ioerange(bool as_out_edges, IERange ie_range, OERange oe_range) {
  return adjlist_undir_ioerange<EDesc, IERange, OERange>(as_out_edges, ie_range, oe_range);
}

template <typename ListS, typename Container, typename EDesc>
struct adjlist_select_out_edge_range {
  template <typename Vertex>
  static auto create_range(Vertex u, Container* cont) {
    return iota_view_for_list_container(cont) | std::views::transform([u, cont](const auto& e_it) {
             return EDesc(u, container_detail::iterator_to_desc(*cont, e_it));
           });
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container& cont) {
    return create_range(u, &cont);
  }
  template <typename Vertex>
  static auto create_empty_range(Vertex u) {
    return create_range(u, nullptr);
  }
};

template <typename Container, typename EDesc>
struct adjlist_select_out_edge_range<vec_s, Container, EDesc> {
  template <typename Vertex>
  static auto create_range(Vertex u, Container* cont) {
    return iota_view_for_container(cont) | std::views::transform([u](const auto& e_id) { return EDesc(u, e_id); });
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container& cont) {
    return create_range(u, &cont);
  }
  template <typename Vertex>
  static auto create_empty_range(Vertex u) {
    return create_range(u, nullptr);
  }
};

template <typename Container, typename EDesc>
struct adjlist_select_out_edge_range<pool_s, Container, EDesc> {
  template <typename Vertex>
  static auto create_range_impl(Vertex u, decltype(Container::m_data)* cont_data) {
    return iota_view_for_container(cont_data) |
           std::views::transform([u](const auto& e_id) { return EDesc(u, e_id); }) |
           std::views::filter(container_detail::edge_is_not_hole(cont_data));
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container* cont) {
    return create_range_impl(u, cont != nullptr ? &cont->m_data : nullptr);
  }
  template <typename Vertex>
  static auto create_range(Vertex u, Container& cont) {
    return create_range_impl(u, &cont.m_data);
  }
  template <typename Vertex>
  static auto create_empty_range(Vertex u) {
    return create_range_impl(u, nullptr);
  }
};

// Simply cast bidirectional edge descriptors to undirected ones for edges iteration.
template <typename BidirEdgeRange>
auto adjlist_undir_eiter_range(BidirEdgeRange&& r) {
  auto e_to_ue = [](const auto& e) { return container_detail::undir_edge_desc(e); };
  return std::ranges::transform_view<BidirEdgeRange, decltype(e_to_ue)>(std::forward<BidirEdgeRange>(r), e_to_ue);
}

}  // namespace bagl::adjlist_detail

#endif  // BAGL_BAGL_ADJLIST_RANGES_H_
