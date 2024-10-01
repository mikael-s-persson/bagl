#ifndef BAGL_BAGL_EDGES_FROM_OUT_EDGES_H_
#define BAGL_BAGL_EDGES_FROM_OUT_EDGES_H_

#include <iterator>
#include <ranges>
#include <type_traits>
#include <utility>
#include <variant>

#include "bagl/partial_range.h"

namespace bagl {

// This is a concatenation of out-edges ranges from vertices range to make a complete range of edges.
// It could be replaced by the following if it worked (join_view is not implemented):
//   vertices(g) | std::views::transform([&g](vertex_descriptor u) { return out_edges(u, g); }) | std::views::join;
template <typename G>
class edges_from_out_edges : public std::ranges::view_interface<edges_from_out_edges<G>> {
 public:
  using self_range = edges_from_out_edges<G>;
  using v_range_t = decltype(vertices(std::declval<G>()));
  using v_iter_t = std::ranges::iterator_t<const v_range_t>;
  using oe_range_t = decltype(out_edges(*std::declval<v_iter_t>(), std::declval<G>()));
  using oe_iter_t = std::ranges::iterator_t<const oe_range_t>;
  using e_desc_t = decltype(*std::declval<oe_iter_t>());

  explicit edges_from_out_edges(const G& g) : p_g(&g), v_range(vertices(*p_g)) {}
  explicit edges_from_out_edges(const G&& g) = delete;
  edges_from_out_edges() = default;

  struct iterator {
    using value_type = e_desc_t;
    using reference = const e_desc_t&;
    using pointer = const e_desc_t*;
    using difference_type = int;
    using iterator_category = std::bidirectional_iterator_tag;

    iterator& operator++() {
      oe_move_to_next();
      while (oe_is_empty() && cur_vi != parent->v_range.end()) {
        cur_oe = parent->get_fresh_eo_range(++cur_vi);
      }
      return *this;
    }
    iterator operator++(int) {
      iterator tmp = *this;
      ++(*this);
      return tmp;
    }
    iterator& operator--() {
      while (oe_is_at_begin() && cur_vi != parent->v_range.begin()) {
        cur_oe = parent->get_fresh_eo_range(--cur_vi);
        oe_move_to_end();
      }
      // If cur_vi is begin and this iterator is also begin, decrementing is undefined behavior.
      oe_move_to_prev();
      return *this;
    }
    iterator operator--(int) {
      iterator tmp = *this;
      --(*this);
      return tmp;
    }
    [[nodiscard]] bool operator==(const iterator& rhs) const {
      return (this->cur_vi == rhs.cur_vi) && ((this->oe_is_empty() && rhs.oe_is_empty()) ||
                                              ((!this->oe_is_empty() && !rhs.oe_is_empty()) &&
                                               std::get<0>(this->cur_oe).begin() == std::get<0>(rhs.cur_oe).begin()));
    }
    [[nodiscard]] bool operator!=(const iterator& rhs) const { return !(*this == rhs); }
    e_desc_t operator*() const { return *std::get<0>(cur_oe).begin(); }
    const e_desc_t* operator->() const { return &*(*this); }

    iterator(const self_range* p, std::true_type /*begin*/)
        : cur_vi(p->v_range.begin()), cur_oe(p->get_fresh_eo_range(cur_vi)), parent(p) {
      // Seek first valid iterator, or end.
      while (oe_is_empty() && cur_vi != parent->v_range.end()) {
        cur_oe = parent->get_fresh_eo_range(++cur_vi);
      }
    }
    iterator(const self_range* p, std::false_type /*not begin*/)
        : cur_vi(p->v_range.end()), cur_oe(std::ranges::empty_view<e_desc_t>()), parent(p) {}

    iterator() : cur_oe(std::ranges::empty_view<e_desc_t>()) {}

    [[nodiscard]] bool oe_is_at_begin() const {
      return (cur_oe.index() == 0 ? (std::get<0>(cur_oe).begin() == std::get<0>(cur_oe).base_begin()) : true);
    }
    [[nodiscard]] bool oe_is_empty() const { return (cur_oe.index() == 0 ? std::get<0>(cur_oe).empty() : true); }
    void oe_move_to_end() {
      if (cur_oe.index() == 0) {
        std::get<0>(cur_oe).move_begin_to_end();
      }
    }
    void oe_move_to_next() {
      if (cur_oe.index() == 0) {
        std::get<0>(cur_oe).move_to_next();
      }
    }
    void oe_move_to_prev() {
      if (cur_oe.index() == 0) {
        std::get<0>(cur_oe).move_to_prev();
      }
    }

    v_iter_t cur_vi;
    std::variant<partial_view<oe_range_t>, std::ranges::empty_view<e_desc_t>> cur_oe;
    const self_range* parent = nullptr;
  };

  std::variant<partial_view<oe_range_t>, std::ranges::empty_view<e_desc_t>> get_fresh_eo_range(v_iter_t vi) const {
    if (vi == v_range.end()) {
      return std::ranges::empty_view<e_desc_t>();  // Empty range.
    }
    return partial_view<oe_range_t>{out_edges(*vi, *p_g)};
  }

  iterator begin() const { return iterator{this, std::true_type{}}; }
  iterator begin() { return iterator{this, std::true_type{}}; }
  iterator end() const { return iterator{this, std::false_type{}}; }
  iterator end() { return iterator{this, std::false_type{}}; }

  const G* p_g = nullptr;
  v_range_t v_range;
};

}  // namespace bagl

#endif  // BAGL_BAGL_EDGES_FROM_OUT_EDGES_H_