// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_DETAIL_VEBL_TREE_RANGES_H_
#define BAGL_BAGL_DETAIL_VEBL_TREE_RANGES_H_

#include <algorithm>
#include <iterator>
#include <limits>
#include <vector>

#include "bagl/detail/bfl_tree_ranges.h"

namespace bagl::bfl_detail {

// This library provides detail-implementations of vEBL-tree ranges.

struct vebl_depth_records {
  std::vector<std::size_t> T;
  std::vector<std::size_t> B;
  std::vector<std::size_t> D;
};

template <std::size_t Base>
struct s_power_array {
  static constexpr auto impl() noexcept {
    constexpr std::size_t max_power =
        static_cast<std::size_t>(std::numeric_limits<std::size_t>::digits * std::numeric_limits<std::size_t>::radix) /
        Base;
    std::array<std::size_t, max_power> vals{};
    vals[0] = 1;
    for (std::size_t i = 1; i < max_power; ++i) {
      vals[i] = vals[i - 1] * Base;
    }
    return vals;
  }

  // Give the joined string static storage
  static constexpr auto values = impl();
};

template <std::size_t Base>
constexpr std::size_t s_power(std::size_t aPower) {
  return s_power_array<Base>::values[aPower];
}

template <std::size_t Arity>
struct s_treesize_array {
  static constexpr auto impl() noexcept {
    constexpr std::size_t max_depth =
        static_cast<std::size_t>(std::numeric_limits<std::size_t>::digits * std::numeric_limits<std::size_t>::radix) /
        Arity;
    std::array<std::size_t, max_depth> vals{};
    vals[0] = 1;
    for (std::size_t i = 1; i < max_depth; ++i) {
      vals[i] = vals[i - 1] * Arity + 1;
    }
    return vals;
  }

  // Give the joined string static storage
  static constexpr auto values = impl();
};

template <std::size_t Arity>
constexpr std::size_t s_treesize(std::size_t aDepth) {
  return s_treesize_array<Arity>::values[aDepth];
}

template <std::size_t Arity>
std::size_t get_bfl_index_depth(std::size_t bfl_i) {
  const auto& s = s_treesize_array<Arity>::values;
  return std::upper_bound(s.begin(), s.end(), bfl_i) - s.begin();
}

template <std::size_t Arity>
std::size_t convert_bfl_to_vebl(std::size_t bfl_i, const vebl_depth_records& rec) {
  std::size_t d = get_bfl_index_depth<Arity>(bfl_i);

  std::size_t vebl_i = 0;
  while (bfl_i > 0) {
    std::size_t d_diff = d - rec.D[d];
    vebl_i += rec.T[d] + rec.B[d] * ((bfl_i + 1) % (s_power<Arity>(d) / s_power<Arity>(rec.D[d])));
    for (std::size_t i = d_diff; i > 0; --i) {
      bfl_i = (bfl_i - 1) / Arity;
    }
    d -= d_diff;
  }

  return vebl_i;
}

template <std::size_t Arity>
void extend_vebl_depth_records(vebl_depth_records& rec) {
  std::size_t cur_d = rec.T.size() - 1;

  // find the first unequal T vs B heights, starting from the bottom layer:
  std::size_t uneq_d = cur_d;
  while ((uneq_d > 0) && (rec.T[uneq_d] <= rec.B[uneq_d])) {
    uneq_d = rec.D[uneq_d];
  }

  if (uneq_d == 0) {
    // this means we did not find an unequal T-B pair.
    // we must create a new top-level layer.
    rec.T.push_back(s_treesize<Arity>(cur_d));
    rec.B.push_back(1);
    rec.D.push_back(0);
    if (rec.B[0] == 1) {
      rec.B[0] = Arity + 1;
    }
  } else {
    rec.T.push_back(s_treesize<Arity>(cur_d - uneq_d));
    // rec.T.push_back(rec.B[uneq_d]);
    rec.B.push_back(1);
    rec.D.push_back(uneq_d);

    while ((uneq_d > 0) && (rec.T[uneq_d] > rec.B[uneq_d])) {
      rec.B[uneq_d] = rec.B[uneq_d] * Arity + 1;
      // HERE IS WHERE COPYING THE SUBTREE IS NEEDED!
      uneq_d = rec.D[uneq_d];
    }
  }
}

template <std::size_t Arity, typename ValueType, typename IdIter>
void rec_copy_vebl_storage(std::vector<ValueType>& cont, vebl_depth_records& rec, IdIter it, IdIter it_end,
                           std::size_t base_i) {
  if (it == it_end) {
    return;
  }
  std::size_t next_b = rec.B[*it] * Arity + 1;
  base_i += rec.T[*it];
  for (std::ptrdiff_t j = (s_power<Arity>(*it) / s_power<Arity>(rec.D[*it])) - 1; j >= 0; --j) {
    if (j > 0) {
      std::copy(std::make_move_iterator(cont.begin() + (base_i + j * rec.B[*it])),
                std::make_move_iterator(cont.begin() + (base_i + (j + 1) * rec.B[*it])),
                cont.begin() + (base_i + j * next_b));
    }
    rec_copy_vebl_storage<Arity>(cont, rec, it + 1, it_end, base_i + j * next_b);
  }
}

template <std::size_t Arity, typename ValueType>
void update_and_copy_vebl_storage(std::vector<ValueType>& cont, vebl_depth_records& rec, std::size_t uneq_d) {
  if ((uneq_d == 0) || (rec.T[uneq_d] <= rec.B[uneq_d])) {
    return;
  }

  std::vector<std::size_t> affected_depths;
  while ((uneq_d > 0) && (rec.T[uneq_d] > rec.B[uneq_d])) {
    affected_depths.push_back(uneq_d);
    uneq_d = rec.D[uneq_d];
  }

  rec_copy_vebl_storage<Arity>(cont, rec, affected_depths.rbegin(), affected_depths.rend(), 0);

  for (std::size_t& affected_depth : affected_depths) {
    rec.B[affected_depth] = rec.B[affected_depth] * Arity + 1;
  }
}

template <std::size_t Arity, typename ValueType>
void extend_vebl_storage(std::vector<ValueType>& cont, vebl_depth_records& rec) {
  std::size_t cur_d = rec.T.size() - 1;

  // find the first unequal T vs B heights, starting from the bottom layer:
  std::size_t uneq_d = cur_d;
  while ((uneq_d > 0) && (rec.T[uneq_d] <= rec.B[uneq_d])) {
    uneq_d = rec.D[uneq_d];
  }

  if (uneq_d == 0) {
    // this means we did not find an unequal T-B pair.
    // we must create a new top-level layer.
    rec.T.push_back(s_treesize<Arity>(cur_d));
    rec.B.push_back(1);
    rec.D.push_back(0);
    if (rec.B[0] == 1) {
      rec.B[0] = Arity + 1;
    }

    std::size_t bottom_layer = s_power<Arity>(cur_d + 1);
    cont.resize(cont.size() + bottom_layer);

  } else {
    rec.T.push_back(s_treesize<Arity>(cur_d - uneq_d));
    // rec.T.push_back(rec.B[uneq_d]);
    rec.B.push_back(1);
    rec.D.push_back(uneq_d);

    std::size_t bottom_layer = s_power<Arity>(cur_d + 1);
    cont.resize(cont.size() + bottom_layer);

    update_and_copy_vebl_storage<Arity>(cont, rec, uneq_d);
  }
}

template <typename Container, std::size_t Arity>
struct vebltree_vertex_validity {
  const Container* p_cont;
  const vebl_depth_records* p_drec;
  explicit vebltree_vertex_validity(const Container* aPCont = nullptr, const vebl_depth_records* aPDRec = nullptr)
      : p_cont(aPCont), p_drec(aPDRec) {}
  bool operator()(std::size_t d) {
    std::size_t d_vebl = convert_bfl_to_vebl<Arity>(d, *p_drec);
    return ((d_vebl < p_cont->size()) && bfltree_is_vertex_valid((*p_cont)[d_vebl]));
  }
};

}  // namespace bagl::bfl_detail

#endif  // BAGL_BAGL_DETAIL_VEBL_TREE_RANGES_H_
