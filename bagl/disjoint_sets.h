// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DISJOINT_SETS_H_
#define BAGL_BAGL_DISJOINT_SETS_H_

#include <vector>

#include "bagl/detail/disjoint_sets.h"
#include "bagl/properties.h"

namespace bagl {

struct find_with_path_halving {
  template <typename ParentPA>
  auto operator()(ParentPA p, property_traits_value_t<ParentPA> v) {
    return disjoint_sets_detail::find_representative_with_path_halving(p, v);
  }
};

struct find_with_full_path_compression {
  template <typename ParentPA>
  auto operator()(ParentPA p, property_traits_value_t<ParentPA> v) {
    return disjoint_sets_detail::find_representative_with_full_compression(p, v);
  }
};

// This is a generalized functor to provide disjoint sets operations
// with "union by rank" and "path compression".  A disjoint-set data
// structure maintains a collection S={S1, S2, ..., Sk} of disjoint
// sets. Each set is identified by a representative, which is some
// member of of the set. Sets are represented by rooted trees. Two
// heuristics: "union by rank" and "path compression" are used to
// speed up the operations.

// Disjoint Set requires two vertex properties for internal use.  A
// RankPA and a ParentPA. The RankPA must map Vertex to some Integral type
// (preferably the size_type associated with Vertex). The ParentPA
// must map Vertex to Vertex.
template <typename RankPA, typename ParentPA, typename FindCompress = find_with_full_path_compression>
class disjoint_sets {
  using self = disjoint_sets;

 public:
  using Vertex = property_traits_value_t<ParentPA>;

  disjoint_sets() = delete;

  disjoint_sets(RankPA r, ParentPA p) : rank_(r), parent_(p) {}

  // Make Set -- Create a singleton set containing vertex x
  void make_set(Vertex x) {
    put(parent_, x, x);
    using R = property_traits_value_t<RankPA>;
    put(rank_, x, R());
  }

  // Link - union the two sets represented by vertex x and y
  void link(Vertex x, Vertex y) { disjoint_sets_detail::link_sets(parent_, rank_, x, y); }

  // Union-Set - union the two sets containing vertex x and y
  void union_set(Vertex x, Vertex y) { link(find_set(x), find_set(y)); }

  // Find-Set - returns the Element representative of the set
  // containing Element x and applies path compression.
  Vertex find_set(Vertex x) { return rep_(parent_, x); }

  template <std::ranges::input_range VRange>
  std::size_t count_sets(VRange& v_rg) {
    std::size_t count = 0;
    for (Vertex v : v_rg) {
      if (get(parent_, v) == v) {
        ++count;
      }
    }
    return count;
  }

  template <std::ranges::input_range VRange>
  void normalize_sets(VRange& v_rg) {
    for (Vertex v : v_rg) {
      disjoint_sets_detail::normalize_node(parent_, v);
    }
  }

  template <std::ranges::input_range VRange>
  void compress_sets(VRange& v_rg) {
    for (Vertex v : v_rg) {
      disjoint_sets_detail::find_representative_with_full_compression(parent_, v);
    }
  }

 protected:
  RankPA rank_;
  ParentPA parent_;
  FindCompress rep_;
};

template <class ID = identity_property_map, class InverseID = identity_property_map,
          class FindCompress = find_with_full_path_compression>
class disjoint_sets_with_storage {
  using Index = property_traits_value_t<ID>;
  using ParentContainer = std::vector<Index>;
  using RankContainer = std::vector<unsigned char>;

 public:
  using size_type = typename ParentContainer::size_type;

  explicit disjoint_sets_with_storage(size_type n = 0, ID id = ID(), InverseID inv = InverseID())
      : id_(id), id_to_vertex_(inv), rank_(n, 0), parent_(n) {
    for (Index i = 0; i < n; ++i) {
      parent_[i] = i;
    }
  }
  // note this is not normally needed
  template <class Element>
  void make_set(Element x) {
    parent_[x] = x;
    rank_[x] = 0;
  }
  template <class Element>
  void link(Element x, Element y) {
    extend_sets(x, y);
    disjoint_sets_detail::link_sets(parent_.data(), rank_.data(), get(id_, x), get(id_, y));
  }
  template <class Element>
  void union_set(Element x, Element y) {
    Element rx = find_set(x);
    Element ry = find_set(y);
    link(rx, ry);
  }
  template <class Element>
  Element find_set(Element x) {
    return id_to_vertex_[rep_(parent_.data(), get(id_, x))];
  }

  template <std::ranges::input_range IRange>
  std::size_t count_sets(IRange& i_rg) {
    std::size_t count = 0;
    for (auto i : i_rg) {
      if (parent_[i] == i) {
        ++count;
      }
    }
    return count;
  }

  template <std::ranges::input_range IRange>
  void normalize_sets(IRange& i_rg) {
    for (auto i : i_rg) {
      disjoint_sets_detail::normalize_node(parent_.data(), i);
    }
  }

  template <std::ranges::input_range IRange>
  void compress_sets(IRange& i_rg) {
    for (auto i : i_rg) {
      disjoint_sets_detail::find_representative_with_full_compression(parent_.data(), i);
    }
  }

  const ParentContainer& parents() { return parent_; }

 protected:
  template <class Element>
  void extend_sets(Element x, Element y) {
    Index needed = get(id_, x) > get(id_, y) ? get(id_, x) + 1 : get(id_, y) + 1;
    if (needed > parent_.size()) {
      rank_.insert(rank_.end(), needed - rank_.size(), 0);
      for (Index k = parent_.size(); k < needed; ++k) {
        parent_.push_back(k);
      }
    }
  }

  ID id_;
  InverseID id_to_vertex_;
  RankContainer rank_;
  ParentContainer parent_;
  FindCompress rep_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_DISJOINT_SETS_H_
