// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_INDIRECT_CMP_H_
#define BAGL_BAGL_INDIRECT_CMP_H_

#include <functional>

#include "bagl/property_map.h"

namespace bagl {

template <class PMap, class Compare>
class indirect_cmp {
 public:
  explicit indirect_cmp(const PMap& df, const Compare& c = Compare()) : d_(df), cmp_(c) {}

  template <class A, class B>
    requires concepts::ReadablePropertyMap<PMap, A> && concepts::ReadablePropertyMap<PMap, B>
  bool operator()(const A& u, const B& v) const {
    const auto& du = get(d_, u);
    const auto& dv = get(d_, v);
    return cmp_(du, dv);
  }

 protected:
  PMap d_;
  Compare cmp_;
};

template <typename Compare, typename PMap>
auto make_indirect_cmp(const Compare& cmp, PMap pmap) {
  return indirect_cmp<PMap, Compare>{pmap, cmp};
}

template <class PMap>
class indirect_pmap {
 public:
  explicit indirect_pmap(const PMap& df) : d_(df) {}

  template <typename K>
    requires concepts::ReadablePropertyMap<PMap, K>
  auto operator()(const K& u) const {
    return get(d_, u);
  }

 protected:
  PMap d_;
};

template <typename PMap>
auto make_indirect_pmap(PMap pmap) {
  return indirect_pmap<PMap>{pmap};
}

}  // namespace bagl

#endif  // BAGL_BAGL_INDIRECT_CMP_H_
