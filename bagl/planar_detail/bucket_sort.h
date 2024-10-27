// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLANAR_DETAIL_BUCKET_SORT_H_
#define BAGL_BAGL_PLANAR_DETAIL_BUCKET_SORT_H_

#include <algorithm>
#include <vector>

#include "bagl/property_map.h"

namespace bagl {

template <typename TupleType, int N, typename PropertyMapWrapper = identity_property_map>
struct property_map_tuple_adaptor
    : public put_get_helper<property_map_tuple_adaptor<TupleType, N, PropertyMapWrapper> > {
  using value_type = typename PropertyMapWrapper::value_type;

  property_map_tuple_adaptor() = default;

  explicit property_map_tuple_adaptor(PropertyMapWrapper wrapper_map) : wrapper_map_(wrapper_map) {}

  template <typename Key>
  value_type operator[](Key&& x) const {
    return get(wrapper_map_, std::get<N>(std::forward<Key>(x)));
  }

  static constexpr int n = N;
  PropertyMapWrapper wrapper_map_;
};

template <typename ForwardIterator, typename ItemToRankMap>
void bucket_sort(ForwardIterator begin, ForwardIterator end, ItemToRankMap rank) {
  std::stable_sort(begin, end, [rank](auto x, auto y) { return get(rank, x) < get(rank, y); });
}

template <typename ForwardIterator>
void bucket_sort(ForwardIterator begin, ForwardIterator end) {
  std::stable_sort(begin, end);
}

}  // namespace bagl

#endif  // BAGL_BAGL_PLANAR_DETAIL_BUCKET_SORT_H_
