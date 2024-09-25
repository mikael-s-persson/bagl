// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_IS_HEAP_H_
#define BAGL_BAGL_IS_HEAP_H_

#include <functional>
#include <iterator>

namespace bagl {

namespace is_heap_detail {

template <class RAIter, class StrictWeakOrdering>
bool is_heap(RAIter first, RAIter last, StrictWeakOrdering comp) {
  using Distance = std::iter_difference_t<RAIter>;
  const Distance n = last - first;

  Distance parent = 0;
  for (Distance child = 1; child < n; ++child) {
    if (comp(first[parent], first[child])) {
      return false;
    }
    if ((child & 1) == 0) {
      ++parent;
    }
  }
  return true;
}

}  // namespace is_heap_detail

template <std::random_access_iterator RAIter>
bool is_heap(RAIter first, RAIter last) {
  return is_heap_detail::is_heap(first, last, std::less<>{});
}

template <std::random_access_iterator RAIter, class StrictWeakOrdering>
bool is_heap(RAIter first, RAIter last, StrictWeakOrdering comp) {
  return is_heap_detail::is_heap(first, last, comp);
}

}  // namespace bagl

#endif  // BAGL_BAGL_IS_HEAP_H_
