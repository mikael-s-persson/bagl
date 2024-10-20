// Copyright 2024 Mikael Persson

#ifndef BAGL_BAGL_ZIP_RANGE_H_
#define BAGL_BAGL_ZIP_RANGE_H_

#include <cassert>
#include <tuple>
#include <limits>
#include <iterator>
#include <ranges>

#include "bagl/graph_traits.h"

namespace bagl {

namespace zip_range_detail {

template <typename Category1, typename Category2>
using least_permissive_pair =
    std::conditional_t<std::is_convertible<Category1, Category2>::value,
                       Category2, Category1>;

template <typename... CategoryTail>
struct least_permissive_impl {
  static_assert(sizeof...(CategoryTail) != 0);
};

template <typename Category1>
struct least_permissive_impl<Category1> {
  using type = Category1;
};

template <typename Category1, typename Category2>
struct least_permissive_impl<Category1, Category2> {
  using type = least_permissive_pair<Category1, Category2>;
};

template <typename Category1, typename Category2, typename... CategoryTail>
struct least_permissive_impl<Category1, Category2, CategoryTail...> {
  using type =
      typename least_permissive_impl<least_permissive_pair<Category1, Category2>,
                                   CategoryTail...>::type;
};

template <typename... CategoryTail>
using least_permissive = typename least_permissive_impl<CategoryTail...>::type;

template <typename Category>
using reduce_to_std_iter_category = std::conditional_t<
    std::is_convertible_v<Category, std::random_access_iterator_tag>,
    std::random_access_iterator_tag,
    std::conditional_t<
        std::is_convertible_v<Category, std::bidirectional_iterator_tag>,
        std::bidirectional_iterator_tag,
        std::conditional_t<
            std::is_convertible_v<Category, std::forward_iterator_tag>,
            std::forward_iterator_tag,
            std::conditional_t<
                std::is_convertible_v<Category, std::input_iterator_tag>,
                std::input_iterator_tag,
                std::conditional_t<
                    std::is_convertible_v<Category, std::output_iterator_tag>,
                    std::output_iterator_tag, void>>>>>;

template <typename... Iters>
using compute_iter_category = least_permissive<reduce_to_std_iter_category<
    typename std::iterator_traits<Iters>::iterator_category>...>;

template <typename T>
T&& variadic_min(T&& val) {
  return std::forward<T>(val);
}

template <typename T0, typename T1, typename... Ts>
auto variadic_min(T0&& val1, T1&& val2, Ts&&... vs) {
  return variadic_min((val1 < val2) ? val1 : val2, std::forward<Ts>(vs)...);
}

}  // namespace zip_range_detail

// A zip_iterator is an iterator that combines multiple underlying iterators
// into a single iterator that produces a tuple of the underlying values when
// dereferenced.
//
// See make_zip_iterator and zip_range functions for convenient ways to create
// zip iterators with template argument deduction.
//
// Example:
// We could copy a vector into another like this:
// void copy(const std::vector<int>& src, std::vector<int>* dest) {
//   for (auto src_dest : zip_range(src, *dest)) {
//     get<1>(src_dest) = get<0>(src_dest);
//   }
// }
//
template <typename... Iters>
class zip_iterator {
 public:
  using iterator_category = zip_range_detail::compute_iter_category<Iters...>;
  using difference_type = int;
  using value_type = std::tuple<decltype(*std::declval<Iters>())...>;
  using pointer = value_type*;
  using reference = value_type;

  using self = zip_iterator<Iters...>;

  // Universal constructor:
  template <typename... OtherIters>
  explicit zip_iterator(OtherIters... it) : it_tuple_(std::move(it)...) {}

  // Default constructor:
  zip_iterator() : it_tuple_() {}

  // Dereference operators:
  reference operator*() const {
    return dereference();
  }
  pointer operator->() const {
    assert(false && "operator->() is only supported at run-time.");
    return &*(*this);
  }

  // Increment / Decrement operators:
  self& operator++() {
    increment();
    return *this;
  }
  self operator++(int) {
    self old_it = *this;
    increment();
    return old_it;  // NRVO
  }
  self& operator--() {
    decrement();
    return *this;
  }
  self operator--(int) {
    self old_it = *this;
    decrement();
    return old_it;  // NRVO
  }

  // Random-access operators:
  self& operator+=(difference_type i) {
    advance(i);
    return *this;
  }
  self& operator-=(difference_type i) {
    advance(-i);
    return *this;
  }

  self operator-(int i) const {
    self result = *this;
    result -= i;
    return result;  // NRVO
  }

  self operator+(int i) {
    self result = *this;
    result += i;
    return result;  // NRVO
  }

  friend self operator+(int i, const self& it) {
    return it + i;
  }

  reference operator[](difference_type i) const {
    self adv_it = *this;
    adv_it.advance(i);
    return adv_it.dereference();
  }

  bool operator==(const self& rhs) const {
    return is_equal(rhs);
  }
  bool operator!=(const self& rhs) const {
    return !is_equal(rhs);
  }

  bool operator<(const self& rhs) const {
    return distance_to(rhs) > 0;
  }
  bool operator>(const self& rhs) const {
    return rhs < *this;
  }
  bool operator<=(const self& rhs) const {
    return !(rhs < *this);
  }
  bool operator>=(const self& rhs) const {
    return !(rhs > *this);
  }

 private:

  // Implementation of the IteratorFacade requirements:
  template <size_t... Ids>
  reference dereference(std::index_sequence<Ids...>  /*ids*/) const {
    return {(*std::get<Ids>(it_tuple_))...};
  }
  template <size_t... Ids>
  void increment(std::index_sequence<Ids...>  /*ids*/) {
    ((void)++std::get<Ids>(it_tuple_), ...);
  }
  template <size_t... Ids>
  void decrement(std::index_sequence<Ids...>  /*ids*/) {
    ((void)--std::get<Ids>(it_tuple_), ...);
  }
  template <size_t... Ids>
  bool is_equal(const self& rhs, std::index_sequence<Ids...>  /*ids*/) const {
    // Only require one pair of iterators to match such that iterations are
    // stopped by the shortest range if all ranges don't match.
    return ((std::get<Ids>(it_tuple_) == std::get<Ids>(rhs.it_tuple_)) || ...);
  }
  template <size_t... Ids>
  int distance_to(const self& rhs,
                 std::index_sequence<Ids...>  /*ids*/) const {
    // Return the smallest distance between iterators because that is where
    // iterations will stop.
    return zip_range_detail::variadic_min(
        (std::get<Ids>(rhs.it_tuple_) - std::get<Ids>(it_tuple_))...);
  }
  template <size_t... Ids>
  void advance(int n, std::index_sequence<Ids...>  /*ids*/) {
    ((void)(std::get<Ids>(it_tuple_) += n), ...);
  }

  using IterIndexSeq = std::make_index_sequence<sizeof...(Iters)>;

  value_type dereference() const { return dereference(IterIndexSeq()); }
  void increment() { increment(IterIndexSeq()); }
  void decrement() { decrement(IterIndexSeq()); }
  bool is_equal(const zip_iterator& rhs) const {
    return is_equal(rhs, IterIndexSeq());
  }
  int distance_to(const zip_iterator& rhs) const {
    return distance_to(rhs, IterIndexSeq());
  }
  void advance(int n) { advance(n, IterIndexSeq()); }

  std::tuple<Iters...> it_tuple_;
};

// Deduction guide
template <typename... OtherIters>
zip_iterator(OtherIters... it) -> zip_iterator<std::decay_t<OtherIters>...>;

// Factory function that conveniently creates a zip iterator object
// using template argument deduction to infer the types of the underlying
// iterators.
// iters: The underlying iterators.
// For example, when zipping two vectors v1 and v2:
//   auto zip_it = make_zip_iterator(v1.begin(), v2.begin());
template <typename... UnderlyingIters>
auto make_zip_iterator(UnderlyingIters&&... iters) {
  return zip_iterator<std::decay_t<UnderlyingIters>...>(
      std::forward<UnderlyingIters>(iters)...);
}

// A zipped_range is a range that combines multiple underlying ranges
// into a single iterator that produces a tuple of the underlying values when
// dereferenced.
template <std::ranges::range... Ranges>
class zipped_range
    : public std::ranges::view_interface<zipped_range<Ranges...>> {
 public:
  using RangeTuple = std::tuple<Ranges...>;
  using ZipIter = zip_iterator<std::ranges::iterator_t<Ranges>...>;

  template <typename... OtherRanges>
  explicit zipped_range(OtherRanges&&... ranges)
      : range_tuple(std::forward<OtherRanges>(ranges)...) {}

  auto begin() const {
    return std::apply(
        [](const Ranges&... base_ranges) {
          using std::begin;
          return ZipIter(begin(base_ranges)...);
        },
        range_tuple);
  }
  auto begin() {
    return std::apply(
        [](Ranges&... base_ranges) {
          using std::begin;
          return ZipIter(begin(base_ranges)...);
        },
        range_tuple);
  }
  auto end() const {
    return std::apply(
        [](const Ranges&... base_ranges) {
          using std::end;
          return ZipIter(end(base_ranges)...);
        },
        range_tuple);
  }
  auto end() {
    return std::apply(
        [](Ranges&... base_ranges) {
          using std::end;
          return ZipIter(end(base_ranges)...);
        },
        range_tuple);
  }

  RangeTuple range_tuple;
};

// Factory function that conveniently creates a zip iterator range
// using template argument deduction to infer the type of the underlying
// ranges.
// ranges: The underlying ranges.
// For example, when zipping two vectors v1 and v2:
//   auto zip_range = zip_range(v1, v2);
template <std::ranges::range... Ranges>
auto zip_range(Ranges&&... ranges) {
  return zipped_range<std::decay_t<Ranges>...>(std::forward<Ranges>(ranges)...);
}
template <std::ranges::range... Ranges>
auto zip_range_ref(Ranges&&... ranges) {
  return zip_range(std::ranges::ref_view(std::forward<Ranges>(ranges))...);
}

// Factory function that conveniently creates an enumerated iterator range
// where each element of the range is paired with its integer index / counter
// represented as an `int` (size of range should be less than max of int).
// For example, when enumerate a vector v:
//   auto enumerated_range = enumerate_range(v);
// or:
//   for (auto [i, x] : enumerate_range(v)) { .. }
template <std::ranges::range Range>
auto enumerate_range(Range&& range) {
  return zip_range(std::ranges::iota_view(0, std::numeric_limits<int>::max()),
                   std::forward<Range>(range));
}

}  // namespace bagl

#endif  // BAGL_BAGL_ZIP_RANGE_H_
