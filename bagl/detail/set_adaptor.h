// Copyright Jeremy Siek 2001.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_SET_ADAPTOR_H_
#define BAGL_BAGL_DETAIL_SET_ADAPTOR_H_

#include <algorithm>
#include <compare>
#include <set>
#include <unordered_set>

namespace bagl {

template <typename K, typename C, typename A, typename T>
[[nodiscard]] bool set_contains(const std::set<K, C, A>& s, const T& x) {
  return s.find(x) != s.end();
}

template <typename K, typename H, typename C, typename A, typename T>
[[nodiscard]] bool set_contains(const std::unordered_set<K, H, C, A>& s, const T& x) {
  return s.find(x) != s.end();
}

template <typename K, typename C, typename A>
[[nodiscard]] bool set_equal(const std::set<K, C, A>& x, const std::set<K, C, A>& y) {
  return x == y;
}

// Not the same as lexicographical_compare_3way applied to std::set.
// this is equivalent semantically to bitset::operator<()
template <typename K, typename C, typename A>
[[nodiscard]] std::strong_ordering set_lex_order(const std::set<K, C, A>& x, const std::set<K, C, A>& y) {
  auto xi = x.begin();
  auto xend = x.end();
  auto yi = y.begin();
  auto yend = y.end();
  for (; xi != xend && yi != yend; ++xi, ++yi) {
    const auto comp_res = *xi <=> *yi;
    if (comp_res == std::strong_ordering::less) {
      return std::strong_ordering::less;
    }
    if (comp_res == std::strong_ordering::greater) {
      return std::strong_ordering::greater;
    }
  }
  if (xi == xend) {
    return (yi == yend) ? std::strong_ordering::equal : std::strong_ordering::greater;
  }
  return std::strong_ordering::less;
}

template <typename K, typename C, typename A>
void set_clear(std::set<K, C, A>& x) {
  x.clear();
}

template <typename K, typename C, typename A>
[[nodiscard]] bool set_empty(const std::set<K, C, A>& x) {
  return x.empty();
}

template <typename K, typename C, typename A, typename T>
void set_insert(std::set<K, C, A>& x, const T& a) {
  x.insert(a);
}

template <typename K, typename C, typename A, typename T>
void set_remove(std::set<K, C, A>& x, const T& a) {
  x.erase(a);
}

template <typename K, typename C, typename A>
void set_intersect(const std::set<K, C, A>& x, const std::set<K, C, A>& y, std::set<K, C, A>& z) {
  z.clear();
  std::set_intersection(x.begin(), x.end(), y.begin(), y.end(), std::inserter(z));
}

template <typename K, typename C, typename A>
void set_union(const std::set<K, C, A>& x, const std::set<K, C, A>& y, std::set<K, C, A>& z) {
  z.clear();
  std::set_union(x.begin(), x.end(), y.begin(), y.end(), std::inserter(z));
}

template <typename K, typename C, typename A>
void set_difference(const std::set<K, C, A>& x, const std::set<K, C, A>& y, std::set<K, C, A>& z) {
  z.clear();
  std::set_difference(x.begin(), x.end(), y.begin(), y.end(), std::inserter(z, z.begin()));
}

template <typename K, typename C, typename A>
[[nodiscard]] bool set_subset(const std::set<K, C, A>& x, const std::set<K, C, A>& y) {
  return std::includes(x.begin(), x.end(), y.begin(), y.end());
}

// Shit, can't implement this without knowing the size of the universe.
template <typename K, typename C, typename A>
void set_compliment(const std::set<K, C, A>& /*x*/, std::set<K, C, A>& z) {
  z.clear();
}

}  // namespace bagl

#endif  // BAGL_BAGL_DETAIL_SET_ADAPTOR_H_
