// Copyright (C) 2005-2010 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock, Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ONE_BIT_COLOR_MAP_H_
#define BAGL_BAGL_ONE_BIT_COLOR_MAP_H_

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>

#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

// One bit per color property map (gray and black are the same, green is not
// supported)
enum class one_bit_color_type : std::uint8_t { white = 0, not_white = 1 };

template <>
struct color_traits<one_bit_color_type> {
  static one_bit_color_type white() { return one_bit_color_type::white; }
  static one_bit_color_type gray() { return one_bit_color_type::not_white; }
  static one_bit_color_type black() { return one_bit_color_type::not_white; }
};

template <typename IndexMap = identity_property_map>
struct one_bit_color_map {
  static constexpr int bits_per_char = std::numeric_limits<std::uint8_t>::digits;
  std::size_t n;
  IndexMap index;
  std::shared_ptr<std::uint8_t[]> data;

  using value_type = one_bit_color_type;
  using reference = void;

  explicit one_bit_color_map(std::size_t a_n, const IndexMap& a_index = IndexMap())
      : n(a_n), index(a_index), data(new std::uint8_t[(n + bits_per_char - 1) / bits_per_char]()) {}
};

template <typename Key, concepts::ReadablePropertyMap<Key> IndexMap>
auto get(const one_bit_color_map<IndexMap>& pm, const Key& key) {
  constexpr int bits_per_char = one_bit_color_map<IndexMap>::bits_per_char;
  auto i = get(pm.index, key);
  assert(i < pm.n);
  return static_cast<one_bit_color_type>((pm.data.get()[i / bits_per_char] >> (i % bits_per_char)) & 1);
}

template <typename Key, concepts::ReadablePropertyMap<Key> IndexMap>
void put(const one_bit_color_map<IndexMap>& pm, const Key& key, one_bit_color_type value) {
  constexpr int bits_per_char = one_bit_color_map<IndexMap>::bits_per_char;
  auto i = get(pm.index, key);
  assert(i < pm.n);
  const std::size_t byte_num = i / bits_per_char;
  const std::size_t bit_position = (i % bits_per_char);
  pm.data.get()[byte_num] = static_cast<std::uint8_t>((pm.data.get()[byte_num] & ~(1 << bit_position)) |
                                                      (static_cast<std::uint8_t>(value) << bit_position));
}

template <typename IndexMap>
auto make_one_bit_color_map(std::size_t n, const IndexMap& index_map) {
  return one_bit_color_map<IndexMap>(n, index_map);
}

}  // namespace bagl

#endif  // BAGL_BAGL_ONE_BIT_COLOR_MAP_H_
