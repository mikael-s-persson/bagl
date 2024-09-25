// Copyright (C) 2005-2006 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock, Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TWO_BIT_COLOR_MAP_H_
#define BAGL_BAGL_TWO_BIT_COLOR_MAP_H_

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>
#include <utility>

#include "bagl/properties.h"
#include "bagl/property_map.h"

// Two bit per color property map

namespace bagl {

enum class two_bit_color_type : std::uint8_t { white = 0, gray = 1, green = 2, black = 3 };

template <>
struct color_traits<two_bit_color_type> {
  static constexpr two_bit_color_type white() { return two_bit_color_type::white; }
  static constexpr two_bit_color_type gray() { return two_bit_color_type::gray; }
  static constexpr two_bit_color_type green() { return two_bit_color_type::green; }
  static constexpr two_bit_color_type black() { return two_bit_color_type::black; }
};

template <typename IndexMap = identity_property_map>
struct two_bit_color_map {
  std::size_t n;
  IndexMap index;
  std::shared_ptr<std::uint8_t[]> data;

  static constexpr int bits_per_char = std::numeric_limits<std::uint8_t>::digits;
  static constexpr int elements_per_char = bits_per_char / 2;
  using key_type = property_traits_key_t<IndexMap>;
  using value_type = two_bit_color_type;
  using reference = void;
  using category = read_write_property_map_tag;

  explicit two_bit_color_map(std::size_t a_n, const IndexMap& a_index = IndexMap())
      : n(a_n), index(a_index), data(new std::uint8_t[(n + elements_per_char - 1) / elements_per_char]()) {}
};

template <typename IndexMap>
two_bit_color_type get(const two_bit_color_map<IndexMap>& pm, property_traits_key_t<IndexMap> key) {
  static constexpr int elements_per_char = two_bit_color_map<IndexMap>::elements_per_char;
  auto i = get(pm.index, key);
  assert(i < pm.n);
  const std::size_t byte_num = i / elements_per_char;
  const std::size_t bit_position = ((i % elements_per_char) * 2);
  return two_bit_color_type((pm.data.get()[byte_num] >> bit_position) & 3);
}

template <typename IndexMap>
void put(const two_bit_color_map<IndexMap>& pm, property_traits_key_t<IndexMap> key, two_bit_color_type value) {
  static constexpr int elements_per_char = two_bit_color_map<IndexMap>::elements_per_char;
  auto i = get(pm.index, key);
  assert(i < pm.n);
  const std::size_t byte_num = i / elements_per_char;
  const std::size_t bit_position = ((i % elements_per_char) * 2);
  pm.data.get()[byte_num] = static_cast<std::uint8_t>((pm.data.get()[byte_num] & ~(3 << bit_position)) |
                                                       (std::to_underlying(value) << bit_position));
}

template <typename IndexMap>
two_bit_color_map<IndexMap> make_two_bit_color_map(std::size_t n, const IndexMap& index_map) {
  return two_bit_color_map<IndexMap>(n, index_map);
}

}  // end namespace bagl

#endif  // BAGL_BAGL_TWO_BIT_COLOR_MAP_H_
