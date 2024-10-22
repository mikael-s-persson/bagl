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

struct two_bit_color_proxy {
  static constexpr std::uint8_t bits_per_char = std::numeric_limits<std::uint8_t>::digits;
  static constexpr std::uint8_t elements_per_char = bits_per_char / 2;
  static constexpr std::size_t byte_index(std::size_t i) { return i / elements_per_char; }
  static constexpr std::uint8_t bit_offset(std::size_t i) {
    return static_cast<std::uint8_t>((i % elements_per_char) * 2);
  }

  operator two_bit_color_type() const { return static_cast<two_bit_color_type>(((*byte_ptr) >> offset) & 3); }
  const two_bit_color_proxy& operator=(two_bit_color_type c) const {
    (*byte_ptr) = static_cast<std::uint8_t>(((*byte_ptr) & ~(3 << offset)) | (static_cast<std::uint8_t>(c) << offset));
    return *this;
  }

  std::uint8_t* byte_ptr = nullptr;
  std::uint8_t offset = 0;
};

template <typename IndexMap = identity_property_map>
struct two_bit_color_map : property_store_helper<two_bit_color_map<IndexMap>> {
  IndexMap index;
  std::vector<std::uint8_t> data;

  using value_type = two_bit_color_type;

  explicit two_bit_color_map(std::size_t a_n, const IndexMap& a_index = IndexMap())
      : index(a_index), data(two_bit_color_proxy::byte_index(a_n) + 1, 0) {}

  template <typename Key>
  two_bit_color_proxy operator[](Key&& k) {
    auto i = get(index, std::forward<Key>(k));
    std::size_t bi = two_bit_color_proxy::byte_index(i);
    if (bi >= data.size()) [[unlikely]] {
      data.resize(bi + 1, 0);
    }
    return two_bit_color_proxy(&data[bi], two_bit_color_proxy::bit_offset(i));
  }
};

}  // namespace bagl

#endif  // BAGL_BAGL_TWO_BIT_COLOR_MAP_H_
