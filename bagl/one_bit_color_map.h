// Copyright (C) 2005-2010 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock, Douglas Gregor, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_ONE_BIT_COLOR_MAP_H_
#define BAGL_BAGL_ONE_BIT_COLOR_MAP_H_

#include <algorithm>
#include <cassert>
#include <cstdint>
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

struct one_bit_color_proxy {
  static constexpr std::uint8_t bits_per_char = std::numeric_limits<std::uint8_t>::digits;
  static constexpr std::size_t byte_index(std::size_t i) { return i / bits_per_char; }
  static constexpr std::uint8_t bit_offset(std::size_t i) { return static_cast<std::uint8_t>(i % bits_per_char); }

  operator one_bit_color_type() const { return static_cast<one_bit_color_type>(((*byte_ptr) >> offset) & 1); }
  const one_bit_color_proxy& operator=(one_bit_color_type c) const {
    (*byte_ptr) = static_cast<std::uint8_t>(((*byte_ptr) & ~(1 << offset)) | (static_cast<std::uint8_t>(c) << offset));
    return *this;
  }

  std::uint8_t* byte_ptr = nullptr;
  std::uint8_t offset = 0;
};

template <typename IndexMap = identity_property_map>
struct one_bit_color_map : property_store_helper<one_bit_color_map<IndexMap>> {
  IndexMap index;
  std::vector<std::uint8_t> data;

  using value_type = one_bit_color_type;

  explicit one_bit_color_map(std::size_t a_n, IndexMap a_index = IndexMap())
      : index(std::move(a_index)), data(one_bit_color_proxy::byte_index(a_n) + 1, 0) {}

  template <typename Key>
  one_bit_color_proxy operator[](Key&& k) {
    auto i = get(index, std::forward<Key>(k));
    std::size_t bi = one_bit_color_proxy::byte_index(i);
    if (bi >= data.size()) [[unlikely]] {
      data.resize(bi + 1, 0);
    }
    return one_bit_color_proxy(&data[bi], one_bit_color_proxy::bit_offset(i));
  }
};

}  // namespace bagl

#endif  // BAGL_BAGL_ONE_BIT_COLOR_MAP_H_
