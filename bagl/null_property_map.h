// Copyright 2007 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_NULL_PROPERTY_MAP_H_
#define BAGL_BAGL_NULL_PROPERTY_MAP_H_

#include "bagl/property_map.h"

namespace bagl {
// A null property is somewhat like the inverse of the constant
// property map except that instead of returning a single value,
// this eats any writes and cannot be read from.

template <typename Key, typename Value>
struct null_property_map {
  using key_type = Key;
  using value_type = Value;
  using reference = void;
  using category = writable_property_map_tag;
};

// The null_property_map<K,V> only has a put() function.
template <typename K, typename V>
void put(const null_property_map<K, V>& /*unused*/, const K& /*unused*/, const V& /*unused*/) {}

// A helper function for intantiating null property maps.
template <typename Key, typename Value>
auto make_null_property() {
  return null_property_map<Key, Value>();
}

}  // namespace bagl

#endif  // BAGL_BAGL_NULL_PROPERTY_MAP_H_
