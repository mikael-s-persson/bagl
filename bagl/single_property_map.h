// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SINGLE_PROPERTY_MAP_H_
#define BAGL_BAGL_SINGLE_PROPERTY_MAP_H_

#include <utility>

#include "bagl/property_map.h"

namespace bagl {

// A single property is one, that regardless of the edge or vertex given,
// will always return a single value.
template <typename Value>
class single_property_map : public put_get_helper<single_property_map<Value>> {
 public:
  using value_type = Value;

  single_property_map() = default;
  explicit single_property_map(value_type value) : value_(std::move(value)) {}

  template <typename T>
  const value_type& operator[](T&& /*unused*/) const {
    return value_;
  }

 private:
  value_type value_;
};

// Same as above, but pretends to be writable as well.
// This class doesn't really make any practical sense, unless you know it won't really be used.
template <typename Value>
class single_writable_property_map {
 public:
  using value_type = Value;

  single_writable_property_map() = default;
  explicit single_writable_property_map(value_type value) : value_(std::move(value)) {}

  template <typename Key>
  friend const Value& get(const single_writable_property_map& me, Key&& /*unused*/) {
    return me.value_;
  }
  template <typename Key>
  friend void put(const single_writable_property_map& /*unused*/, Key&& /*unused*/, const Value& /*unused*/) {}

 private:
  value_type value_;
};

// A property map that always returns a reference to the same object.
template <typename ValueType>
class single_ref_property_map : public put_get_helper<single_ref_property_map<ValueType>> {
 public:
  using value_type = ValueType;

  explicit single_ref_property_map(ValueType& v) : value_(&v) {}

  template <typename T>
  ValueType& operator[](T&& /*unused*/) const {
    return *value_;
  }

 private:
  ValueType* value_;
};

}  // namespace bagl

#endif  // BAGL_BAGL_SINGLE_PROPERTY_MAP_H_
