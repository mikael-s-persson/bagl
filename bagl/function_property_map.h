// Copyright 2012, Philipp Moeller
// Author: Philipp Moeller
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_FUNCTION_PROPERTY_MAP_H_
#define BAGL_BAGL_FUNCTION_PROPERTY_MAP_H_

#include <type_traits>
#include <utility>

#include "bagl/property_map.h"

namespace bagl {

template <typename Func, typename Key, typename Ret = decltype(std::declval<Func>()(std::declval<Key>()))>
class function_property_map : public put_get_helper<function_property_map<Func, Key, Ret> > {
 public:
  using key_type = Key;
  using reference = Ret;
  using value_type = std::decay_t<Ret>;

  using category = std::conditional_t<std::is_reference_v<Ret> && !std::is_const_v<Ret>, lvalue_property_map_tag,
                                      readable_property_map_tag>;

  explicit function_property_map(Func f = Func()) : f_(f) {}

  reference operator[](const Key& k) const { return f_(k); }

 private:
  Func f_;
};

template <typename Key, typename Func>
auto make_function_property_map(const Func& f) {
  return function_property_map<Func, Key>(f);
}

template <typename Key, typename Ret, typename Func>
auto make_function_property_map(const Func& f) {
  return function_property_map<Func, Key, Ret>(f);
}

}  // namespace bagl

#endif  // BAGL_BAGL_FUNCTION_PROPERTY_MAP_H_
