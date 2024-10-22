// Copyright 2012, Philipp Moeller
// Author: Philipp Moeller
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_TRANSFORM_VALUE_PROPERTY_MAP_H_
#define BAGL_BAGL_TRANSFORM_VALUE_PROPERTY_MAP_H_

#include <type_traits>
#include <utility>

#include "bagl/property_map.h"

namespace bagl {

template <typename Func, typename PM>
class transform_value_property_map : public put_get_helper<transform_value_property_map<Func, PM>> {
 public:
  using self = transform_value_property_map<Func, PM>;

  transform_value_property_map(Func f, PM pm) : f_(f), pm_(pm) {}

  template <typename Key>
  decltype(auto) operator[](const Key& k) const {
    return f_(get(pm_, k));
  }

  using reference = decltype(std::declval<Func>()(std::declval<property_traits_value_t<PM>>()));
  using value_type = std::remove_cv_t<std::remove_reference_t<reference>>;

 private:
  Func f_;
  PM pm_;
};

template <typename Key, typename Func>
auto transformed_property_map(Func&& f) {
  return transform_value_property_map{std::forward<Func>(f), typed_identity_property_map<Key>{}};
}

}  // namespace bagl

#endif  // BAGL_BAGL_TRANSFORM_VALUE_PROPERTY_MAP_H_
