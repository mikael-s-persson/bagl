// Copyright 2003 Vladimir Prus.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_VECTOR_PROPERTY_MAP_H_
#define BAGL_BAGL_VECTOR_PROPERTY_MAP_H_

#include <iterator>
#include <memory>
#include <ranges>
#include <utility>
#include <vector>

#include "bagl/property_map.h"

namespace bagl {
template <typename T, typename IndexMap = identity_property_map>
class vector_property_map : public property_store_helper<vector_property_map<T, IndexMap>> {
 public:
  using value_type = T;

  explicit vector_property_map(std::size_t initial_size, IndexMap index = IndexMap(), T default_value = T{})
      : store_(initial_size + 1, default_value), index_(std::move(index)) {}

  explicit vector_property_map(IndexMap index = IndexMap(), T default_value = T{})
      : vector_property_map(0, std::move(index), std::move(default_value)) {}

  auto storage_range() { return std::ranges::ref_view(store_); }
  auto storage_range() const { return std::ranges::ref_view(store_); }

  IndexMap& get_index_map() { return index_; }
  const IndexMap& get_index_map() const { return index_; }

  template <typename Key>
  decltype(auto) operator[](Key&& k) {
    auto i = static_cast<std::size_t>(get(index_, std::forward<Key>(k)));
    if (i >= store_.size()) [[unlikely]] {
      T default_value = store_.back();
      store_.resize(i + 1, default_value);
    }
    return store_[i];
  }

 private:
  std::vector<T> store_;
  IndexMap index_;
};

template <typename T, typename IndexMap>
auto make_vector_property_map(IndexMap index) {
  return vector_property_map<T, IndexMap>(std::move(index));
}

template <typename T, typename IndexMap>
auto make_vector_property_map(std::size_t initial_size, IndexMap index) {
  return vector_property_map<T, IndexMap>(initial_size, std::move(index));
}

template <typename T, typename IndexMap>
auto make_vector_property_map(std::size_t initial_size, IndexMap index, const T& default_value) {
  return vector_property_map<T, IndexMap>(initial_size, std::move(index), default_value);
}

}  // namespace bagl

#endif  // BAGL_BAGL_VECTOR_PROPERTY_MAP_H_
