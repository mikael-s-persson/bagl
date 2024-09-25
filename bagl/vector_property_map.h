// Copyright 2003 Vladimir Prus.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_VECTOR_PROPERTY_MAP_H_
#define BAGL_BAGL_VECTOR_PROPERTY_MAP_H_

#include <iterator>
#include <memory>
#include <vector>

#include "bagl/property_map.h"

namespace bagl {
template <typename T, typename IndexMap = identity_property_map>
class vector_property_map : public put_get_helper<vector_property_map<T, IndexMap>> {
 public:
  using key_type = property_traits_key_t<IndexMap>;
  using value_type = T;
  using reference = decltype(std::declval<std::vector<T>>()[std::declval<property_traits_value_t<IndexMap>>()]);
  using category = lvalue_property_map_tag;

  explicit vector_property_map(IndexMap index = IndexMap())
      : store_(std::make_shared<std::vector<T>>()), index_(std::move(index)) {}

  explicit vector_property_map(std::size_t initial_size, IndexMap index = IndexMap())
      : store_(std::make_shared<std::vector<T>>(initial_size)), index_(std::move(index)) {}

  auto storage_begin() { return store_->begin(); }

  auto storage_end() { return store_->end(); }

  auto storage_begin() const { return store_->cbegin(); }

  auto storage_end() const { return store_->cend(); }

  IndexMap& get_index_map() { return index_; }
  const IndexMap& get_index_map() const { return index_; }

  reference operator[](const key_type& v) const {
    auto i = get(index_, v);
    if (static_cast<std::size_t>(i) >= store_->size()) {
      store_->resize(i + 1, T{});
    }
    return (*store_)[i];
  }

 private:
  // Conceptually, we have a vector of infinite size. For practical
  // purposes, we start with an empty vector and grow it as needed.
  // Note that we cannot store pointer to vector here -- we cannot
  // store pointer to data, because if copy of property map resizes
  // the vector, the pointer to data will be invalidated.
  // I wonder if class 'pmap_ref' is simply needed.
  std::shared_ptr<std::vector<T>> store_;
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
}  // namespace bagl

#endif  // BAGL_BAGL_VECTOR_PROPERTY_MAP_H_
