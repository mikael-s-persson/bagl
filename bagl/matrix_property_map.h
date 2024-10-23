// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MATRIX_PROPERTY_MAP_H_
#define BAGL_BAGL_MATRIX_PROPERTY_MAP_H_

#include <vector>

#include "bagl/property_map.h"

namespace bagl {

// This property map is built specifically for property maps over
// matrices. Like the basic property map over a container, this builds
// the property abstraction over a matrix (vector of vectors)
// and returns property maps over the nested containers.
template <typename Matrix, typename IndexMap = identity_property_map>
struct matrix_property_map : put_get_helper<matrix_property_map<Matrix, IndexMap>> {
  using value_type =
      iterator_property_map<decltype(std::declval<Matrix>()[std::declval<property_traits_value_t<IndexMap>>()].begin()),
                            IndexMap>;

  matrix_property_map() = default;
  explicit matrix_property_map(Matrix* matrix, IndexMap index = IndexMap())
      : matrix_(matrix), index_(std::move(index)) {}

  template <typename Key>
  auto operator[](Key&& k) const {
    auto i = get(index_, std::forward<Key>(k));
    if constexpr (!std::is_const_v<Matrix>) {
      if (i >= matrix_->size()) [[unlikely]] {
        matrix_->resize(i);
        if (matrix_[0].empty()) {
          for (auto& v : (*matrix_)) {
            v.resize(i + 1);
          }
        } else {
          auto default_value = (*matrix_)[0].back();
          for (auto& v : (*matrix_)) {
            v.resize(i + 1, default_value);
          }
        }
      }
    }
    return iterator_property_map((*matrix_)[i].begin(), index_);
  }

 private:
  Matrix* matrix_ = nullptr;
  IndexMap index_;
};

template <typename Value, typename IndexMap = identity_property_map>
struct matrix_property_store : property_store_helper<matrix_property_store<Value, IndexMap>> {
  explicit matrix_property_store(std::size_t initial_size, IndexMap index = IndexMap(), Value default_value = Value{})
      : matrix_((initial_size == 0 ? std::size_t{1} : initial_size),
                std::vector<Value>(initial_size + 1, default_value)),
        index_(std::move(index)) {}

  explicit matrix_property_store(IndexMap index = IndexMap(), Value default_value = Value{})
      : matrix_property_store(0, std::move(index), std::move(default_value)) {}

  auto& storage() { return matrix_; }
  const auto& storage() const { return matrix_; }

  IndexMap& get_index_map() { return index_; }
  const IndexMap& get_index_map() const { return index_; }

  auto ref() { return matrix_property_map(&matrix_, index_); }
  auto ref() const { return matrix_property_map(&matrix_, index_); }

  template <typename Key>
  decltype(auto) operator[](Key&& k) {
    return ref()[std::forward<Key>(k)];
  }

 private:
  std::vector<std::vector<Value>> matrix_;
  IndexMap index_;
};

}  // namespace bagl

#endif
