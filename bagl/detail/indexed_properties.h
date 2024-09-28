// Copyright 2005 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

// Indexed properties -- used for CSR and CSR-like graphs

#ifndef BAGL_BAGL_DETAIL_INDEXED_PROPERTIES_H_
#define BAGL_BAGL_DETAIL_INDEXED_PROPERTIES_H_

#include <algorithm>
#include <climits>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"

namespace bagl::indexed_detail {

template <typename Derived, typename Property, typename Descriptor, typename IndexMap>
class indexed_vertex_properties {
 public:
  using vertex_property_type = no_property;
  using vertex_bundled = Property;

  // Directly access a vertex or edge bundle
  Property& operator[](Descriptor v) { return vertex_properties_[get(vertex_index, derived(), v)]; }

  const Property& operator[](Descriptor v) const { return vertex_properties_[get(vertex_index, derived(), v)]; }

  auto get_vertex_bundle(const IndexMap& index_map = IndexMap()) {
    return make_iterator_property_map(vertex_properties_.begin(), index_map);
  }

  auto get_vertex_bundle(const IndexMap& index_map = IndexMap()) const {
    return make_iterator_property_map(vertex_properties_.begin(), index_map);
  }

 protected:
  // Default-construct with no property values
  indexed_vertex_properties() = default;

  // Initialize with n default-constructed property values
  explicit indexed_vertex_properties(std::size_t n) : vertex_properties_(n) {}

 public:
  // Clear the properties vector
  void clear() { vertex_properties_.clear(); }

  // Resize the properties vector
  void resize(std::size_t n) { vertex_properties_.resize(n); }

  // Reserve space in the vector of properties
  void reserve(std::size_t n) { vertex_properties_.reserve(n); }

  // Add a new property value to the back
  void push_back(const Property& prop) { vertex_properties_.push_back(prop); }

  // Write an element by raw index
  void write_by_index(std::size_t idx, const Property& prop) { vertex_properties_[idx] = prop; }

  // Access to the derived object
  Derived& derived() { return *static_cast<Derived*>(this); }
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

 private:  // should be private, but friend templates not portable
  std::vector<Property> vertex_properties_;
};

template <typename Derived, typename Descriptor, typename IndexMap>
class indexed_vertex_properties<Derived, void, Descriptor, IndexMap> {
  struct secret {};

 public:
  using vertex_property_type = no_property;
  using vertex_bundled = void;

  secret operator[](secret /*unused*/) { return secret(); }

  auto get_vertex_bundle() const { return secret(); }

 protected:
  // All operations do nothing.
  indexed_vertex_properties() = default;
  explicit indexed_vertex_properties(std::size_t /*unused*/) {}

 public:
  void clear() {}
  void resize(std::size_t /*unused*/) {}
  void reserve(std::size_t /*unused*/) {}
};

template <typename Derived, typename Property, typename Descriptor, typename IndexMap>
class indexed_edge_properties {
 public:
  using edge_property_type = no_property;
  using edge_bundled = Property;
  using edge_push_back_type = Property;

  // Directly access a edge or edge bundle
  Property& operator[](Descriptor v) { return edge_properties_[get(edge_index, derived(), v)]; }

  const Property& operator[](Descriptor v) const { return edge_properties_[get(edge_index, derived(), v)]; }

  auto get_edge_bundle(const IndexMap& index_map = IndexMap()) {
    return make_iterator_property_map(edge_properties_.begin(), index_map);
  }

  auto get_edge_bundle(const IndexMap& index_map = IndexMap()) const {
    return make_iterator_property_map(edge_properties_.begin(), index_map);
  }

 protected:
  // Default-construct with no property values
  indexed_edge_properties() = default;

  // Initialize with n default-constructed property values
  explicit indexed_edge_properties(std::size_t n) : edge_properties_(n) {}

  // Get the size of the properties vector
  [[nodiscard]] std::size_t size() const { return edge_properties_.size(); }

  // Clear the properties vector
  void clear() { edge_properties_.clear(); }

  // Resize the properties vector
  void resize(std::size_t n) { edge_properties_.resize(n); }

  // Reserve space in the vector of properties
  void reserve(std::size_t n) { edge_properties_.reserve(n); }

  // Write an element by raw index
  void write_by_index(std::size_t idx, const Property& prop) { edge_properties_[idx] = prop; }

 public:
  // Add a new property value to the back
  void push_back(const Property& prop) { edge_properties_.push_back(prop); }

  // Move range of properties backwards
  void move_range(std::size_t src_begin, std::size_t src_end, std::size_t dest_begin) {
    std::copy_backward(edge_properties_.begin() + src_begin, edge_properties_.begin() + src_end,
                       edge_properties_.begin() + dest_begin + (src_end - src_begin));
  }

  auto as_range() const { return std::ranges::ref_view(edge_properties_); }

 private:
  // Access to the derived object
  Derived& derived() { return *static_cast<Derived*>(this); }
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

 public:  // should be private, but friend templates not portable
  std::vector<Property> edge_properties_;
};

template <typename Derived, typename Descriptor, typename IndexMap>
class indexed_edge_properties<Derived, void, Descriptor, IndexMap> {
  struct secret {};

 public:
  using edge_property_type = no_property;
  using edge_bundled = void;
  using edge_push_back_type = void*;

  secret operator[](secret) { return secret(); }
  void write_by_index(std::size_t /*idx*/, const no_property& /*prop*/) {}

  secret get_edge_bundle(const IndexMap& /*unused*/ = IndexMap()) const { return secret(); }

 protected:
  // All operations do nothing.
  indexed_edge_properties() = default;
  explicit indexed_edge_properties(std::size_t /*unused*/) {}
  [[nodiscard]] std::size_t size() const { return 0; }
  void clear() {}
  void resize(std::size_t /*unused*/) {}
  void reserve(std::size_t /*unused*/) {}

 public:
  void push_back(const edge_push_back_type& /*unused*/) {}
  void move_range(std::size_t /*src_begin*/, std::size_t /*src_end*/, std::size_t /*dest_begin*/) {}

  auto as_range() const { return std::ranges::empty_view<no_property>(); }
};

}  // namespace bagl::indexed_detail

#endif  // BAGL_BAGL_DETAIL_INDEXED_PROPERTIES_H_
