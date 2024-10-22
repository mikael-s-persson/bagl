// Copyright 1999-2001 Jeremy Siek.
// Copyright 2006 Trustees of Indiana University
// Authors: Douglas Gregor and Jeremy Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PROPERTY_MAP_H_
#define BAGL_BAGL_PROPERTY_MAP_H_

#include <cstddef>
#include <iterator>
#include <type_traits>

#include "bagl/has_trait_member.h"

namespace bagl {

//=========================================================================
// property_traits class

BAGL_GRAPH_HAS_TRAIT_MEMBER(value_type, void)

template <typename PA>
struct property_traits {
  using value_type = get_value_type_or_not<PA>;
};

template <typename PA>
using property_traits_value_t = typename property_traits<PA>::value_type;

//=========================================================================
// property_traits specialization for pointers

template <class T>
struct property_traits<T*> {
  using value_type = T;
};
template <class T>
struct property_traits<const T*> {
  using value_type = T;
};

// V must be convertible to T
template <typename T, typename V>
void put(T* pa, std::ptrdiff_t k, V&& val) {
  pa[k] = std::forward<V>(val); // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

template <typename T>
const T& get(const T* pa, std::ptrdiff_t k) {
  return pa[k]; // NOLINT(cppcoreguidelines-pro-bounds-pointer-arithmetic)
}

//=========================================================================
// concept checks for property maps

namespace concepts {

template <typename PMap, typename Key>
concept ReadablePropertyMap = requires(const PMap& pmap, const Key& k) {
  { get(pmap, k) } -> std::convertible_to<property_traits_value_t<PMap>>;
};

template <typename PMap, typename Key>
concept WritablePropertyMap = requires(PMap& pmap, const Key& k) {
  put(pmap, k, std::declval<property_traits_value_t<PMap>>());
};

template <typename PMap, typename Key>
concept ReadWritePropertyMap = ReadablePropertyMap<PMap, Key> && WritablePropertyMap<PMap, Key>;

template <typename PMap, typename Key>
concept LvaluePropertyMap = ReadablePropertyMap<PMap, Key> && requires(PMap& pmap, const Key& k) {
  { pmap[k] } -> std::convertible_to<const property_traits_value_t<PMap>&>;
};

template <typename PMap, typename Key>
concept MutableLvaluePropertyMap = ReadWritePropertyMap<PMap, Key> && requires(PMap& pmap, const Key& k) {
  { pmap[k] } -> std::same_as<property_traits_value_t<PMap>&>;
};

template <typename Func, typename PMap>
concept PropertyCombinator = requires(Func f, const property_traits_value_t<PMap>& value) {
  { f(value, value) } -> std::convertible_to<property_traits_value_t<PMap>>;
};

template <typename Func, typename PMap>
concept PropertyComparator = requires(Func f, const property_traits_value_t<PMap>& value) {
  { f(value, value) } -> std::same_as<bool>;
};

} // namespace concepts

template <typename KeyArchetype, typename ValueArchetype>
struct readable_property_map_archetype {
  using value_type = ValueArchetype;
};
template <typename K, typename V>
V get(const readable_property_map_archetype<K, V>& /*unused*/, const K& /*unused*/) {
  return V{};
}

template <typename KeyArchetype, typename ValueArchetype>
struct writable_property_map_archetype {
  using value_type = ValueArchetype;
};
template <typename K, typename V, std::convertible_to<V> U>
void put(const writable_property_map_archetype<K, V>& /*unused*/, const K& /*unused*/, U&& /*unused*/) {}

template <typename KeyArchetype, typename ValueArchetype>
struct read_write_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype>,
                                           writable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using value_type = ValueArchetype;
};

template <typename KeyArchetype, typename ValueArchetype>
struct lvalue_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using value_type = ValueArchetype;
  const value_type& operator[](const KeyArchetype& /*unused*/) const { return value; }
  value_type value{};
};

template <typename KeyArchetype, typename ValueArchetype>
struct mutable_lvalue_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype>,
                                               writable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using value_type = ValueArchetype;
  value_type& operator[](const KeyArchetype& /*unused*/) const { return value; }
  mutable value_type value{};
};

// A helper class for constructing a property map
// from a class that implements operator[]

template <typename LvaluePropertyMap>
struct put_get_helper {};

template <typename PropertyMap, typename K>
decltype(auto) get(const put_get_helper<PropertyMap>& pa, K&& k) {
  return static_cast<const PropertyMap&>(pa)[std::forward<K>(k)];
}
template <typename PropertyMap, typename K, typename U>
void put(const put_get_helper<PropertyMap>& pa, K&& k, U&& u) {
  static_cast<const PropertyMap&>(pa)[std::forward<K>(k)] = std::forward<U>(u);
}

//=========================================================================
// Adapter to turn an iterator into a property map

template <std::random_access_iterator RAIter, typename IndexMap>
class iterator_property_map : public put_get_helper<iterator_property_map<RAIter, IndexMap>> {
 public:
  using value_type = std::decay_t<decltype(*std::declval<RAIter>())>;

  explicit iterator_property_map(RAIter iter = RAIter(), IndexMap index = IndexMap())
      : iter_(std::move(iter)), index_(std::move(index)) {}

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    return *(iter_ + get(index_, std::forward<Key>(k)));
  }

 protected:
  RAIter iter_;
  IndexMap index_;
};

template <std::random_access_iterator RAIter, typename IndexMap>
class safe_iterator_property_map : public put_get_helper<safe_iterator_property_map<RAIter, IndexMap>> {
 public:
  using value_type = std::decay_t<decltype(*std::declval<RAIter>())>;

  explicit safe_iterator_property_map(RAIter first = RAIter(), std::size_t n = 0,
                             IndexMap index = IndexMap())
      : iter_(std::move(first)), n_(n), index_(std::move(index)) {}

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    assert(get(index_, k) < n_);
    return *(iter_ + get(index_, std::forward<Key>(k)));
  }
  auto size() const { return n_; }

 protected:
  RAIter iter_;
  std::size_t n_ = 0;
  IndexMap index_;
};

//=========================================================================
// An adaptor to turn a Unique Pair Associative Container like std::map or
// std::unordered_map into an Lvalue Property Map.

template <typename UniquePairAssociativeContainer>
class associative_property_map : public put_get_helper<associative_property_map<UniquePairAssociativeContainer>> {
  using C = UniquePairAssociativeContainer;
  using NC = std::remove_cv_t<C>;

 public:
  using value_type = typename NC::mapped_type;
  associative_property_map() = default;
  explicit associative_property_map(C& c) : c_(&c) {}

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    if constexpr (std::is_const_v<C>) {
      return c_->find(std::forward<Key>(k))->second;
    } else {
      return (*c_)[std::forward<Key>(k)];
    }
  }

 private:
  C* c_ = nullptr;
};

//=========================================================================
// A generalized identity property map
template <typename T = void>
struct typed_identity_property_map : public put_get_helper<typed_identity_property_map<T>> {
  using value_type = T;

  template <typename U>
  decltype(auto) operator[](U&& p) const {
    return std::forward<U>(p);
  }
};

//=========================================================================
// A property map that applies the identity function to integers
using identity_property_map = typed_identity_property_map<std::size_t>;

//======== Data-member property-map ==========

// This property-map class can be used to map an object (by reference) to one
// of its data-members (by reference).
template <typename T, typename PropertyType>
class data_member_property_map : public put_get_helper<data_member_property_map<T, PropertyType>> {
 public:
  using member_ptr_type = T PropertyType::*;
  using self = data_member_property_map<T, PropertyType>;

 private:
  member_ptr_type mem_ptr_ = nullptr;

 public:
  using value_type = T;

  explicit data_member_property_map(member_ptr_type aMemPtr) : mem_ptr_(aMemPtr) {}
  data_member_property_map() = default;
  template <typename OtherProp>
  decltype(auto) operator[](OtherProp& p) const {
    return p.*mem_ptr_;
  }
};

//======== Composite property-map ==========

// This property-map allows for the composition of two property-maps.
// The resulting property-map has the value-type of the output-map and
// the key-type of the input-map.
template <typename OutputMap, typename InputMap>
class composite_property_map : public put_get_helper<composite_property_map<OutputMap, InputMap>> {
 private:
  OutputMap prop_out_;
  InputMap prop_in_;

 public:
  using value_type = property_traits_value_t<OutputMap>;

  explicit composite_property_map(OutputMap prop_out, InputMap prop_in = {}) : prop_out_(prop_out), prop_in_(prop_in) {}
  composite_property_map() = default;

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    return prop_out_[prop_in_[std::forward<Key>(k)]];
  }
};

//======== Property-map reference-wrapper ==========

// This property-map is like std::reference_wrapper for property-maps.
template <typename UnderlyingMap>
class property_map_ref : public put_get_helper<property_map_ref<UnderlyingMap>> {
 private:
  UnderlyingMap* prop_ = nullptr;

 public:
  using value_type = property_traits_value_t<UnderlyingMap>;

  explicit property_map_ref(UnderlyingMap& prop) : prop_(&prop) {}
  property_map_ref() = default;

  template <typename Key>
  decltype(auto) operator[](Key&& k) const {
    return (*prop_)[std::forward<Key>(k)];
  }
};

//=========================================================================
// A property map that does not do anything, for
// when you have to supply a property map, but don't need it.
namespace property_map_detail {
struct dummy_pmap_proxy {
  template <typename T>
  dummy_pmap_proxy& operator=(const T& /*unused*/) {
    return *this;
  }
  template <typename T>
  operator T() { return T{}; } // NOLINT
};
}  // namespace property_map_detail

class dummy_property_map : public put_get_helper<dummy_property_map> {
 public:
  using value_type = int;

  template <typename T>
  property_map_detail::dummy_pmap_proxy& operator[](T&& /*unused*/) const {
    return dummy_lvalue_;
  }
 private:
  mutable property_map_detail::dummy_pmap_proxy dummy_lvalue_{};
};

// Convert a Readable property map into a function object
template <typename PropMap>
class property_map_function {
  PropMap pm_;

 public:
  explicit property_map_function(const PropMap& pm) : pm_(pm) {}
  template <typename T>
  decltype(auto) operator()(T&& k) const {
    return get(pm_, std::forward<T>(k));
  }
};

template <typename PropMap>
property_map_function<PropMap> make_property_map_function(const PropMap& pm) {
  return property_map_function<PropMap>(pm);
}

}  // namespace bagl

#include "bagl/single_property_map.h"
#include "bagl/vector_property_map.h"

#endif  // BAGL_BAGL_PROPERTY_MAP_H_
