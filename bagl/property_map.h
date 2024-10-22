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
  using reference = value_type&;
};
template <class T>
struct property_traits<const T*> {
  using value_type = T;
  using reference = const value_type&;
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
  using reference = const ValueArchetype&;
};
template <typename K, typename V>
V get(const readable_property_map_archetype<K, V>& /*unused*/, const K& /*unused*/) {
  return V{};
}

template <typename KeyArchetype, typename ValueArchetype>
struct writable_property_map_archetype {
  using value_type = ValueArchetype;
  using reference = void;
};
template <typename K, typename V, std::convertible_to<V> U>
void put(const writable_property_map_archetype<K, V>& /*unused*/, const K& /*unused*/, U&& /*unused*/) {}

template <typename KeyArchetype, typename ValueArchetype>
struct read_write_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype>,
                                           writable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using value_type = ValueArchetype;
  using reference = const ValueArchetype&;
};

template <typename KeyArchetype, typename ValueArchetype>
struct lvalue_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using value_type = ValueArchetype;
  using reference = const ValueArchetype&;
  const value_type& operator[](const KeyArchetype& /*unused*/) const { return value; }
  value_type value{};
};

template <typename KeyArchetype, typename ValueArchetype>
struct mutable_lvalue_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype>,
                                               writable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using value_type = ValueArchetype;
  using reference = ValueArchetype&;
  value_type& operator[](const KeyArchetype& /*unused*/) const { return value; }
  mutable value_type value{};
};

template <typename T>
struct typed_identity_property_map;

// A helper class for constructing a property map
// from a class that implements operator[]

template <typename LvaluePropertyMap>
struct put_get_helper {};

template <typename PropertyMap, typename K>
decltype(auto) get(const put_get_helper<PropertyMap>& pa, const K& k) {
  return static_cast<const PropertyMap&>(pa)[k];
}
template <typename PropertyMap, typename K, typename U>
void put(const put_get_helper<PropertyMap>& pa, const K& k, U&& u) {
  static_cast<const PropertyMap&>(pa)[k] = std::forward<U>(u);
}

//=========================================================================
// Adapter to turn a RandomAccessIterator into a property map

template <std::random_access_iterator RandomAccessIterator, typename IndexMap>
class iterator_property_map : public put_get_helper<iterator_property_map<RandomAccessIterator, IndexMap>> {
 public:
  using reference = decltype(*std::declval<RandomAccessIterator>());
  using value_type = std::decay_t<reference>;

  explicit iterator_property_map(RandomAccessIterator iter = RandomAccessIterator(), IndexMap index = IndexMap())
      : iter_(std::move(iter)), index_(std::move(index)) {}

  template <typename Key>
  reference operator[](const Key& k) const {
    return *(iter_ + get(index_, k));
  }

 protected:
  RandomAccessIterator iter_;
  IndexMap index_;
};

template <std::random_access_iterator RAIter, typename IndexMap>
auto make_iterator_property_map(RAIter iter, IndexMap id) {
  return iterator_property_map<RAIter, IndexMap>(std::move(iter), std::move(id));
}

template <std::random_access_iterator RandomAccessIterator, typename IndexMap>
class safe_iterator_property_map : public put_get_helper<safe_iterator_property_map<RandomAccessIterator, IndexMap>> {
 public:
  using reference = decltype(*std::declval<RandomAccessIterator>());
  using value_type = std::decay_t<reference>;

  explicit safe_iterator_property_map(RandomAccessIterator first, property_traits_value_t<IndexMap> n = 0,
                             IndexMap index = IndexMap())
      : iter_(std::move(first)), n_(n), index_(std::move(index)) {}
  safe_iterator_property_map() = default;

  template <typename Key>
  reference operator[](const Key& k) const {
    assert(get(index_, k) < n_);
    return *(iter_ + get(index_, k));
  }
  auto size() const { return n_; }

 protected:
  RandomAccessIterator iter_;
  property_traits_value_t<IndexMap> n_;
  IndexMap index_;
};

template <std::random_access_iterator RAIter, typename IndexMap>
auto make_safe_iterator_property_map(RAIter iter, property_traits_value_t<IndexMap> n, IndexMap id) {
  return safe_iterator_property_map<RAIter, IndexMap>(std::move(iter), n, std::move(id));
}

//=========================================================================
// An adaptor to turn a Unique Pair Associative Container like std::map or
// std::unordered_map into an Lvalue Property Map.

template <typename UniquePairAssociativeContainer>
class associative_property_map : public put_get_helper<associative_property_map<UniquePairAssociativeContainer>> {
  using C = UniquePairAssociativeContainer;

 public:
  using value_type = typename C::mapped_type;
  using reference = value_type&;
  associative_property_map() : c_(nullptr) {}
  explicit associative_property_map(C& c) : c_(&c) {}

  template <typename Key>
  reference operator[](const Key& k) const {
    return (*c_)[k];
  }

 private:
  C* c_;
};

template <typename UniquePairAssociativeContainer>
auto make_assoc_property_map(UniquePairAssociativeContainer& c) {
  return associative_property_map<UniquePairAssociativeContainer>(c);
}

template <typename UniquePairAssociativeContainer>
class const_associative_property_map
    : public put_get_helper<const_associative_property_map<UniquePairAssociativeContainer>> {
  using C = UniquePairAssociativeContainer;

 public:
  using value_type = typename C::mapped_type;
  using reference = const value_type&;
  const_associative_property_map() : c_(nullptr) {}
  explicit const_associative_property_map(const C& c) : c_(&c) {}

  template <typename Key>
  reference operator[](const Key& k) const {
    return c_->find(k)->second;
  }

 private:
  const C* c_;
};

template <typename UniquePairAssociativeContainer>
auto make_assoc_property_map(const UniquePairAssociativeContainer& c) {
  return const_associative_property_map<UniquePairAssociativeContainer>(c);
}

//=========================================================================
// A property map that always returns the same object by value.
//
template <typename ValueType>
class static_property_map : public put_get_helper<static_property_map<ValueType>> {
  ValueType value_;

 public:
  using value_type = ValueType;
  using reference = const ValueType&;
  explicit static_property_map(ValueType v) : value_(std::move(v)) {}

  template <typename T>
  reference operator[](T&& /*unused*/) const {
    return value_;
  }
};

template <typename U>
auto make_static_property_map(U&& v) {
  return static_property_map<std::decay_t<U>>(std::forward<U>(v));
}

//=========================================================================
// A property map that always returns a reference to the same object.
//
template <typename ValueType>
class ref_property_map : public put_get_helper<ref_property_map<ValueType>> {
  ValueType* value_;

 public:
  using value_type = ValueType;
  using reference = ValueType&;
  explicit ref_property_map(ValueType& v) : value_(&v) {}

  template <typename T>
  ValueType& operator[](T&& /*unused*/) const {
    return *value_;
  }
};

//=========================================================================
// A generalized identity property map
template <typename T>
struct typed_identity_property_map : public put_get_helper<typed_identity_property_map<T>> {
  using value_type = T;
  using reference = T;

  T operator[](T v) const { return v; }
};

//=========================================================================
// A property map that applies the identity function to integers
using identity_property_map = typed_identity_property_map<std::size_t>;

//=========================================================================
// A property map that does not do anything, for
// when you have to supply a property map, but don't need it.
namespace property_map_detail {
struct dummy_pmap_reference {
  template <class T>
  dummy_pmap_reference& operator=(const T& /*unused*/) {
    return *this;
  }
  operator int() { return 0; } // NOLINT
};
}  // namespace property_map_detail

class dummy_property_map : public put_get_helper<dummy_property_map> {
 public:
  using value_type = int;
  using reference = property_map_detail::dummy_pmap_reference;
  template <typename T>
  reference operator[](T&& /*unused*/) const {
    return {};
  }
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

#include "bagl/vector_property_map.h"

#endif  // BAGL_BAGL_PROPERTY_MAP_H_
