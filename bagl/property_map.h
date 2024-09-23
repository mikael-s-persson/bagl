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

BAGL_GRAPH_HAS_TRAIT_MEMBER(key_type, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(value_type, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(reference, void)
BAGL_GRAPH_HAS_TRAIT_MEMBER(category, void)

template <class PA>
constexpr bool is_property_map_v = has_key_type_v<PA>&& has_value_type_v<PA>&& has_reference_v<PA>&& has_category_v<PA>;

template <typename PA>
struct default_property_traits {
  using key_type = get_key_type_or_not<PA>;
  using value_type = get_value_type_or_not<PA>;
  using reference = get_reference_or_not<PA>;
  using category = get_category_or_not<PA>;
};

struct null_property_traits {};

template <typename PA>
struct property_traits : std::conditional_t<is_property_map_v<PA>, default_property_traits<PA>, null_property_traits> {
};

template <typename PA>
using property_traits_key_t = typename property_traits<PA>::key_type;
template <typename PA>
using property_traits_value_t = typename property_traits<PA>::value_type;
template <typename PA>
using property_traits_reference_t = typename property_traits<PA>::reference;
template <typename PA>
using property_traits_category_t = typename property_traits<PA>::category;

//=========================================================================
// property_traits category tags

namespace property_map_detail {
enum ePropertyMapID { READABLE_PA, WRITABLE_PA, READ_WRITE_PA, LVALUE_PA, OP_BRACKET_PA, RAND_ACCESS_ITER_PA, LAST_PA };
}
struct readable_property_map_tag {
  static constexpr auto id = property_map_detail::READABLE_PA;
};
struct writable_property_map_tag {
  static constexpr auto id = property_map_detail::WRITABLE_PA;
};
struct read_write_property_map_tag : public readable_property_map_tag, public writable_property_map_tag {
  static constexpr auto id = property_map_detail::READ_WRITE_PA;
};

struct lvalue_property_map_tag : public read_write_property_map_tag {
  static constexpr auto id = property_map_detail::LVALUE_PA;
};

//=========================================================================
// property_traits specialization for pointers

template <class T>
struct property_traits<T*> {
  using value_type = T;
  using reference = value_type&;
  using key_type = std::ptrdiff_t;
  using category = lvalue_property_map_tag;
};
template <class T>
struct property_traits<const T*> {
  using value_type = T;
  using reference = const value_type&;
  using key_type = std::ptrdiff_t;
  using category = lvalue_property_map_tag;
};

// V must be convertible to T
template <typename T, typename V>
void put(T* pa, std::ptrdiff_t k, V&& val) {
  pa[k] = std::forward<V>(val);
}

template <typename T>
const T& get(const T* pa, std::ptrdiff_t k) {
  return pa[k];
}

//=========================================================================
// concept checks for property maps

namespace concepts {

template <typename PMap, typename Key>
concept ReadablePropertyMap = std::convertible_to<property_traits_category_t<PMap>, readable_property_map_tag> &&
    requires(const PMap& pmap, const Key& k) {
  { get(pmap, k) } -> std::convertible_to<property_traits_value_t<PMap>>;
};

template <typename PMap, typename Key>
concept WritablePropertyMap = std::convertible_to<property_traits_category_t<PMap>, writable_property_map_tag> &&
    requires(PMap& pmap, const Key& k) {
  put(pmap, k, std::declval<property_traits_value_t<PMap>>());
};

template <typename PMap, typename Key>
concept ReadWritePropertyMap = std::convertible_to<property_traits_category_t<PMap>, read_write_property_map_tag> &&
    ReadablePropertyMap<PMap, Key> && WritablePropertyMap<PMap, Key>;

template <typename PMap, typename Key>
concept LvaluePropertyMap = std::convertible_to<property_traits_category_t<PMap>, lvalue_property_map_tag> &&
    ReadablePropertyMap<PMap, Key> && requires(PMap& pmap, const Key& k) {
  { pmap[k] } -> std::same_as<property_traits_value_t<PMap>&> || std::same_as<const property_traits_value_t<PMap>&>;
};

template <typename PMap, typename Key>
concept Mutable_LvaluePropertyMap = std::convertible_to<property_traits_category_t<PMap>, lvalue_property_map_tag> &&
    ReadWritePropertyMap<PMap, Key> && requires(PMap& pmap, const Key& k) {
  { pmap[k] } -> std::same_as<property_traits_value_t<PMap>&>;
};

} // namespace concepts

template <typename KeyArchetype, typename ValueArchetype>
struct readable_property_map_archetype {
  using key_type = KeyArchetype;
  using value_type = ValueArchetype;
  using reference = const ValueArchetype&;
  using category = readable_property_map_tag;
};
template <typename K, typename V>
const V& get(const readable_property_map_archetype<K, V>&, const K&) {
  return std::declval<V>();
}

template <typename KeyArchetype, typename ValueArchetype>
struct writable_property_map_archetype {
  using key_type = KeyArchetype;
  using value_type = ValueArchetype;
  using reference = void;
  using category = writable_property_map_tag;
};
template <typename K, typename V, std::convertible_to<V> U>
void put(const writable_property_map_archetype<K, V>&, const K&, U&&) {}

template <typename KeyArchetype, typename ValueArchetype>
struct read_write_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype>,
                                           writable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using key_type = KeyArchetype;
  using value_type = ValueArchetype;
  using reference = const ValueArchetype&;
  using category = read_write_property_map_tag;
};

template <typename KeyArchetype, typename ValueArchetype>
struct lvalue_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using key_type = KeyArchetype;
  using value_type = ValueArchetype;
  using reference = const ValueArchetype&;
  using category = lvalue_property_map_tag;
  const value_type& operator[](const key_type&) const { return std::declval<value_type>(); }
};

template <typename KeyArchetype, typename ValueArchetype>
struct mutable_lvalue_property_map_archetype : readable_property_map_archetype<KeyArchetype, ValueArchetype>,
                                               writable_property_map_archetype<KeyArchetype, ValueArchetype> {
  using key_type = KeyArchetype;
  using value_type = ValueArchetype;
  using reference = ValueArchetype&;
  using category = lvalue_property_map_tag;
  value_type& operator[](const key_type&) const { return std::declval<value_type>(); }
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
requires ReadablePropertyMap<IndexMap, property_traits_key_t<IndexMap>>
class iterator_property_map : public put_get_helper<iterator_property_map<RandomAccessIterator, IndexMap>> {
 public:
  using key_type = property_traits_key_t<IndexMap>;
  using reference = decltype(*std::declval<RandomAccessIterator>());
  using value_type = std::decay_t<reference>;
  using category = lvalue_property_map_tag;

  iterator_property_map(RandomAccessIterator iter = RandomAccessIterator(), IndexMap index = IndexMap())
      : iter_(std::move(iter)), index_(std::move(index)) {}
  reference operator[](const key_type& k) const { return *(iter_ + get(index_, k)); }

 protected:
  RandomAccessIterator iter_;
  IndexMap index_;
};

template <std::random_access_iterator RAIter, typename IndexMap>
requires ReadablePropertyMap<IndexMap, property_traits_key_t<IndexMap>>
auto make_iterator_property_map(RAIter iter, IndexMap id) {
  return iterator_property_map<RAIter, IndexMap>(std::move(iter), std::move(id));
}

template <std::random_access_iterator RandomAccessIterator, typename IndexMap>
requires ReadablePropertyMap<IndexMap, property_traits_key_t<IndexMap>>
class safe_iterator_property_map : public put_get_helper<safe_iterator_property_map<RandomAccessIterator, IndexMap>> {
 public:
  using key_type = property_traits_key_t<IndexMap>;
  using reference = decltype(*std::declval<RandomAccessIterator>());
  using value_type = std::decay_t<reference>;
  using category = lvalue_property_map_tag;

  safe_iterator_property_map(RandomAccessIterator first, property_traits_value_t<IndexMap> n = 0,
                             IndexMap index = IndexMap())
      : iter_(std::move(first)), n_(n), index_(std::move(index)) {}
  safe_iterator_property_map() {}
  reference operator[](const key_type& k) const {
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
requires ReadablePropertyMap<IndexMap, property_traits_key_t<IndexMap>>
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
  using key_type = typename C::key_type;
  using reference = decltype(std::declval<C>()[std::declval<key_type>()]);
  using value_type = std::decay_t<reference>;
  using category = lvalue_property_map_tag;
  associative_property_map() : c_(nullptr) {}
  explicit associative_property_map(C& c) : c_(&c) {}
  reference operator[](const key_type& k) const { return (*c_)[k]; }

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
  using key_type = typename C::key_type;
  using reference = decltype(std::declval<const C&>().find(std::declval<key_type>())->second);
  using value_type = std::decay_t<reference>;
  using category = lvalue_property_map_tag;
  const_associative_property_map() : c_(nullptr) {}
  const_associative_property_map(const C& c) : c_(&c) {}
  reference operator[](const key_type& k) const { return c_->find(k)->second; }

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
template <typename ValueType, typename KeyType = void>
class static_property_map : public put_get_helper<static_property_map<ValueType>> {
  ValueType value_;

 public:
  using key_type = KeyType;
  using value_type = ValueType;
  using reference = ValueType;
  using category = readable_property_map_tag;
  explicit static_property_map(ValueType v) : value_(std::move(v)) {}

  template <typename T>
  reference operator[](T&& /*unused*/) const {
    return value_;
  }
};

template <typename KeyType, typename U>
auto make_static_property_map(U&& v) {
  return static_property_map<std::decay_t<U>, KeyType>(std::forward<U>(v));
}

//=========================================================================
// A property map that always returns a reference to the same object.
//
template <typename KeyType, typename ValueType>
class ref_property_map : public put_get_helper<ref_property_map<KeyType, ValueType>> {
  ValueType* value_;

 public:
  using key_type = KeyType;
  using value_type = ValueType;
  using reference = ValueType&;
  using category = lvalue_property_map_tag;
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
  using key_type = T;
  using value_type = T;
  using reference = T;
  using category = readable_property_map_tag;

  reference operator[](key_type v) const { return v; }
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
  dummy_pmap_reference& operator=(const T&) {
    return *this;
  }
  operator int() { return 0; }
};
}  // namespace property_map_detail

class dummy_property_map : public put_get_helper<dummy_property_map> {
 public:
  using key_type = void;
  using value_type = int;
  using reference = property_map_detail::dummy_pmap_reference;
  using category = read_write_property_map_tag;
  template <typename T>
  reference operator[](T&& /*unused*/) const {
    return reference();
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
