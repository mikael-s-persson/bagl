// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_DETAIL_CONTAINER_GENERATORS_H_
#define BAGL_BAGL_DETAIL_CONTAINER_GENERATORS_H_

#include <cstdint>
#include <limits>
#include <list>
#include <map>
#include <ranges>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

#include "bagl/graph_selectors.h"

namespace bagl {

namespace container_detail {

/*************************************************************************
 *        pool-container values: validity and retrieval
 * **********************************************************************/

struct hole_desc {
  std::size_t value;
  explicit hole_desc(std::size_t aValue) : value(aValue) {}

  hole_desc() : hole_desc((std::numeric_limits<std::size_t>::max)()) {}

  bool operator==(hole_desc rhs) const { return (rhs.value == this->value); }
  bool operator!=(hole_desc rhs) const { return (rhs.value != this->value); }
};

template <typename ValueType>
bool is_elem_valid(const ValueType& /*unused*/) {
  return true;
}

template <typename ValueType>
bool is_elem_valid(const std::variant<ValueType, hole_desc>& vp) {
  return (vp.index() == 0);
}

template <typename Container>
struct is_not_hole {
  const Container* p_cont;
  explicit is_not_hole(const Container* cont) : p_cont(cont) {}
  is_not_hole() : is_not_hole(nullptr) {}
  template <typename Desc>
  bool operator()(Desc d) const {
    return ((*p_cont)[d].index() == 0);
  }
};

template <typename Container>
struct edge_is_not_hole {
  const Container* p_cont;
  explicit edge_is_not_hole(const Container* cont) : p_cont(cont) {}
  edge_is_not_hole() : edge_is_not_hole(nullptr) {}
  template <typename Desc>
  bool operator()(Desc e) const {
    return ((*p_cont)[e.edge_id].index() == 0);
  }
};

template <typename ValueType>
ValueType& get_value(ValueType& vp) {
  return vp;
}

template <typename ValueType>
ValueType& get_value(std::variant<ValueType, hole_desc>& vp) {
  return std::get<ValueType>(vp);
}

template <typename ValueType>
const ValueType& get_value(const std::variant<ValueType, hole_desc>& vp) {
  return std::get<ValueType>(vp);
}

// A pooled vector is basically an arena allocator implemented as a std::vector
// of elements that can be emptied and later reused.
template <typename ValueType>
struct pooled_vector {
  using value_type = ValueType;
  using size_type = std::size_t;
  using stored_type = std::variant<ValueType, hole_desc>;
  using container_type = std::vector<stored_type>;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

  std::vector<stored_type> m_data;
  hole_desc m_first_hole;
  std::size_t m_num_elements{0};

  pooled_vector() : m_data() {}

  void swap(pooled_vector& rhs) {
    using std::swap;
    m_data.swap(rhs.m_data);
    swap(m_first_hole, rhs.m_first_hole);
    swap(m_num_elements, rhs.m_num_elements);
  }

  std::size_t size() const { return m_num_elements; }
  std::size_t capacity() const { return m_data.capacity(); }
  void clear() {
    m_data.clear();
    m_first_hole = hole_desc();
    m_num_elements = 0;
  }
  iterator begin() { return m_data.begin(); }
  const_iterator begin() const { return m_data.begin(); }
  iterator end() { return m_data.end(); }
  const_iterator end() const { return m_data.end(); }
};

template <typename DescType>
struct null_desc {
  static constexpr DescType value() { return DescType{}; }
};

template <typename DescType>
struct null_desc<DescType*> {
  static constexpr DescType* value() { return nullptr; }
};

template <>
struct null_desc<std::size_t> {
  static constexpr std::size_t value() { return (std::numeric_limits<std::size_t>::max)(); }
};

template <typename Key, typename Value>
struct null_desc<std::pair<Key, Value*>> {
  static constexpr std::pair<Key, Value*> value() { return {null_desc<Key>::value(), nullptr}; }
};

// vector/pool containers.
inline bool desc_less_than(std::size_t lhs, std::size_t rhs) { return (lhs < rhs); }

// node containers.
template <typename Iter>
bool desc_less_than(Iter lhs, Iter rhs) {
  return (&(*lhs) < &(*rhs));
}

// hash containers.
template <typename Key, typename Value>
bool desc_less_than(std::pair<Key, Value*> lhs, std::pair<Key, Value*> rhs) {
  return (desc_less_than(lhs.first, rhs.first) || ((lhs.first == rhs.first) && (lhs.second < rhs.second)));
}

/*************************************************************************
 *                      indexable descriptor
 * **********************************************************************/

// This is just a simple wrappper for a descriptor that gives it ordering and hashing.
// This is so they can be used as keys in the out-edge containers, we could provide
// hash/compare functors, but this is just simpler since we don't have to have the hash/comp
// template parameters everywhere, and we template on key types already.
template <typename RawDesc>
struct indexable_desc {
  RawDesc value;

  // Implicitly convertible.
  template <typename... Args>
  indexable_desc(Args&&... args) : value(std::forward<Args>(args)...) {}
  operator RawDesc&() { return value; }
  operator const RawDesc&() const { return value; }

  bool operator==(const indexable_desc& rhs) const { return value == rhs.value; }
  bool operator!=(const indexable_desc& rhs) const { return value != rhs.value; }
  bool operator<(const indexable_desc& rhs) const { return desc_less_than(value, rhs.value); }
  bool operator>(const indexable_desc& rhs) const { return desc_less_than(rhs.value, value); }
  bool operator<=(const indexable_desc& rhs) const { return !desc_less_than(rhs.value, value); }
  bool operator>=(const indexable_desc& rhs) const { return !desc_less_than(value, rhs.value); }
};

template <typename RawDesc>
using indexable_desc_t = std::conditional_t<std::is_integral_v<RawDesc>, RawDesc, indexable_desc<RawDesc>>;

}  // namespace container_detail

template <typename Selector, typename ValueType>
struct container_gen {};

template <typename ValueType>
struct container_gen<vec_s, ValueType> {
  using type = std::vector<ValueType>;
};

template <typename ValueType>
struct container_gen<pool_s, ValueType> {
  using type = container_detail::pooled_vector<ValueType>;
};

template <typename ValueType>
struct container_gen<list_s, ValueType> {
  using type = std::list<ValueType>;
};

template <typename Key, typename Value>
struct container_gen<set_s, std::pair<Key, Value>> {
  using type = std::map<container_detail::indexable_desc_t<Key>, Value>;
};

template <typename Key, typename Value>
struct container_gen<multiset_s, std::pair<Key, Value>> {
  using type = std::multimap<container_detail::indexable_desc_t<Key>, Value>;
};

template <typename Key, typename Value>
struct container_gen<unordered_set_s, std::pair<Key, Value>> {
  using type = std::unordered_map<container_detail::indexable_desc_t<Key>, Value>;
};

template <typename Key, typename Value>
struct container_gen<unordered_multiset_s, std::pair<Key, Value>> {
  using type = std::unordered_multimap<container_detail::indexable_desc_t<Key>, Value>;
};

namespace container_detail {

template <typename Selector>
constexpr bool is_random_access_v = std::is_same_v<Selector, vec_s>;

/*************************************************************************
 *                      edge descriptors
 * **********************************************************************/

template <typename Vertex, typename EdgeRawDesc>
struct edge_desc {
  using edge_id_type = EdgeRawDesc;
  using source_descriptor = Vertex;

  Vertex source;
  edge_id_type edge_id;

  explicit edge_desc(Vertex src, edge_id_type edgeid = null_desc<edge_id_type>::value())
      : source(src), edge_id(edgeid) {}

  edge_desc() : edge_desc(null_desc<Vertex>::value()) {}

  static edge_desc null_value() { return edge_desc(); }
};

template <typename Vertex, typename EdgeRawDesc>
bool operator==(const edge_desc<Vertex, EdgeRawDesc>& lhs, const edge_desc<Vertex, EdgeRawDesc>& rhs) {
  if ((lhs.source == null_desc<Vertex>::value()) && (rhs.source == null_desc<Vertex>::value())) {
    return true;
  }
  return (lhs.source == rhs.source) && (lhs.edge_id == rhs.edge_id);
}
template <typename Vertex, typename EdgeRawDesc>
bool operator!=(const edge_desc<Vertex, EdgeRawDesc>& lhs, const edge_desc<Vertex, EdgeRawDesc>& rhs) {
  if ((lhs.source == null_desc<Vertex>::value()) && (rhs.source == null_desc<Vertex>::value())) {
    return false;
  }
  return (lhs.source != rhs.source) || (lhs.edge_id != rhs.edge_id);
}
template <typename Vertex, typename EdgeRawDesc>
bool operator>(const edge_desc<Vertex, EdgeRawDesc>& lhs, const edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (desc_less_than(rhs.source, lhs.source) ||
          ((lhs.source == rhs.source) && (desc_less_than(rhs.edge_id, lhs.edge_id))));
}
template <typename Vertex, typename EdgeRawDesc>
bool operator>=(const edge_desc<Vertex, EdgeRawDesc>& lhs, const edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (desc_less_than(rhs.source, lhs.source) ||
          ((lhs.source == rhs.source) && (!desc_less_than(lhs.edge_id, rhs.edge_id))));
}

template <typename Vertex, typename EdgeRawDesc>
bool operator<(const edge_desc<Vertex, EdgeRawDesc>& lhs, const edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (desc_less_than(lhs.source, rhs.source) ||
          ((lhs.source == rhs.source) && (desc_less_than(lhs.edge_id, rhs.edge_id))));
}
template <typename Vertex, typename EdgeRawDesc>
bool operator<=(const edge_desc<Vertex, EdgeRawDesc>& lhs, const edge_desc<Vertex, EdgeRawDesc>& rhs) {
  return (desc_less_than(lhs.source, rhs.source) ||
          ((lhs.source == rhs.source) && (!desc_less_than(rhs.edge_id, lhs.edge_id))));
}

template <typename EdgeDesc>
struct undir_edge_desc : EdgeDesc {
  bool is_reversed = false;

  explicit undir_edge_desc(const EdgeDesc& base, bool isreversed = false) : EdgeDesc(base), is_reversed(isreversed) {}

  undir_edge_desc() : undir_edge_desc(EdgeDesc()) {}

  static undir_edge_desc<EdgeDesc> null_value() { return undir_edge_desc<EdgeDesc>(EdgeDesc::null_value()); }
};

/*************************************************************************
 *        iterator / descriptor translation functions for all container types
 * **********************************************************************/

template <typename ValueType>
std::size_t iterator_to_desc(const std::vector<ValueType>& c, typename std::vector<ValueType>::iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t iterator_to_desc(const pooled_vector<ValueType>& c, typename pooled_vector<ValueType>::iterator it) {
  return it - c.begin();
}

template <typename Key, typename Value>
std::pair<Key, Value*> iterator_to_desc(const std::unordered_map<Key, Value>& c,
                                        typename std::unordered_map<Key, Value>::iterator it) {
  return {it->first, &it->second};
}

template <typename Key, typename Value>
std::pair<Key, Value*> iterator_to_desc(const std::unordered_multimap<Key, Value>& c,
                                        typename std::unordered_multimap<Key, Value>::iterator it) {
  return {it->first, &it->second};
}

template <typename Container, typename Iter>
Iter iterator_to_desc(const Container& /*unused*/, Iter it) {
  return it;
}

/*************************************************************************
 *        descriptor / value translation functions for all container types
 * **********************************************************************/

template <typename ValueType>
ValueType& desc_to_value(std::vector<ValueType>& c, std::size_t d) {
  return c[d];
}
template <typename ValueType>
const ValueType& desc_to_value(const std::vector<ValueType>& c, std::size_t d) {
  return c[d];
}

template <typename ValueType>
ValueType& desc_to_value(pooled_vector<ValueType>& c, std::size_t d) {
  return std::get<ValueType>(c.m_data[d]);
}
template <typename ValueType>
const ValueType& desc_to_value(const pooled_vector<ValueType>& c, std::size_t d) {
  return std::get<ValueType>(c.m_data[d]);
}

template <typename Container>
auto& desc_to_value(const Container& /*unused*/, typename Container::iterator it) {
  return *it;
}

template <typename Key, typename Value>
std::pair<Key, Value&> desc_to_value(std::unordered_map<Key, Value>& /*unused*/, std::pair<Key, Value*> p) {
  return {p.first, *p.second};
}
template <typename Key, typename Value>
std::pair<Key, const Value&> desc_to_value(const std::unordered_map<Key, Value>& /*unused*/, std::pair<Key, Value*> p) {
  return {p.first, *p.second};
}

template <typename Key, typename Value>
std::pair<Key, Value&> desc_to_value(std::unordered_multimap<Key, Value>& /*unused*/, std::pair<Key, Value*> p) {
  return {p.first, *p.second};
}
template <typename Key, typename Value>
std::pair<Key, const Value&> desc_to_value(const std::unordered_multimap<Key, Value>& /*unused*/,
                                           std::pair<Key, Value*> p) {
  return {p.first, *p.second};
}

template <typename Container>
struct select_descriptor {
  using type = typename Container::iterator;
};

template <typename Key, typename Value>
struct select_descriptor<std::unordered_map<Key, Value>> {
  using type = std::pair<Key, Value*>;
};

template <typename Key, typename Value>
struct select_descriptor<std::unordered_multimap<Key, Value>> {
  using type = std::pair<Key, Value*>;
};

template <typename Container>
struct select_descriptor<Container*> : select_descriptor<Container> {};

template <typename ValueType>
struct select_descriptor<std::vector<ValueType>> {
  using type = std::size_t;
};

template <typename ValueType>
struct select_descriptor<pooled_vector<ValueType>> {
  using type = std::size_t;
};

// Dummy "ignore" output iterator that is used in remove-branch functions.
struct ignore_output_iter {
  struct value_type {
    template <typename T>
    value_type& operator=(const T& /*unused*/) {
      return *this;
    }
    template <typename T>
    value_type& operator=(T&& /*unused*/) {
      return *this;
    }
  };
  using difference_type = std::ptrdiff_t;
  using pointer = value_type*;
  using reference = value_type&;
  using iterator_category = std::bidirectional_iterator_tag;

  mutable value_type empty_obj;

  ignore_output_iter& operator++() { return *this; }
  ignore_output_iter operator++(int) { return *this; }
  ignore_output_iter& operator--() { return *this; }
  ignore_output_iter operator--(int) { return *this; }

  reference operator*() const { return empty_obj; }

  // Never equal, aka infinite range.
  bool operator==(const ignore_output_iter& /*unused*/) const { return false; }
  bool operator!=(const ignore_output_iter& /*unused*/) const { return true; }
};

inline auto ignore_output_range() { return std::ranges::subrange<ignore_output_iter, ignore_output_iter>(); }

// We use Jon Maiga's implementation from
// http://jonkagstrom.com/mx3/mx3_rev2.html
inline std::uint64_t desc_get_hash(std::uint64_t x) {
  const std::uint64_t m = 0xe9846af9b1a615d;
  x ^= x >> 32;
  x *= m;
  x ^= x >> 32;
  x *= m;
  x ^= x >> 28;
  return x;
}

template <std::indirectly_readable Iter>
std::uint64_t desc_get_hash(Iter it) {
  return desc_get_hash(reinterpret_cast<std::uint64_t>(&(*it)));
}

// Not really needed.
template <typename Key, typename Value>
std::uint64_t desc_get_hash(const std::pair<Key, Value*>& p) {
  return desc_get_hash(p.first);
}

template <typename... Args>
std::uint64_t desc_combined_hash(const Args&... args) {
  std::uint64_t seed = 0;
  ((seed = desc_get_hash(seed + 0x9e3779b9 + desc_get_hash(args))), ...);
  return seed;
}

template <typename Vertex, typename EdgeRawDesc>
std::uint64_t desc_get_hash(const bagl::container_detail::edge_desc<Vertex, EdgeRawDesc>& x) {
  return bagl::container_detail::desc_combined_hash(x.source, x.edge_id);
}

template <typename EdgeRawDesc>
std::uint64_t desc_get_hash(const bagl::container_detail::undir_edge_desc<EdgeRawDesc>& x) {
  return bagl::container_detail::desc_combined_hash(x.is_reversed ? 1 : 0, static_cast<const EdgeRawDesc&>(x));
}

struct desc_hasher {
  template <typename T>
  std::uint64_t operator()(const T& v) const {
    return desc_get_hash(v);
  }
};

}  // namespace container_detail

}  // namespace bagl

template <typename Vertex, typename EdgeRawDesc>
struct std::hash<bagl::container_detail::edge_desc<Vertex, EdgeRawDesc>> {
  std::size_t operator()(const bagl::container_detail::edge_desc<Vertex, EdgeRawDesc>& x) const {
    return bagl::container_detail::desc_get_hash(x);
  }
};

template <typename EdgeRawDesc>
struct std::hash<bagl::container_detail::undir_edge_desc<EdgeRawDesc>> {
  std::size_t operator()(const bagl::container_detail::undir_edge_desc<EdgeRawDesc>& x) const {
    return bagl::container_detail::desc_get_hash(x);
  }
};

template <typename RawDesc>
struct std::hash<bagl::container_detail::indexable_desc<RawDesc>> {
  std::size_t operator()(const bagl::container_detail::indexable_desc<RawDesc>& x) const {
    return bagl::container_detail::desc_get_hash(x.value);
  }
};

#endif  // BAGL_BAGL_DETAIL_CONTAINER_GENERATORS_H_
