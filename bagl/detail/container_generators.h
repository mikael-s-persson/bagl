// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_DETAIL_CONTAINER_GENERATORS_H_
#define BAGL_BAGL_DETAIL_CONTAINER_GENERATORS_H_

#include <limits>
#include <list>
#include <set>
#include <ranges>
#include <type_traits>
#include <unordered_set>
#include <variant>
#include <vector>

#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"

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

template <typename Container>
std::size_t get_size(const Container& cont) {
  return cont.size();
}

template <typename Container>
std::size_t get_size(Container* cont) {
  return cont->size();
}

template <typename Container>
std::size_t get_capacity(const Container& cont) {
  return cont.capacity();
}

template <typename Container>
std::size_t get_capacity(Container* cont) {
  return cont->capacity();
}

template <typename Container>
void clear_all(Container& cont) {
  cont.clear();
}

template <typename Container>
void clear_all(Container* cont) {
  cont->clear();
}

}  // namespace container_detail

struct vec_s {};
struct pool_s {};
struct list_s {};
struct set_s {};
struct multiset_s {};
struct map_s {};
struct multimap_s {};
struct unordered_set_s {};
struct unordered_multiset_s {};
struct unordered_map_s {};
struct unordered_multimap_s {};

template <class Selector, class ValueType>
struct container_gen {};

template <class ValueType>
struct container_gen<vec_s, ValueType> {
  using type = std::vector<ValueType>;
};

template <class ValueType>
struct container_gen<pool_s, ValueType> {
  using type = container_detail::pooled_vector<ValueType>;
};

template <class ValueType>
struct container_gen<list_s, ValueType> {
  using type = std::list<ValueType>;
};

template <class ValueType>
struct container_gen<set_s, ValueType> {
  using type = std::set<ValueType>;
};

template <class ValueType>
struct container_gen<multiset_s, ValueType> {
  using type = std::multiset<ValueType>;
};

template <class ValueType>
struct container_gen<map_s, ValueType> {
  using type = std::set<ValueType>;
};

template <class ValueType>
struct container_gen<multimap_s, ValueType> {
  using type = std::multiset<ValueType>;
};

template <class ValueType>
struct container_gen<unordered_set_s, ValueType> {
  using type = std::unordered_set<ValueType>;
};

template <class ValueType>
struct container_gen<unordered_multiset_s, ValueType> {
  using type = std::unordered_multiset<ValueType>;
};

template <class ValueType>
struct container_gen<unordered_map_s, ValueType> {
  using type = std::unordered_set<ValueType>;
};

template <class ValueType>
struct container_gen<unordered_multimap_s, ValueType> {
  using type = std::unordered_multiset<ValueType>;
};

namespace container_detail {

template <class Selector>
struct parallel_edge_traits {
  using type = allow_parallel_edge_tag;
};

template <>
struct parallel_edge_traits<set_s> {
  using type = disallow_parallel_edge_tag;
};

template <>
struct parallel_edge_traits<unordered_set_s> {
  using type = disallow_parallel_edge_tag;
};

// map_s is obsolete, replaced with set_s
template <>
struct parallel_edge_traits<map_s> {
  using type = disallow_parallel_edge_tag;
};

template <>
struct parallel_edge_traits<unordered_map_s> {
  using type = disallow_parallel_edge_tag;
};

template <typename Selector>
struct is_random_access : std::false_type {};

template <>
struct is_random_access<vec_s> : std::true_type {};

template <typename DescType>
struct null_desc {
  static DescType value() { return DescType{}; }
};

template <typename DescType>
struct null_desc<DescType*> {
  static DescType* value() { return nullptr; }
};

template <>
struct null_desc<std::size_t> {
  static std::size_t value() { return (std::numeric_limits<std::size_t>::max)(); }
};

inline bool desc_less_than(std::size_t lhs, std::size_t rhs) { return (lhs < rhs); }

template <typename Iter>
bool desc_less_than(Iter lhs, Iter rhs) {
  return (&(*lhs) < &(*rhs));
}

inline std::size_t desc_get_hash(std::size_t d) {
  std::hash<std::size_t> hasher;
  return hasher(d);
}

template <typename Iter>
std::size_t desc_get_hash(Iter it) {
  using ValueType = typename Iter::value_type;
  std::hash<ValueType*> hasher;
  return hasher(&(*it));
}

struct desc_hasher {
  std::size_t operator()(std::size_t d) const {
    std::hash<std::size_t> hasher;
    return hasher(d);
  }
  template <typename Iter>
  std::size_t operator()(Iter it) const {
    using ValueType = typename Iter::value_type;
    std::hash<ValueType*> hasher;
    return hasher(&(*it));
  }
};

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

template <typename Container>
auto desc_to_iterator(Container& c, std::size_t d) {
  return c.begin() + d;
}

template <typename Container>
auto desc_to_iterator(const Container& c, std::size_t d) {
  return c.begin() + d;
}

template <typename Container, typename Iter>
Iter desc_to_iterator(const Container& /*unused*/, Iter it) {
  return it;
}

template <typename ValueType>
std::size_t iterator_to_desc(const std::vector<ValueType>& c, typename std::vector<ValueType>::iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t iterator_to_desc(const std::vector<ValueType>& c, typename std::vector<ValueType>::const_iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t iterator_to_desc(const pooled_vector<ValueType>& c, typename pooled_vector<ValueType>::iterator it) {
  return it - c.begin();
}

template <typename ValueType>
std::size_t iterator_to_desc(const pooled_vector<ValueType>& c, typename pooled_vector<ValueType>::const_iterator it) {
  return it - c.begin();
}

template <typename Container, typename Iter>
Iter iterator_to_desc(const Container& /*unused*/, Iter it) {
  return it;
}

template <typename Container>
auto get_begin_iter(Container& c) {
  return c.begin();
}
template <typename Container>
auto get_begin_iter(Container* c) {
  return get_begin_iter(*c);
}

template <typename ValueType>
std::size_t get_begin_desc(const std::vector<ValueType>& c) {
  return 0;
}
template <typename ValueType>
std::size_t get_begin_desc(const pooled_vector<ValueType>& c) {
  return 0;
}
template <typename Container>
auto get_begin_desc(Container& c) {
  return c.begin();
}
template <typename Container>
auto get_begin_desc(Container* c) {
  return get_begin_desc(*c);
}

template <typename ValueType>
std::size_t get_end_desc(const std::vector<ValueType>& c) {
  return c.size();
}
template <typename ValueType>
std::size_t get_end_desc(const pooled_vector<ValueType>& c) {
  return c.m_data.size();
}
template <typename Container>
auto get_end_desc(Container& c) {
  return c.end();
}
template <typename Container>
auto get_end_desc(Container* c) {
  return get_end_desc(*c);
}
template <typename Container>
auto get_end_iter(Container& c) {
  return c.end();
}
template <typename Container>
auto get_end_iter(Container* c) {
  return get_end_iter(*c);
}

template <typename Container>
auto get_range(Container& c) {
  return std::views::all(c);
}
template <typename Container>
auto get_range(Container* c) {
  return get_range(*c);
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

template <typename Container, typename Desc>
auto& desc_to_value(Container* p_c, Desc d) {
  return desc_to_value(*p_c, d);
}

template <typename Container>
struct select_descriptor {
  using type = typename Container::iterator;
};

template <typename Container>
struct select_descriptor<Container*> {
  using type = typename Container::iterator;
};

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

  ignore_output_iter& operator++() { return *this; }
  ignore_output_iter& operator++(int) { return *this; }
  ignore_output_iter& operator--() { return *this; }
  ignore_output_iter& operator--(int) { return *this; }

  value_type operator*() const { return {}; }
};

}  // namespace container_detail

}  // namespace bagl

#endif  // BAGL_BAGL_DETAIL_CONTAINER_GENERATORS_H_
