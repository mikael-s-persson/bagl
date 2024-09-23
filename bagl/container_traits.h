// Copyright 2004 Jeremy Siek
// Copyright 2010 Thomas Claveirole
// Copyright 2010 Ignacy Gawedzki
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CONTAINER_TRAITS_H_
#define BAGL_BAGL_CONTAINER_TRAITS_H_

// Sure would be nice to be able to forward declare these
// instead of pulling in all the headers. Too bad that
// is not legal. There ought to be a standard <stlfwd> header. -JGS

#include <algorithm>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

// The content of this file is in 'container_traits_detail' because
// it could name-clash easily.
namespace bagl::container_traits_detail {

//======================================================================
// Container Category Tags
//
//   They use virtual inheritance because there are lots of
//   inheritance diamonds.

struct container_tag {};
struct forward_container_tag : virtual public container_tag {};
struct reversible_container_tag : virtual public forward_container_tag {};
struct random_access_container_tag : virtual public reversible_container_tag {};

struct sequence_tag : virtual public forward_container_tag {};

struct associative_container_tag : virtual public forward_container_tag {};

struct sorted_associative_container_tag : virtual public associative_container_tag,
                                          virtual public reversible_container_tag {};

struct front_insertion_sequence_tag : virtual public sequence_tag {};
struct back_insertion_sequence_tag : virtual public sequence_tag {};

struct unique_associative_container_tag : virtual public associative_container_tag {};
struct multiple_associative_container_tag : virtual public associative_container_tag {};
struct simple_associative_container_tag : virtual public associative_container_tag {};
struct pair_associative_container_tag : virtual public associative_container_tag {};

//======================================================================
// Iterator Stability Tags
//
// Do mutating operations such as insert/erase/resize invalidate all
// outstanding iterators?

struct stable_tag {};
struct unstable_tag {};

//======================================================================
// Container Traits Class and container_category() function

// std::vector
struct vector_tag : virtual public random_access_container_tag, virtual public back_insertion_sequence_tag {};

template <class T, class Alloc>
vector_tag container_category(const std::vector<T, Alloc>&) {
  return {};
}

template <class T, class Alloc>
unstable_tag iterator_stability(const std::vector<T, Alloc>&) {
  return {};
}

// std::list
struct list_tag : virtual public reversible_container_tag,
                  virtual public back_insertion_sequence_tag
// this causes problems for push_dispatch...
//    virtual public front_insertion_sequence_tag
{};

template <class T, class Alloc>
list_tag container_category(const std::list<T, Alloc>&) {
  return {};
}

template <class T, class Alloc>
stable_tag iterator_stability(const std::list<T, Alloc>&) {
  return {};
}

// std::set
struct set_tag : virtual public sorted_associative_container_tag,
                 virtual public simple_associative_container_tag,
                 virtual public unique_associative_container_tag {};

template <class Key, class Cmp, class Alloc>
set_tag container_category(const std::set<Key, Cmp, Alloc>&) {
  return {};
}

template <class Key, class Cmp, class Alloc>
stable_tag iterator_stability(const std::set<Key, Cmp, Alloc>&) {
  return {};
}

// std::multiset
struct multiset_tag : virtual public sorted_associative_container_tag,
                      virtual public simple_associative_container_tag,
                      virtual public multiple_associative_container_tag {};

template <class Key, class Cmp, class Alloc>
multiset_tag container_category(const std::multiset<Key, Cmp, Alloc>&) {
  return {};
}

template <class Key, class Cmp, class Alloc>
stable_tag iterator_stability(const std::multiset<Key, Cmp, Alloc>&) {
  return {};
}

// deque

// std::map
struct map_tag : virtual public sorted_associative_container_tag,
                 virtual public pair_associative_container_tag,
                 virtual public unique_associative_container_tag {};

template <class Key, class T, class Cmp, class Alloc>
map_tag container_category(const std::map<Key, T, Cmp, Alloc>&) {
  return {};
}

template <class Key, class T, class Cmp, class Alloc>
stable_tag iterator_stability(const std::map<Key, T, Cmp, Alloc>&) {
  return {};
}

// std::multimap
struct multimap_tag : virtual public sorted_associative_container_tag,
                      virtual public pair_associative_container_tag,
                      virtual public multiple_associative_container_tag {};

template <class Key, class T, class Cmp, class Alloc>
multimap_tag container_category(const std::multimap<Key, T, Cmp, Alloc>&) {
  return {};
}

template <class Key, class T, class Cmp, class Alloc>
stable_tag iterator_stability(const std::multimap<Key, T, Cmp, Alloc>&) {
  return {};
}

// hash_set, hash_map

struct unordered_set_tag : virtual public simple_associative_container_tag,
                           virtual public unique_associative_container_tag {};

struct unordered_multiset_tag : virtual public simple_associative_container_tag,
                                virtual public multiple_associative_container_tag {};

struct unordered_map_tag : virtual public pair_associative_container_tag,
                           virtual public unique_associative_container_tag {};

struct unordered_multimap_tag : virtual public pair_associative_container_tag,
                                virtual public multiple_associative_container_tag {};

template <class Key, class Eq, class Hash, class Alloc>
unordered_set_tag container_category(const std::unordered_set<Key, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class T, class Eq, class Hash, class Alloc>
unordered_map_tag container_category(const std::unordered_map<Key, T, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class Eq, class Hash, class Alloc>
unstable_tag iterator_stability(const std::unordered_set<Key, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class T, class Eq, class Hash, class Alloc>
unstable_tag iterator_stability(const std::unordered_map<Key, T, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class Eq, class Hash, class Alloc>
unordered_multiset_tag container_category(const std::unordered_multiset<Key, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class T, class Eq, class Hash, class Alloc>
unordered_multimap_tag container_category(const std::unordered_multimap<Key, T, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class Eq, class Hash, class Alloc>
unstable_tag iterator_stability(const std::unordered_multiset<Key, Eq, Hash, Alloc>&) {
  return {};
}

template <class Key, class T, class Eq, class Hash, class Alloc>
unstable_tag iterator_stability(const std::unordered_multimap<Key, T, Eq, Hash, Alloc>&) {
  return {};
}

// don't use this unless there is partial specialization
template <class Container>
struct container_traits {
  using category = decltype(container_category(std::declval<Container>()));
  using stability = decltype(iterator_stability(std::declval<Container>()));
};

//===========================================================================
// Generalized Container Functions

// Erase
template <class Sequence, class T>
void erase_dispatch(Sequence& c, const T& x, sequence_tag) {
  c.erase(std::remove(c.begin(), c.end(), x), c.end());
}

template <class AssociativeContainer, class T>
void erase_dispatch(AssociativeContainer& c, const T& x, associative_container_tag) {
  c.erase(x);
}
template <class Container, class T>
void erase(Container& c, const T& x) {
  erase_dispatch(c, x, container_category(c));
}

// Erase If
template <class Sequence, class Predicate, class IteratorStability>
void erase_if_dispatch(Sequence& c, Predicate p, sequence_tag, IteratorStability) {
  if (!c.empty()) {
    c.erase(std::remove_if(c.begin(), c.end(), p), c.end());
  }
}
template <class AssociativeContainer, class Predicate>
void erase_if_dispatch(AssociativeContainer& c, Predicate p, associative_container_tag, stable_tag) {
  for (auto i = c.begin(), next = c.begin(); next != c.end(); i = next) {
    ++next;
    if (p(*i)) {
      c.erase(i);
    }
  }
}
template <class AssociativeContainer, class Predicate>
void erase_if_dispatch(AssociativeContainer& c, Predicate p, associative_container_tag, unstable_tag) {
  // This method is really slow, so hopefully we won't have any
  // associative containers with unstable iterators!
  // Is there a better way to do this?
  auto n = c.size();
  while (n--) {
    for (auto i = c.begin(); i != c.end(); ++i) {
      if (p(*i)) {
        c.erase(i);
        break;
      }
    }
  }
}
template <class Container, class Predicate>
void erase_if(Container& c, Predicate p) {
  erase_if_dispatch(c, p, container_category(c), iterator_stability(c));
}

// Push
template <class Container, class T>
auto push_dispatch(Container& c, T&& v, back_insertion_sequence_tag) {
  c.push_back(std::forward<T>(v));
  return std::make_pair(std::prev(c.end()), true);
}

template <class Container, class T>
auto push_dispatch(Container& c, T&& v, front_insertion_sequence_tag) {
  c.push_front(std::forward<T>(v));
  return std::make_pair(c.begin(), true);
}

template <class AssociativeContainer, class T>
auto push_dispatch(AssociativeContainer& c, T&& v, unique_associative_container_tag) {
  return c.insert(std::forward<T>(v));
}

template <class AssociativeContainer, class T>
auto push_dispatch(AssociativeContainer& c, T&& v, multiple_associative_container_tag) {
  return std::make_pair(c.insert(std::forward<T>(v)), true);
}

template <class Container, class T>
auto push(Container& c, T&& v) {
  return push_dispatch(c, std::forward<T>(v), container_category(c));
}

// Find
template <class Container, class Value>
auto find_dispatch(Container& c, const Value& value, container_tag) {
  return std::find(c.begin(), c.end(), value);
}

template <class AssociativeContainer, class Value>
auto find_dispatch(AssociativeContainer& c, const Value& value, associative_container_tag) {
  return c.find(value);
}

template <class Container, class Value>
auto find(Container& c, const Value& value) {
  return find_dispatch(c, value, container_category(c));
}

// Find (const versions)
template <class Container, class Value>
auto find_dispatch(const Container& c, const Value& value, container_tag) {
  return std::find(c.begin(), c.end(), value);
}

template <class AssociativeContainer, class Value>
auto find_dispatch(const AssociativeContainer& c, const Value& value, associative_container_tag) {
  return c.find(value);
}

template <class Container, class Value>
auto find(const Container& c, const Value& value) {
  return find_dispatch(c, value, container_category(c));
}

// Equal range

template <class AssociativeContainer, class Value>
auto equal_range_dispatch(AssociativeContainer& c, const Value& value, associative_container_tag) {
  return c.equal_range(value);
}

template <class Container, class Value>
auto equal_range(Container& c, const Value& value) {
  return equal_range_dispatch(c, value, container_category(c));
}

}  // namespace bagl::container_traits_detail

#endif  // BAGL_BAGL_CONTAINER_TRAITS_H_
