// Copyright 2010 Daniel Trebbien.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BUFFER_CONCEPTS_H_
#define BAGL_BAGL_BUFFER_CONCEPTS_H_

#include <queue>
#include <type_traits>
#include <vector>

#include "bagl/has_trait_member.h"
#include "bagl/property_map.h"

namespace bagl {
namespace concepts {

namespace buffer_concepts_detail {

struct invalid_type {};
BAGL_GRAPH_HAS_TRAIT_MEMBER(value_type, invalid_type)
BAGL_GRAPH_HAS_TRAIT_MEMBER(key_map, invalid_type)

}  // namespace buffer_concepts_detail

template <typename B>
concept Buffer = requires(B& buf, const buffer_concepts_detail::get_value_type_or_not<B>& value) {
  buf.push(value);
  buf.pop();
  { buf.top() } -> std::convertible_to<buffer_concepts_detail::get_value_type_or_not<B>&>;
}
&&requires(const B& buf) {
  { buf.top() } -> std::convertible_to<const buffer_concepts_detail::get_value_type_or_not<B>&>;
  { buf.empty() } -> std::convertible_to<bool>;
};

template <typename Q>
concept UpdatableQueue = Buffer<Q> && requires(Q& q, const buffer_concepts_detail::get_value_type_or_not<Q>& value) {
  q.update(value);
} && requires(const Q& q, const buffer_concepts_detail::get_value_type_or_not<Q>& value) {
  { q.contains(value) } -> std::convertible_to<bool>;
};

template <typename Q>
concept KeyedUpdatableQueue = UpdatableQueue<Q> &&
    ReadWritePropertyMap<buffer_concepts_detail::get_key_map_or_not<Q>,
                         buffer_concepts_detail::get_value_type_or_not<Q>> &&
    requires(const Q& q, const buffer_concepts_detail::get_value_type_or_not<Q>& value) {
  { q.keys() } -> std::convertible_to<buffer_concepts_detail::get_key_map_or_not<Q>>;
};

}  // namespace concepts

template <typename T>
class buffer_queue : public std::queue<T> {
 public:
  using std::queue<T>::queue;
  decltype(auto) top() const { return this->front(); }
  decltype(auto) top() { return this->front(); }
};
template <typename T>
class buffer_vector : public std::vector<T> {
 public:
  using std::vector<T>::vector;
  decltype(auto) top() const { return this->front(); }
  decltype(auto) top() { return this->front(); }
  void push(const T& t) { this->push_back(t); }
  void pop() { this->erase(this->begin()); }
};

}  // namespace bagl

#endif  // BAGL_BAGL_BUFFER_CONCEPTS_H_
