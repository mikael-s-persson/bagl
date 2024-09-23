// Copyright 2010 Daniel Trebbien.
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BUFFER_CONCEPTS_H_
#define BAGL_BAGL_BUFFER_CONCEPTS_H_

#include <type_traits>

#include "bagl/has_trait_member.h"
#include "bagl/property_map.h"

namespace bagl {

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
  { buf.size() } -> std::integral;
  { buf.empty() } -> std::convertible_to<bool>;
};

template <typename Q>
concept UpdatableQueue = Buffer<Q> && requires(Q& q, const buffer_concepts_detail::get_value_type_or_not<Q>& value) {
  q.update(value);
} && requires(const Q& q, const buffer_concepts_detail::get_value_type_or_not<Q>& value) {
  { q.contains() } -> std::convertible_to<bool>;
};

template <typename Q>
concept KeyedUpdatableQueue = UpdatableQueue<Q> &&
    ReadWritePropertyMap<buffer_concepts_detail::get_key_map_or_not<Q>,
                         buffer_concepts_detail::get_value_type_or_not<Q>> &&
    requires(const Q& q, const buffer_concepts_detail::get_value_type_or_not<Q>& value) {
  { q.keys() } -> std::convertible_to<buffer_concepts_detail::get_key_map_or_not<Q>>;
};

}  // namespace bagl

#endif  // BAGL_BAGL_BUFFER_CONCEPTS_H_
