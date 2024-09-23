// Copyright 2024 Mikael Persson

#ifndef BAGL_BAGL_HAS_TRAIT_MEMBER_H_
#define BAGL_BAGL_HAS_TRAIT_MEMBER_H_

#include <type_traits>

#define BAGL_GRAPH_HAS_TRAIT_MEMBER(name, default_type)                          \
  template <typename T>                                                          \
  using has_##name##_detect = typename T::name;                                  \
  template <typename T, typename = std::void_t<>>                                \
  struct has_##name : std::false_type {};                                        \
  template <typename T>                                                          \
  struct has_##name<T, std::void_t<has_##name##_detect<T>>> : std::true_type {}; \
  template <typename T>                                                          \
  constexpr bool has_##name##_v = has_##name<T>::value;                   \
  template <typename T>                                                          \
  using get_##name##_or_not = std::conditional_t<has_##name##_v<T>, typename T::name, default_type>;

#define BAGL_GRAPH_HAS_NONTYPE_MEMBER(name)                          \
  template <typename T>                                                          \
  using has_##name##_detect = decltype(&T::name);                                \
  template <typename T, typename = std::void_t<>>                                \
  struct has_##name : std::false_type {};                                        \
  template <typename T>                                                          \
  struct has_##name<T, std::void_t<has_##name##_detect<T>>> : std::true_type {}; \
  template <typename T>                                                          \
  constexpr bool has_##name##_v = has_##name<T>::value;

#define BAGL_GRAPH_HAS_MEMBER_FUNCTION(name)                          \
  template <typename T>                                                          \
  using has_##name##_detect = decltype(&T::name);                                \
  template <typename T, typename = std::void_t<>>                                \
  struct has_##name : std::false_type {};                                        \
  template <typename T>                                                          \
  struct has_##name<T, std::void_t<has_##name##_detect<T>>> : std::true_type {}; \
  template <typename T>                                                          \
  constexpr bool has_##name##_v = has_##name<T>::value; \
  template <typename T, typename... Args> \
  auto invoke_##name(T& t, Args&&... args) { \
    if constexpr ( has_##name##_v<std::decay_t<T>> ) { \
      return t.name(std::forward<Args>(args)...); \
    } \
  } \
  template <typename T, typename... Args> \
  void invoke_##name##_on_all(T& t, Args&&... args) { \
    std::apply([&args...](auto&&... ts) { \
      (invoke_##name(ts, std::forward<Args>(args)...),...); \
    }, t); \
  }

#endif // BAGL_BAGL_HAS_TRAIT_MEMBER_H_
