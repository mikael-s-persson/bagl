// Copyright 2004 Jeremy Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PROPERTY_H_
#define BAGL_BAGL_PROPERTY_H_

#include <string_view>
#include <type_traits>
#include <utility>

#include "bagl/has_trait_member.h"

namespace bagl {

struct no_property {};

template <typename Tag, typename T, typename Base = no_property>
struct property {
  using next_type = Base;
  using tag_type = Tag;
  using value_type = T;

  using self = property<Tag, T, Base>;

  property() = default;

  template <typename U, typename... BaseUs>
  requires std::constructible_from<T, U&&>
  explicit property(U&& v, BaseUs&&... bs) : m_value(std::forward<U>(v)), m_base(std::forward<BaseUs>(bs)...) {}

  template <typename U, typename BaseU>
  // NOLINTNEXTLINE(google-explicit-constructor) Generic copy-constructor.
  property(const property<Tag, U, BaseU>& rhs) : m_value(rhs.m_value), m_base(rhs.m_base) {}

  template <typename U, typename BaseU>
  // NOLINTNEXTLINE(google-explicit-constructor) Generic move-constructor.
  property(property<Tag, U, BaseU>&& rhs) noexcept : m_value(std::move(rhs.m_value)), m_base(std::move(rhs.m_base)) {}

  template <typename U, typename BaseU>
  self& operator=(const property<Tag, U, BaseU>& rhs) {
    if (this == &rhs) {
      return *this;
    }
    m_value = rhs.m_value;
    m_base = rhs.m_base;
    return *this;
  }
  template <typename U, typename BaseU>
  self& operator=(property<Tag, U, BaseU>&& rhs) {
    if (this == &rhs) {
      return *this;
    }
    m_value = std::move(rhs.m_value);
    m_base = std::move(rhs.m_base);
    return *this;
  }

  T m_value;
  [[no_unique_address]] Base m_base;
};

template <typename Tag, typename T, typename Base>
Base& remove_top_property(property<Tag, T, Base>& p) {
  return p.m_base;
}
template <typename Tag, typename T, typename Base>
const Base& remove_top_property(const property<Tag, T, Base>& p) {
  return p.m_base;
}
struct top_property_remover {
  template <typename Tag, typename T, typename Base>
  Base& operator()(property<Tag, T, Base>& p) const {
    return p.m_base;
  }
  template <typename Tag, typename T, typename Base>
  const Base& operator()(const property<Tag, T, Base>& p) const {
    return p.m_base;
  }
};

template <typename Tag, typename T, typename Base>
auto add_top_property(Tag /*unused*/, T&& top_p, Base&& base_p) {
  return property<Tag, std::decay_t<T>, std::decay_t<Base>>(std::forward<T>(top_p), std::forward<Base>(base_p));
}
struct top_property_adder {
  template <typename Tag, typename T, typename Base>
  auto operator()(Tag /*unused*/, T&& top_p, Base&& base_p) const {
    return add_top_property(Tag{}, std::forward<T>(top_p), std::forward<Base>(base_p));
  }
};
template <typename Tag, typename T>
struct top_property_value_adder {
  T value;
  template <typename Base>
  auto operator()(Base&& base_p) const {
    return add_top_property(Tag{}, value, std::forward<Base>(base_p));
  }
};

// Kinds of properties
struct graph_property_tag {};
struct vertex_property_tag {};
struct edge_property_tag {};

// Kinds of properties
namespace property_detail {
BAGL_GRAPH_HAS_TRAIT_MEMBER(kind, void)
}  // namespace property_detail

// Having a default is to make this trait work for any type, not just valid
// properties, to work around VC++ <= 10 bugs related to SFINAE in
// compressed_sparse_row_graph's get functions and similar
template <typename PropertyTag>
struct property_kind {
    using type = property_detail::get_kind_or_not<PropertyTag>;
};
template <typename PropertyTag>
using property_kind_t = typename property_kind<PropertyTag>::type;

template <typename PropertyTag>
constexpr bool is_vertex_property_kind_v = std::is_same_v<property_kind_t<PropertyTag>, vertex_property_tag>;
template <typename PropertyTag>
constexpr bool is_edge_property_kind_v = std::is_same_v<property_kind_t<PropertyTag>, edge_property_tag>;
template <typename PropertyTag>
constexpr bool is_graph_property_kind_v = std::is_same_v<property_kind_t<PropertyTag>, graph_property_tag>;

// Some standard properties defined independently:
struct vertex_all_t {
  using kind = vertex_property_tag;
  static constexpr std::string_view name = "vertex_all";
};
constexpr vertex_all_t vertex_all = {};
struct edge_all_t {
  using kind = edge_property_tag;
  static constexpr std::string_view name = "edge_all";
};
constexpr edge_all_t edge_all = {};
struct graph_all_t {
  using kind = graph_property_tag;
  static constexpr std::string_view name = "graph_all";
};
constexpr graph_all_t graph_all = {};

// These tags are used for property bundles
struct vertex_bundle_t {
  using kind = vertex_property_tag;
  static constexpr std::string_view name = "vertex_bundle";
};
constexpr vertex_bundle_t vertex_bundle = {};
struct edge_bundle_t {
  using kind = edge_property_tag;
  static constexpr std::string_view name = "edge_bundle";
};
constexpr edge_bundle_t edge_bundle = {};
struct graph_bundle_t {
  using kind = graph_property_tag;
  static constexpr std::string_view name = "graph_bundle";
};
constexpr graph_bundle_t graph_bundle = {};

// Code to look up one property in a property list:
template <typename PList, typename PropName, typename Enable = void>
struct lookup_one_property_internal {
  static constexpr bool found = false;
  using type = no_property;
};

// Special-case properties (vertex_all, edge_all, graph_all)
// This "lookup" unconditionally finds the top-level property as matching the tag.
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) Immediately undef'd
#define BAGL_PROPERTY_ALL_PROP(tag)                                                        \
  template <typename T>                                                                    \
  struct lookup_one_property_internal<T, tag> {                                            \
    static constexpr bool found = true;                                                    \
    using type = T;                                                                        \
    static T& lookup(T& x, tag) { return x; }                                              \
    static const T& lookup(const T& x, tag) { return x; }                                  \
  };                                                                                       \
  template <typename Tag, typename T, typename Base>                                       \
  struct lookup_one_property_internal<property<Tag, T, Base>, tag> { /* Avoid ambiguity */ \
    static constexpr bool found = true;                                                    \
    using type = property<Tag, T, Base>;                                                   \
    static type& lookup(type& x, tag) { return x; }                                        \
    static const type& lookup(const type& x, tag) { return x; }                            \
  };

BAGL_PROPERTY_ALL_PROP(vertex_all_t)
BAGL_PROPERTY_ALL_PROP(edge_all_t)
BAGL_PROPERTY_ALL_PROP(graph_all_t)
#undef BAGL_PROPERTY_ALL_PROP

// *_bundled; these need to be macros rather than inheritance to resolve
// ambiguities
// This lookup finds any non-property element as being the bundled property.
//  - lookup<Foo, *_bundle_t> finds 'Foo'
//  - lookup<property<Tag1, T1, property<Tag2, T2, Foo>>, *_bundle_t> finds 'Foo'
// NOLINTNEXTLINE(cppcoreguidelines-macro-usage) Immediately undef'd
#define BAGL_PROPERTY_DO_ONE_BUNDLE_TYPE(kind)                                                \
  template <typename T>                                                                       \
  struct lookup_one_property_internal<T, kind##_bundle_t> {                                   \
    static constexpr bool found = true;                                                       \
    using type = T;                                                                           \
    static T& lookup(T& x, kind##_bundle_t) { return x; }                                     \
    static const T& lookup(const T& x, kind##_bundle_t) { return x; }                         \
  };                                                                                          \
                                                                                              \
  template <typename Tag, typename T, typename Base>                                          \
  struct lookup_one_property_internal<property<Tag, T, Base>, kind##_bundle_t>                \
      : lookup_one_property_internal<Base, kind##_bundle_t> {                                 \
   private:                                                                                   \
    using base_type = lookup_one_property_internal<Base, kind##_bundle_t>;                    \
                                                                                              \
   public:                                                                                    \
    template <typename BundleTag>                                                             \
    static std::enable_if_t<(base_type::found && std::is_same_v<BundleTag, kind##_bundle_t>), \
                            typename base_type::type&>                                        \
    lookup(property<Tag, T, Base>& p, BundleTag) {                                            \
      return base_type::lookup(p.m_base, kind##_bundle_t());                                  \
    }                                                                                         \
    template <typename BundleTag>                                                             \
    static std::enable_if_t<(base_type::found && std::is_same_v<BundleTag, kind##_bundle_t>), \
                            const typename base_type::type&>                                  \
    lookup(const property<Tag, T, Base>& p, BundleTag) {                                      \
      return base_type::lookup(p.m_base, kind##_bundle_t());                                  \
    }                                                                                         \
  };

BAGL_PROPERTY_DO_ONE_BUNDLE_TYPE(vertex)
BAGL_PROPERTY_DO_ONE_BUNDLE_TYPE(edge)
BAGL_PROPERTY_DO_ONE_BUNDLE_TYPE(graph)
#undef BAGL_PROPERTY_DO_ONE_BUNDLE_TYPE

// Normal old-style properties; second case also handles chaining of bundled
// property accesses
template <typename Tag, typename T, typename Base>
struct lookup_one_property_internal<property<Tag, T, Base>, Tag> {
  static constexpr bool found = true;
  using prop = property<Tag, T, Base>;
  using type = T;
  static T& lookup(prop& p, const Tag& /*unused*/) {
    return p.m_value;
  }
  static const T& lookup(const prop& p, const Tag& /*unused*/) {
    return p.m_value;
  }
};

template <typename Tag, typename T, typename Base, typename PropName>
struct lookup_one_property_internal<property<Tag, T, Base>, PropName> : lookup_one_property_internal<Base, PropName> {
 private:
  using base_type = lookup_one_property_internal<Base, PropName>;

 public:
  static typename base_type::type& lookup(property<Tag, T, Base>& prop, const PropName& tag) {
    return base_type::lookup(prop.m_base, tag);
  }
  static const typename base_type::type& lookup(const property<Tag, T, Base>& prop, const PropName& tag) {
    return base_type::lookup(prop.m_base, tag);
  }
};

// Pointer-to-member access to bundled properties
template <typename T, typename TMaybeBase, typename R>
struct lookup_one_property_internal<T, R TMaybeBase::*, std::enable_if_t<std::is_base_of_v<TMaybeBase, T>>> {
  static constexpr bool found = true;
  using type = R;
  static R& lookup(T& x, R TMaybeBase::*ptr) { return x.*ptr; }
  static const R& lookup(const T& x, R TMaybeBase::*ptr) { return x.*ptr; }
};

// Version of above handling const property lists properly
template <typename T, typename Tag>
struct lookup_one_property : lookup_one_property_internal<T, Tag> {};

template <typename T, typename Tag>
struct lookup_one_property<const T, Tag> {
  static constexpr bool found = lookup_one_property_internal<T, Tag>::found;
  using type = const typename lookup_one_property_internal<T, Tag>::type;
  static auto& lookup(const T& p, Tag tag) {
    return lookup_one_property_internal<T, Tag>::lookup(p, tag);
  }
};

template <typename T, typename Tag>
using lookup_one_property_t = typename lookup_one_property<T, Tag>::type;

// The BGL properties specialize property_kind and
// property_num, and use enum's for the Property type (see
// graph/properties.hpp), but the user may want to use a class
// instead with a nested kind type and num.  Also, we may want to
// switch BGL back to using class types for properties at some point.

template <typename P>
constexpr bool has_property_v = !std::is_same_v<P, no_property>;

template <class PropertyList, class Tag>
struct property_value : lookup_one_property<PropertyList, Tag> {};

template <class PropertyList, class Tag>
auto& get_property_value(PropertyList& p, Tag tag) {
  return lookup_one_property<PropertyList, Tag>::lookup(p, tag);
}

template <class PropertyList, class Tag>
const auto& get_property_value(const PropertyList& p, Tag tag) {
  return lookup_one_property<PropertyList, Tag>::lookup(p, tag);
}

namespace property_detail {

/** This trait returns true if T is no_property. */
template <typename T>
constexpr bool is_no_property_v = std::is_same_v<T, no_property>;

template <typename PList, typename Tag>
class lookup_one_property_f {
  Tag tag_;

 public:
  explicit lookup_one_property_f(Tag tag) : tag_(tag) {}

  decltype(auto) operator()(PList& pl) const {
    return lookup_one_property<PList, Tag>::lookup(pl, tag_);
  }
};

}  // namespace property_detail

}  // namespace bagl

#endif  // BAGL_BAGL_PROPERTY_H_
