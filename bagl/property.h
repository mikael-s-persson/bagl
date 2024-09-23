// Copyright 2004 Jeremy Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PROPERTY_H_
#define BAGL_BAGL_PROPERTY_H_

#include <type_traits>

#include "bagl/has_trait_member.h"

namespace bagl {

struct no_property {};

template <typename Tag, typename T, typename Base = no_property>
struct property {
  typedef Base next_type;
  typedef Tag tag_type;
  typedef T value_type;
  property(const T& v = T()) : m_value(v) {}
  property(const T& v, const Base& b) : m_value(v), m_base(b) {}
  // copy constructor and assignment operator will be generated by compiler

  T m_value;
  BOOST_ATTRIBUTE_NO_UNIQUE_ADDRESS Base m_base;
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


// Some standard properties defined independently:
struct vertex_all_t { using kind = vertex_property_tag; };
constexpr vertex_all_t vertex_all = {};
struct edge_all_t { using kind = edge_property_tag; };
constexpr edge_all_t edge_all = {};
struct graph_all_t { using kind = graph_property_tag; };
constexpr graph_all_t graph_all = {};

// These tags are used for property bundles
struct vertex_bundle_t { using kind = vertex_property_tag; };
constexpr vertex_bundle_t vertex_bundle = {};
struct edge_bundle_t { using kind = edge_property_tag; };
constexpr edge_bundle_t edge_bundle = {};
struct graph_bundle_t { using kind = graph_property_tag; };
constexpr graph_bundle_t graph_bundle = {};

// Code to look up one property in a property list:
template <typename PList, typename PropName, typename Enable = void>
struct lookup_one_property_internal {
  static constexpr bool found = false;
  using type = void;
};

// Special-case properties (vertex_all, edge_all, graph_all)
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
                            std::add_reference_t<typename base_type::type>>                   \
    lookup(property<Tag, T, Base>& p, BundleTag) {                                            \
      return base_type::lookup(p.m_base, kind##_bundle_t());                                  \
    }                                                                                         \
    template <typename BundleTag>                                                             \
    static std::enable_if_t<(base_type::found && std::is_same_v<BundleTag, kind##_bundle_t>), \
                            std::add_reference_t<const typename base_type::type>>             \
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
  template <typename U>
  static std::enable_if_t<std::is_same_v<prop, U>, T&> lookup(U& prop, const Tag&) {
    return prop.m_value;
  }
  template <typename U>
  static std::enable_if_t<std::is_same_v<prop, U>, const T&> lookup(const U& prop, const Tag&) {
    return prop.m_value;
  }
};

template <typename Tag, typename T, typename Base, typename PropName>
struct lookup_one_property_internal<property<Tag, T, Base>, PropName> : lookup_one_property_internal<Base, PropName> {
 private:
  using base_type = lookup_one_property_internal<Base, PropName>;

 public:
  template <typename PL>
  static std::enable_if_t<std::is_same_v<PL, property<Tag, T, Base>>, std::add_reference_t<typename base_type::type>>
  lookup(PL& prop, const PropName& tag) {
    return base_type::lookup(prop.m_base, tag);
  }
  template <typename PL>
  static std::enable_if_t<std::is_same_v<PL, property<Tag, T, Base>>,
                          std::add_reference_t<const typename base_type::type>>
  lookup(const PL& prop, const PropName& tag) {
    return base_type::lookup(prop.m_base, tag);
  }
};

// Pointer-to-member access to bundled properties
#ifndef BOOST_GRAPH_NO_BUNDLED_PROPERTIES
template <typename T, typename TMaybeBase, typename R>
struct lookup_one_property_internal<T, R TMaybeBase::*, std::enable_if_t<std::is_base_of_v<TMaybeBase, T>>> {
  static constexpr bool found = true;
  using type = R;
  static R& lookup(T& x, R TMaybeBase::*ptr) { return x.*ptr; }
  static const R& lookup(const T& x, R TMaybeBase::*ptr) { return x.*ptr; }
};
#endif

// Version of above handling const property lists properly
template <typename T, typename Tag>
struct lookup_one_property : lookup_one_property_internal<T, Tag> {};

template <typename T, typename Tag>
struct lookup_one_property<const T, Tag> {
  static constexpr bool found = lookup_one_property_internal<T, Tag>::found;
  using type = const typename lookup_one_property_internal<T, Tag>::type;
  template <typename U>
  static std::enable_if_t<std::is_same_v<T, U>,
                          std::add_reference_t<const typename lookup_one_property_internal<T, Tag>::type>>
  lookup(const U& p, Tag tag) {
    return lookup_one_property_internal<T, Tag>::lookup(p, tag);
  }
};

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
inline typename lookup_one_property<PropertyList, Tag>::type& get_property_value(PropertyList& p, Tag tag) {
  return lookup_one_property<PropertyList, Tag>::lookup(p, tag);
}

template <class PropertyList, class Tag>
inline const typename lookup_one_property<PropertyList, Tag>::type& get_property_value(const PropertyList& p, Tag tag) {
  return lookup_one_property<PropertyList, Tag>::lookup(p, tag);
}

namespace property_detail {

/** This trait returns true if T is no_property. */
template <typename T>
constexpr bool is_no_property_v std::is_same_v<T, no_property>;

template <typename PList, typename Tag>
class lookup_one_property_f;

template <typename PList, typename Tag, typename F>
struct lookup_one_property_f_result;

template <typename PList, typename Tag>
struct lookup_one_property_f_result<PList, Tag, const lookup_one_property_f<PList, Tag>(PList)> {
  using type = typename lookup_one_property<PList, Tag>::type;
};

template <typename PList, typename Tag>
struct lookup_one_property_f_result<PList, Tag, const lookup_one_property_f<PList, Tag>(PList&)> {
  using type = typename lookup_one_property<PList, Tag>::type&;
};

template <typename PList, typename Tag>
struct lookup_one_property_f_result<PList, Tag, const lookup_one_property_f<PList, Tag>(const PList&)> {
  using type = const typename lookup_one_property<PList, Tag>::type&;
};

template <typename PList, typename Tag>
class lookup_one_property_f {
  Tag tag;

 public:
  lookup_one_property_f(Tag tag) : tag(tag) {}
  template <typename F>
  struct result : lookup_one_property_f_result<PList, Tag, F> {};

  typename lookup_one_property_f_result<PList, Tag, const lookup_one_property_f(PList&)>::type operator()(
      PList& pl) const {
    return lookup_one_property<PList, Tag>::lookup(pl, tag);
  }
};

}  // namespace property_detail

#if 0
namespace property_detail
{
    // Stuff for directed_graph and undirected_graph to skip over their first
    // vertex_index and edge_index properties when providing vertex_all and
    // edge_all; make sure you know the exact structure of your properties if
    // you use there.
    struct remove_first_property
    {
        template < typename F > struct result
        {
            typedef typename boost::function_traits< F >::arg1_type a1;
            typedef typename boost::remove_reference< a1 >::type non_ref;
            typedef typename non_ref::next_type nx;
            typedef typename boost::mpl::if_< boost::is_const< non_ref >,
                boost::add_const< nx >, nx >::type with_const;
            typedef typename boost::add_reference< with_const >::type type;
        };
        template < typename Prop >
        typename Prop::next_type& operator()(Prop& p) const
        {
            return p.m_base;
        }
        template < typename Prop >
        const typename Prop::next_type& operator()(const Prop& p) const
        {
            return p.m_base;
        }
    };
}
#endif

}  // namespace bagl

#endif  // BAGL_BAGL_PROPERTY_H_
