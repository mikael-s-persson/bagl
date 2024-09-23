// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DETAIL_EDGE_H_
#define BAGL_BAGL_DETAIL_EDGE_H_

#include <iosfwd>
#include <utility>

namespace bagl::detail {

template <typename Directed, typename Vertex>
struct edge_base {
  inline edge_base() = default;
  inline edge_base(Vertex s, Vertex d) : m_source(s), m_target(d) {}
  Vertex m_source;
  Vertex m_target;
};

template <typename Directed, typename Vertex>
class edge_desc_impl : public edge_base<Directed, Vertex> {
  using self = edge_desc_impl;
  using Base = edge_base<Directed, Vertex>;

 public:
  using property_type = void;

  inline edge_desc_impl() : m_eproperty(0) {}

  inline edge_desc_impl(Vertex s, Vertex d, const property_type* eplug)
      : Base(s, d), m_eproperty(const_cast<property_type*>(eplug)) {}

  property_type* get_property() { return m_eproperty; }
  const property_type* get_property() const { return m_eproperty; }

  //  protected:
  property_type* m_eproperty;
};

template <class D, class V>
inline bool operator==(const edge_desc_impl<D, V>& a, const edge_desc_impl<D, V>& b) {
  return a.get_property() == b.get_property();
}
template <class D, class V>
inline bool operator!=(const edge_desc_impl<D, V>& a, const edge_desc_impl<D, V>& b) {
  return !(a.get_property() == b.get_property());
}

// Order edges according to the address of their property object
template <class D, class V>
inline bool operator<(const edge_desc_impl<D, V>& a, const edge_desc_impl<D, V>& b) {
  return a.get_property() < b.get_property();
}
template <class D, class V>
inline bool operator<=(const edge_desc_impl<D, V>& a, const edge_desc_impl<D, V>& b) {
  return a.get_property() <= b.get_property();
}
template <class D, class V>
inline bool operator>(const edge_desc_impl<D, V>& a, const edge_desc_impl<D, V>& b) {
  return a.get_property() > b.get_property();
}
template <class D, class V>
inline bool operator>=(const edge_desc_impl<D, V>& a, const edge_desc_impl<D, V>& b) {
  return a.get_property() >= b.get_property();
}

template <class Char, class Traits, class D, class V>
std::basic_ostream<Char, Traits>& operator<<(std::basic_ostream<Char, Traits>& os, const edge_desc_impl<D, V>& e) {
  return os << "(" << e.m_source << "," << e.m_target << ")";
}

}  // namespace bagl::detail

// Standard hash
template <typename D, typename V>
struct std::hash<bagl::detail::edge_desc_impl<D, V> > {
  std::size_t operator()(const bagl::detail::edge_desc_impl<D, V>& x) const { return hash_value(x.get_property()); }
};

#endif  // BAGL_BAGL_DETAIL_EDGE_H_
