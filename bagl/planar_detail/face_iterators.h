// Copyright (c) Aaron Windsor 2007
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLANAR_DETAIL_FACE_ITERATORS_H_
#define BAGL_BAGL_PLANAR_DETAIL_FACE_ITERATORS_H_

#include <iterator>
#include <ranges>
#include <type_traits>

#include "bagl/graph_traits.h"

namespace bagl {

// tags for defining traversal properties

// VisitorType
struct lead_visitor {};
struct follow_visitor {};

// TraversalType
struct single_side {};
struct both_sides {};

// TraversalSubType
struct first_side {};   // for single_side
struct second_side {};  // for single_side
struct alternating {};  // for both_sides

// Time
struct current_iteration {};
struct previous_iteration {};

// Why TraversalType AND TraversalSubType? TraversalSubType is a function
// template parameter passed in to the constructor of the face iterator,
// whereas TraversalType is a class template parameter. This lets us decide
// at runtime whether to move along the first or second side of a bicomp (by
// assigning a face_iterator that has been constructed with TraversalSubType
// = first_side or second_side to a face_iterator variable) without any of
// the virtual function overhead that comes with implementing this
// functionality as a more structured form of type erasure. It also allows
// a single face_iterator to be the end iterator of two iterators traversing
// both sides of a bicomp.

// ValueType is either graph_vertex_descriptor_t<Graph> or graph_edge_descriptor_t<Graph>

// forward declaration (defining defaults)
template <typename Graph, typename FaceHandlesMap, typename ValueType, typename TraversalType = single_side,
          typename VisitorType = lead_visitor, typename Time = current_iteration>
class face_iterator;

// specialization for TraversalType = traverse_vertices
template <typename Graph, typename FaceHandlesMap, typename ValueType, typename TraversalType, typename VisitorType,
          typename Time>
class face_iterator {
 public:
  using Vertex = graph_vertex_descriptor_t<Graph>;
  using edge_t = graph_edge_descriptor_t<Graph>;
  using self = face_iterator<Graph, FaceHandlesMap, ValueType, TraversalType, VisitorType, Time>;
  using face_handle_t = typename FaceHandlesMap::value_type;

  face_iterator() : lead_(graph_traits<Graph>::null_vertex()), follow_(graph_traits<Graph>::null_vertex()) {}

  template <typename TraversalSubType>
  face_iterator(face_handle_t anchor_handle, FaceHandlesMap face_handles, TraversalSubType traversal_type)
      : follow_(anchor_handle.get_anchor()), face_handles_(face_handles) {
    set_lead_dispatch(anchor_handle, traversal_type);
  }

  template <typename TraversalSubType>
  face_iterator(Vertex anchor, FaceHandlesMap face_handles, TraversalSubType traversal_type)
      : follow_(anchor), face_handles_(face_handles) {
    set_lead_dispatch(face_handles_[anchor], traversal_type);
  }

  using value_type = ValueType;
  using reference = ValueType;
  using pointer = const ValueType*;
  using difference_type = int;
  using iterator_category = std::forward_iterator_tag;

  bool operator==(const self& rhs) { return equal(rhs); }
  bool operator!=(const self& rhs) { return !equal(rhs); }

  self& operator++() {
    increment();
    return *this;
  }
  self operator++(int) {
    self tmp = *this;
    increment();
    return tmp;
  }

  auto operator*() const { return dereference(); }
  auto operator->() const { return &dereference(); }

 private:
  static Vertex get_second_vertex(face_handle_t anchor_handle) {
    if constexpr (std::is_same_v<Time, previous_iteration>) {
      return anchor_handle.old_second_vertex();
    } else {
      return anchor_handle.second_vertex();
    }
  }

  static Vertex get_first_vertex(face_handle_t anchor_handle) {
    if constexpr (std::is_same_v<Time, previous_iteration>) {
      return anchor_handle.old_first_vertex();
    } else {
      return anchor_handle.first_vertex();
    }
  }

  static auto get_edge_to_second(face_handle_t anchor_handle) {
    if constexpr (std::is_same_v<ValueType, edge_t>) {
      if constexpr (std::is_same_v<Time, previous_iteration>) {
        return anchor_handle.old_second_edge();
      } else {
        return anchor_handle.second_edge();
      }
    } else {
      return false;
    }
  }

  static auto get_edge_to_first(face_handle_t anchor_handle) {
    if constexpr (std::is_same_v<ValueType, edge_t>) {
      if constexpr (std::is_same_v<Time, previous_iteration>) {
        return anchor_handle.old_first_edge();
      } else {
        return anchor_handle.first_edge();
      }
    } else {
      return false;
    }
  }

  void set_lead_dispatch(face_handle_t anchor_handle, first_side /*unused*/) {
    lead_ = get_first_vertex(anchor_handle);
    edge_ = get_edge_to_first(anchor_handle);
  }

  void set_lead_dispatch(face_handle_t anchor_handle, second_side /*unused*/) {
    lead_ = get_second_vertex(anchor_handle);
    edge_ = get_edge_to_second(anchor_handle);
  }

  void increment() {
    face_handle_t curr_face_handle(face_handles_[lead_]);
    Vertex first = get_first_vertex(curr_face_handle);
    Vertex second = get_second_vertex(curr_face_handle);
    if (first == follow_) {
      follow_ = lead_;
      edge_ = get_edge_to_second(curr_face_handle);
      lead_ = second;
    } else if (second == follow_) {
      follow_ = lead_;
      edge_ = get_edge_to_first(curr_face_handle);
      lead_ = first;
    } else {
      lead_ = follow_ = graph_traits<Graph>::null_vertex();
    }
  }

  bool equal(self const& other) const { return lead_ == other.lead_ && follow_ == other.follow_; }

  auto dereference() const {
    if constexpr (std::is_same_v<ValueType, edge_t>) {
      return edge_;
    } else if constexpr (std::is_same_v<VisitorType, lead_visitor>) {
      return lead_;
    } else {
      return follow_;
    }
  }

  Vertex lead_;
  Vertex follow_;
  std::conditional_t<std::is_same_v<ValueType, edge_t>, edge_t, bool> edge_;
  FaceHandlesMap face_handles_;
};

template <typename Graph, typename FaceHandlesMap, typename ValueType, typename VisitorType, typename Time>
class face_iterator<Graph, FaceHandlesMap, ValueType, both_sides, VisitorType, Time> {
 public:
  using self = face_iterator<Graph, FaceHandlesMap, ValueType, both_sides, VisitorType, Time>;
  using Vertex = typename graph_traits<Graph>::vertex_descriptor;
  using face_handle_t = typename FaceHandlesMap::value_type;

  face_iterator() = default;

  face_iterator(face_handle_t anchor_handle, FaceHandlesMap face_handles)
      : first_itr_(anchor_handle, face_handles, first_side()),
        second_itr_(anchor_handle, face_handles, second_side()) {}

  face_iterator(Vertex anchor, FaceHandlesMap face_handles)
      : first_itr_(face_handles[anchor], face_handles, first_side()),
        second_itr_(face_handles[anchor], face_handles, second_side()) {}

  using value_type = ValueType;
  using reference = ValueType;
  using pointer = const ValueType*;
  using difference_type = int;
  using iterator_category = std::forward_iterator_tag;

  bool operator==(const self& rhs) { return equal(rhs); }
  bool operator!=(const self& rhs) { return !equal(rhs); }

  self& operator++() {
    increment();
    return *this;
  }
  self operator++(int) {
    self tmp = *this;
    increment();
    return tmp;
  }

  auto operator*() const { return dereference(); }
  auto operator->() const { return &dereference(); }

 private:
  using inner_itr_t = face_iterator<Graph, FaceHandlesMap, ValueType, single_side, follow_visitor, Time>;

  void increment() {
    if (first_increment_) {
      ++first_itr_;
      ++second_itr_;
      first_increment_ = false;
    } else if (first_is_active_) {
      ++first_itr_;
    } else {
      ++second_itr_;
    }
    first_is_active_ = !first_is_active_;
  }

  bool equal(self const& other) const {
    // Want this iterator to be equal to the "end" iterator when at least
    // one of the iterators has reached the root of the current bicomp.
    // This isn't ideal, but it works.

    return (first_itr_ == other.first_itr_ || second_itr_ == other.second_itr_);
  }

  ValueType dereference() const { return first_is_active_ ? *first_itr_ : *second_itr_; }

  inner_itr_t first_itr_;
  inner_itr_t second_itr_;
  inner_itr_t face_end_;
  bool first_is_active_ = true;
  bool first_increment_ = true;
};

}  // namespace bagl

#endif  // BAGL_BAGL_PLANAR_DETAIL_FACE_ITERATORS_H_
