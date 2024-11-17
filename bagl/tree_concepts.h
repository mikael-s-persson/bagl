// Copyright 2013-2024 Mikael Persson

#ifndef BAGL_BAGL_TREE_CONCEPTS_H_
#define BAGL_BAGL_TREE_CONCEPTS_H_

#include <utility>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/tree_adaptor.h"
#include "bagl/tree_traits.h"

namespace bagl::concepts {

// This library provides concept classes to verify that a data structure has a tree interface,
// i.e., models the concepts of trees as used in the Boost Graph Library. The tree concepts
// should be regarded as a special kind of graph and reuses much of the general graph traits.

template <typename T>
concept Tree = IncidenceGraph<T> && requires(const T& t, graph_vertex_descriptor_t<T> v) {
  { tree_root(t) } -> std::convertible_to<graph_vertex_descriptor_t<T>>;
  { children(v, t) } -> VertexRange<T>;
};

template <typename T>
concept BidirectionalTree = Tree<T> && requires(const T& t, graph_vertex_descriptor_t<T> v) {
  { parent(v, t) } -> std::convertible_to<graph_vertex_descriptor_t<T>>;
};

template <typename T>
concept MutableTree = Tree<T> && requires(T& t, graph_vertex_descriptor_t<T> v) {
  { create_root(t) } -> std::convertible_to<graph_vertex_descriptor_t<T>>;
  {
    add_child(v, t)
    } -> std::convertible_to<std::tuple<graph_vertex_descriptor_t<T>, graph_edge_descriptor_t<T>, bool>>;
  remove_branch(v, t);   // Removes `v` and all sub-tree
  clear_children(v, t);  // Removes sub-tree below `v`, but keeps `v`
};

// This concept defines the requirements to fulfill in order to model a mutable property-tree.
// A mutable property-tree is a mutable tree whose mutating functions take or deliver the
// vertex- or edge- property values associated with the vertices or edges. During removal
// of a branch, all the vertex-properties are collected into an output iterator (e.g., back-inserter).
// During additions of child nodes, the corresponding vertex-properties can be used to initialize
// the new vertex and edge directly. This not only makes such a tree easier to use (not having to
// manually collect or set vertex properties before or after the mutation), but it can also be
// necessary in some situations. One typical use-case is when re-balancing a branch of a tree,
// which often results in collecting properties of vertices (to preserve them), clearing the
// branch, re-adding the vertices in a new arrangement, and restoring their original properties.
//
// Valid expressions (vertex u, v; edge e; TreeType tree; vertex_property_type vp; edge_property_type ep; ):
//
// The root vertex of the tree can be created with a given vertex-property value (vp):
//   v = create_root(vp, g);
//   v = create_root(std::move(vp), g);
//
// A child vertex (v) and its incident edge (e) can be created for a given parent vertex (u) and with the properties:
//   tie(v,e,b) = add_child(u, vp, g);
//   tie(v,e,b) = add_child(u, vp, ep, g);
//   tie(v,e,b) = add_child(u, std::move(vp), g);
//   tie(v,e,b) = add_child(u, std::move(vp), std::move(ep), g);
//
// The sub-tree below a given vertex, INCLUDING the vertex, can be removed from the tree, with properies moved to output
// iterators:
//   remove_branch(v, vp_out_iter, g);
//   remove_branch(v, vp_out_iter, ep_out_iter, g);
//
// The sub-tree below a given vertex, EXCLUDING the vertex, can be removed from the tree, with properies moved to output
// iterators:
//   clear_children(v, vp_out_iter, g);
//   clear_children(v, vp_out_iter, Ep_out_iter, g);
template <typename T>
concept MutablePropertyTree = Tree<T> &&
    requires(T& t, graph_vertex_descriptor_t<T> v, vertex_property_type<T> vp, edge_property_type<T> ep) {
  { create_root(t, vp) } -> std::convertible_to<graph_vertex_descriptor_t<T>>;
  {
    add_child(v, t, vp)
    } -> std::convertible_to<std::tuple<graph_vertex_descriptor_t<T>, graph_edge_descriptor_t<T>, bool>>;
  {
    add_child(v, t, vp, ep)
    } -> std::convertible_to<std::tuple<graph_vertex_descriptor_t<T>, graph_edge_descriptor_t<T>, bool>>;
  { remove_branch(v, t, &vp) } -> std::convertible_to<vertex_property_type<T>*>;
  { remove_branch(v, t, &vp, &ep) } -> std::convertible_to<std::pair<vertex_property_type<T>*, edge_property_type<T>*>>;
  { clear_children(v, t, &vp) } -> std::convertible_to<vertex_property_type<T>*>;
  {
    clear_children(v, t, &vp, &ep)
    } -> std::convertible_to<std::pair<vertex_property_type<T>*, edge_property_type<T>*>>;
};

}  // namespace bagl::concepts

#endif  // BAGL_BAGL_TREE_CONCEPTS_H_
