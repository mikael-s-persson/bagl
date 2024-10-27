// Copyright 1997-2001 University of Notre Dame.
// Copyright 2009 Trustees of Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek, Michael Hansen
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_INCREMENTAL_COMPONENTS_H_
#define BAGL_BAGL_INCREMENTAL_COMPONENTS_H_

#include <iterator>
#include <memory>
#include <tuple>

#include "bagl/detail/incremental_components.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "detail/disjoint_sets.h"

namespace bagl {

// A connected component algorithm for the case when dynamically
// adding (but not removing) edges is common.  The
// incremental_components() function is a preparing operation. Call
// same_component to check whether two vertices are in the same
// component, or use disjoint_set::find_set to determine the
// representative for a vertex.

// This version of connected components does not require a full
// Graph. Instead, it just needs an edge list, where the vertices of
// each edge need to be of integer type. The edges are assumed to
// be undirected. The other difference is that the result is stored in
// a container, instead of just a decorator.  The container should be
// empty before the algorithm is called. It will grow during the
// course of the algorithm. The container must be a model of
// BackInsertionSequence and RandomAccessContainer
// (std::vector is a good choice). After running the algorithm the
// index container will map each vertex to the representative
// vertex of the component to which it belongs.
//
// Adapted from an implementation by Alex Stepanov. The disjoint
// sets data structure is from Tarjan's "Data Structures and Network
// Algorithms", and the application to connected components is
// similar to the algorithm described in Ch. 22 of "Intro to
// Algorithms" by Cormen, et. all.
//

// An implementation of disjoint sets can be found in
// boost/pending/disjoint_sets.hpp

template <concepts::EdgeListGraph G, class DisjointSets>
void incremental_components(const G& g, DisjointSets& ds) {
  for (auto e : edges(g)) {
    ds.union_set(source(e, g), target(e, g));
  }
}

template <std::ranges::random_access_range ParentRange>
void compress_components(const ParentRange& p_rg) {
  for (auto current = p_rg.begin(); current != p_rg.end(); ++current) {
    disjoint_sets_detail::find_representative_with_full_compression(p_rg.begin(), current - p_rg.begin());
  }
}

template <std::ranges::random_access_range ParentRange>
auto component_count(const ParentRange& p_rg) {
  std::ptrdiff_t count = 0;
  for (auto current = p_rg.begin(); current != p_rg.end(); ++current) {
    if (*current == current - p_rg.begin()) {
      ++count;
    }
  }
  return count;
}

// This algorithm can be applied to the result container of the
// connected_components algorithm to normalize
// the components.
template <std::ranges::random_access_range ParentRange>
void normalize_components(const ParentRange& p_rg) {
  for (auto current = p_rg.begin(); current != p_rg.end(); ++current) {
    disjoint_sets_detail::normalize_node(p_rg.begin(), current - p_rg.begin());
  }
}

template <concepts::VertexListGraph G, class DisjointSets>
void initialize_incremental_components(G& g, DisjointSets& ds) {
  for (auto v : vertices(g)) {
    ds.make_set(v);
  }
}

template <class Vertex, class DisjointSet>
bool same_component(Vertex u, Vertex v, DisjointSet& ds) {
  return ds.find_set(u) == ds.find_set(v);
}

// Class that builds a quick-access indexed linked list that allows
// for fast iterating through a parent component's children.
class component_index {
 private:
  using IndexContainer = std::vector<std::size_t>;

 public:
  using component_iterator = inc_comp_detail::component_index_iterator<typename IndexContainer::const_iterator>;

  template <typename ParentMap, typename VertexIndexMap, typename VertexRange>
  component_index(ParentMap parent, VertexIndexMap vindex, VertexRange vrange) {
    build_index_lists(parent, vindex, vrange);

  }  // component_index

  template <typename ParentMap, typename VertexRange>
  component_index(ParentMap parent, VertexRange vrange) {
    build_index_lists(parent, identity_property_map(), vrange);

  }  // component_index

  // Returns the number of components
  [[nodiscard]] std::size_t size() const { return (components_.size()); }

  // Range of component indices
  auto components() const { return std::ranges::iota_view(std::size_t{0}, std::size_t{components_.size()}); }

  // Returns a range for the child elements of component [component_index].
  auto component(std::size_t component_index) const {
    std::size_t first_index = components_[component_index];
    return std::ranges::subrange(component_iterator(index_list_.begin(), first_index),
                                 component_iterator(index_list_.size()));
  }

 private:
  template <typename ParentMap, typename VertexIndexMap, typename VertexRange>
  void build_index_lists(ParentMap parent, VertexIndexMap vindex, VertexRange vrange) {
    std::size_t v_sz = std::distance(vrange.begin(), vrange.end());
    index_list_.resize(v_sz);

    // First pass - find root elements, construct index list
    for (auto v : vrange) {
      auto v_index = get(vindex, v);
      auto u = parent[v];
      auto u_index = get(vindex, u);

      if (v_index != u_index) {
        index_list_[v_index] = u_index;
      } else {
        components_.push_back(v_index);
        // size is the linked list terminator
        index_list_[v_index] = v_sz;
      }
    }

    // Second pass - build linked list
    for (auto v : vrange) {
      auto v_index = get(vindex, v);
      auto u = parent[v];
      auto u_index = get(vindex, u);

      if (v_index != u_index) {
        // Follow list until a component parent is found
        while (index_list_[u_index] != v_sz) {
          u_index = index_list_[u_index];
        }

        // Push element to the front of the linked list
        index_list_[v_index] = index_list_[u_index];
        index_list_[u_index] = v_index;
      }
    }
  }

 protected:
  std::vector<std::size_t> components_;
  std::vector<std::size_t> index_list_;

};  // class component_index

}  // namespace bagl

#endif  // BAGL_BAGL_INCREMENTAL_COMPONENTS_H_
