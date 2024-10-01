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
template <typename IndexType>
class component_index {
 private:
  using IndexContainer = std::vector<IndexType>;

 public:
  using component_iterator = inc_comp_detail::component_index_iterator<typename IndexContainer::iterator>;

  template <typename ParentIterator, typename ElementIndexMap>
  component_index(ParentIterator parent_start, ParentIterator parent_end, const ElementIndexMap& index_map)
      : num_elements_(std::distance(parent_start, parent_end)),
        components_(std::make_shared<IndexContainer>()),
        index_list_(std::make_shared<IndexContainer>(num_elements_)) {
    build_index_lists(parent_start, index_map);

  }  // component_index

  template <typename ParentIterator>
  component_index(ParentIterator parent_start, ParentIterator parent_end)
      : num_elements_(std::distance(parent_start, parent_end)),
        components_(std::make_shared<IndexContainer>()),
        index_list_(std::make_shared<IndexContainer>(num_elements_)) {
    build_index_lists(parent_start, identity_property_map());

  }  // component_index

  // Returns the number of components
  [[nodiscard]] std::size_t size() const { return (components_->size()); }

  // Range of component indices
  auto components() const { return std::ranges::iota_view(IndexType{0}, IndexType{components_->size()}); }

  // Returns a range for the child elements of component [component_index].
  auto component(IndexType component_index) const {
    IndexType first_index = (*components_)[component_index];
    return std::ranges::subrange(component_iterator(index_list_->begin(), first_index),
                                 component_iterator(num_elements_));
  }

 private:
  template <typename ParentIterator, typename ElementIndexMap>
  void build_index_lists(ParentIterator parent_start, const ElementIndexMap& index_map) {
    auto index_list = index_list_->begin();

    // First pass - find root elements, construct index list
    for (IndexType element_index = 0; element_index < num_elements_; ++element_index) {
      auto parent_element = parent_start[element_index];
      IndexType parent_index = get(index_map, parent_element);

      if (element_index != parent_index) {
        index_list[element_index] = parent_index;
      } else {
        components_->push_back(element_index);

        // m_num_elements is the linked list terminator
        index_list[element_index] = num_elements_;
      }
    }

    // Second pass - build linked list
    for (IndexType element_index = 0; element_index < num_elements_; ++element_index) {
      auto parent_element = parent_start[element_index];
      IndexType parent_index = get(index_map, parent_element);

      if (element_index != parent_index) {
        // Follow list until a component parent is found
        while (index_list[parent_index] != num_elements_) {
          parent_index = index_list[parent_index];
        }

        // Push element to the front of the linked list
        index_list[element_index] = index_list[parent_index];
        index_list[parent_index] = element_index;
      }
    }

  }  // build_index_lists

 protected:
  IndexType num_elements_;
  std::shared_ptr<IndexContainer> components_;
  std::shared_ptr<IndexContainer> index_list_;

};  // class component_index

}  // namespace bagl

#endif  // BAGL_BAGL_INCREMENTAL_COMPONENTS_H_
