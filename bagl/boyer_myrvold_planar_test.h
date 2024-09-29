// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BOYER_MYRVOLD_PLANAR_TEST_H_
#define BAGL_BAGL_BOYER_MYRVOLD_PLANAR_TEST_H_

#include <type_traits>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/planar_detail/boyer_myrvold_impl.h"
#include "bagl/property_map.h"

namespace bagl {

// Basic versions.

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, VertexIndexMap v_index) {
  boyer_myrvold_impl<G, VertexIndexMap, planar_detail::no_old_handles, planar_detail::no_embedding> planarity_tester(
      g, v_index);

  return planarity_tester.is_planar();
}

template <typename G>
bool boyer_myrvold_planarity_test(const G& g) {
  return boyer_myrvold_planarity_test(g, get(vertex_index, g));
}

// Versions that gather the kuratowski subgraph.

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator,
          concepts::ReadableEdgePropertyMap<G> EdgeIndexMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, VertexIndexMap v_index, OutputIterator kuratowski_subgraph,
                                  EdgeIndexMap e_index) {
  boyer_myrvold_impl<G, VertexIndexMap, planar_detail::store_old_handles, planar_detail::no_embedding> planarity_tester(
      g, v_index);

  if (planarity_tester.is_planar()) {
    return true;
  }
  planarity_tester.extract_kuratowski_subgraph(kuratowski_subgraph, e_index);
  return false;
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, VertexIndexMap v_index, OutputIterator kuratowski_subgraph) {
  return boyer_myrvold_planarity_test(g, v_index, kuratowski_subgraph, get(edge_index, g));
}

template <concepts::VertexAndEdgeListGraph G, std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator,
          concepts::ReadableEdgePropertyMap<G> EdgeIndexMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, OutputIterator kuratowski_subgraph, EdgeIndexMap e_index) {
  return boyer_myrvold_planarity_test(g, get(vertex_index, g), kuratowski_subgraph, e_index);
}

template <concepts::VertexAndEdgeListGraph G, std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, OutputIterator kuratowski_subgraph) {
  return boyer_myrvold_planarity_test(g, get(vertex_index, g), kuratowski_subgraph, get(edge_index, g));
}

// Versions that gather the permutation map.

template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::MutableLvalueVertexPropertyMap<G> PermutationMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, VertexIndexMap v_index, PermutationMap embedding) {
  boyer_myrvold_impl<G, VertexIndexMap, planar_detail::no_old_handles, planar_detail::std_list> planarity_tester(
      g, v_index);

  if (planarity_tester.is_planar()) {
    planarity_tester.make_edge_permutation(embedding);
    return true;
  }
  return false;
}

template <concepts::VertexListGraph G, concepts::MutableLvalueVertexPropertyMap<G> PermutationMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, PermutationMap embedding) {
  return boyer_myrvold_planarity_test(g, get(vertex_index, g), embedding);
}

// Versions that gather both the kuratowski subgraph and permutation map.

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::MutableLvalueVertexPropertyMap<G> PermutationMap,
          std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator,
          concepts::ReadableEdgePropertyMap<G> EdgeIndexMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, VertexIndexMap v_index, PermutationMap embedding,
                                  OutputIterator kuratowski_subgraph, EdgeIndexMap e_index) {
  boyer_myrvold_impl<G, VertexIndexMap, planar_detail::store_old_handles, planar_detail::std_list> planarity_tester(
      g, v_index);

  if (planarity_tester.is_planar()) {
    planarity_tester.make_edge_permutation(embedding);
    return true;
  }
  planarity_tester.extract_kuratowski_subgraph(kuratowski_subgraph, e_index);
  return false;
}

template <concepts::VertexAndEdgeListGraph G, concepts::MutableLvalueVertexPropertyMap<G> PermutationMap,
          std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator,
          concepts::ReadableEdgePropertyMap<G> EdgeIndexMap>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, PermutationMap embedding, OutputIterator kuratowski_subgraph,
                                  EdgeIndexMap e_index) {
  return boyer_myrvold_planarity_test(g, get(vertex_index, g), embedding, kuratowski_subgraph, e_index);
}

template <concepts::VertexAndEdgeListGraph G, concepts::ReadableVertexPropertyMap<G> VertexIndexMap,
          concepts::MutableLvalueVertexPropertyMap<G> PermutationMap,
          std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, VertexIndexMap v_index, PermutationMap embedding,
                                  OutputIterator kuratowski_subgraph) {
  return boyer_myrvold_planarity_test(g, v_index, embedding, kuratowski_subgraph, get(edge_index, g));
}

template <concepts::VertexAndEdgeListGraph G, concepts::MutableLvalueVertexPropertyMap<G> PermutationMap,
          std::output_iterator<graph_edge_descriptor_t<G>> OutputIterator>
requires concepts::IncidenceGraph<G>
bool boyer_myrvold_planarity_test(const G& g, PermutationMap embedding, OutputIterator kuratowski_subgraph) {
  return boyer_myrvold_planarity_test(g, get(vertex_index, g), embedding, kuratowski_subgraph, get(edge_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_BOYER_MYRVOLD_PLANAR_TEST_H_
