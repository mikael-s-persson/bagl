// Copyright 2010 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_RANDOM_SPANNING_TREE_H_
#define BAGL_BAGL_RANDOM_SPANNING_TREE_H_

#include <cassert>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/loop_erased_random_walk.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace random_detail {
// Use Wilson's algorithm (based on loop-free random walks) to generate a
// random spanning tree.  The distribution of edges used is controlled by
// the next_edge() function, so this version allows either weighted or
// unweighted selection of trees.
// Algorithm is from http://en.wikipedia.org/wiki/Uniform_spanning_tree
template <concepts::VertexListGraph Graph, concepts::WritableVertexPropertyMap<Graph> PredMap,
          concepts::ReadWriteVertexPropertyMap<Graph> ColorMap, typename NextEdge>
void random_spanning_tree_internal(const Graph& g, graph_vertex_descriptor_t<Graph> s, PredMap pred, ColorMap color,
                                   NextEdge next_edge) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // g must also be undirected (or symmetric) and connected
  assert(num_vertices(g) >= 1);

  using Color = color_traits<property_traits_value_t<ColorMap>>;
  for (auto v : vertices(g)) {
    put(color, v, Color::white());
  }

  std::vector<Vertex> path;

  put(color, s, Color::black());
  put(pred, s, graph_traits<Graph>::null_vertex());

  for (auto v : vertices(g)) {
    if (get(color, v) != Color::white()) {
      continue;
    }
    loop_erased_random_walk(g, v, next_edge, color, path);
    for (auto i = path.rbegin(); std::next(i) != path.rend(); ++i) {
      auto j = std::next(i);
      assert(get(color, *j) == Color::gray());
      put(color, *j, Color::black());
      put(pred, *j, *i);
    }
  }
}
}  // namespace random_detail

// Compute a uniformly-distributed spanning tree on a graph.  Use Wilson's
// algorithm:
// @inproceedings{wilson96generating,
//    author = {Wilson, David Bruce},
//    title = {Generating random spanning trees more quickly than the cover
//    time}, booktitle = {STOC '96: Proceedings of the twenty-eighth annual ACM
//    symposium on Theory of computing}, year = {1996}, isbn = {0-89791-785-5},
//    pages = {296--303},
//    location = {Philadelphia, Pennsylvania, United States},
//    doi = {http://doi.acm.org/10.1145/237814.237880},
//    publisher = {ACM},
//    address = {New York, NY, USA},
//  }
//
template <concepts::VertexListGraph Graph, std::uniform_random_bit_generator Gen,
          concepts::WritableVertexPropertyMap<Graph> PredMap, concepts::ReadWriteVertexPropertyMap<Graph> ColorMap>
void random_spanning_tree(const Graph& g, Gen& gen, graph_vertex_descriptor_t<Graph> root, PredMap pred,
                          ColorMap color) {
  unweighted_random_out_edge_gen<Graph, Gen> random_oe(gen);
  random_detail::random_spanning_tree_internal(g, root, pred, color, random_oe);
}

template <concepts::VertexListGraph Graph, std::uniform_random_bit_generator Gen,
          concepts::WritableVertexPropertyMap<Graph> PredMap>
void random_spanning_tree(const Graph& g, Gen& gen, graph_vertex_descriptor_t<Graph> root, PredMap pred) {
  random_spanning_tree(
      g, gen, root, pred,
      vector_property_map(num_vertices(g), get(vertex_index, g), default_color_type::white_color).ref());
}

// Compute a weight-distributed spanning tree on a graph.
template <concepts::VertexListGraph Graph, std::uniform_random_bit_generator Gen,
          concepts::WritableVertexPropertyMap<Graph> PredMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap,
          concepts::ReadWriteVertexPropertyMap<Graph> ColorMap>
void random_spanning_tree(const Graph& g, Gen& gen, graph_vertex_descriptor_t<Graph> root, PredMap pred,
                          WeightMap weight, ColorMap color) {
  weighted_random_out_edge_gen<Graph, WeightMap, Gen> random_oe(weight, gen);
  random_detail::random_spanning_tree_internal(g, root, pred, color, random_oe);
}
template <concepts::VertexListGraph Graph, std::uniform_random_bit_generator Gen,
          concepts::WritableVertexPropertyMap<Graph> PredMap, concepts::ReadableEdgePropertyMap<Graph> WeightMap>
void random_spanning_tree(const Graph& g, Gen& gen, graph_vertex_descriptor_t<Graph> root, PredMap pred,
                          WeightMap weight) {
  random_spanning_tree(
      g, gen, root, pred, weight,
      vector_property_map(num_vertices(g), get(vertex_index, g), default_color_type::white_color).ref());
}

}  // namespace bagl

#endif  // BAGL_BAGL_RANDOM_SPANNING_TREE_H_
