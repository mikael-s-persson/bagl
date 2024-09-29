// Copyright 2010 Matthias Walter (xammy@xammy.homelinux.net)
// Authors: Matthias Walter
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_BIPARTITE_H_
#define BAGL_BAGL_BIPARTITE_H_

#include <exception>
#include <iterator>
#include <ranges>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/depth_first_search.h"
#include "bagl/one_bit_color_map.h"
#include "bagl/properties.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace bipartite_detail {

// The bipartite_visitor_error is thrown if an edge cannot be colored.
// The witnesses are the edges incident vertices.
template <typename Vertex>
struct bipartite_visitor_error : std::exception {
  std::pair<Vertex, Vertex> witnesses;

  bipartite_visitor_error(Vertex a, Vertex b) : witnesses(a, b) {}

  [[nodiscard]] const char* what() const noexcept override { return "Graph is not bipartite."; }
};

// Functor which colors edges to be non-monochromatic.
template <typename PartitionMap>
struct bipartition_visitor {
  explicit bipartition_visitor(PartitionMap partition_map) : partition_map_(partition_map) {}

  template <typename Vertex, typename Graph>
  void start_vertex(Vertex v, const Graph&) {
    using color_traits = color_traits<property_traits_value_t<PartitionMap>>;
    put(partition_map_, v, color_traits::white());
  }

  template <typename Edge, typename Graph>
  void tree_edge(Edge e, const Graph& g) {
    using color_traits = color_traits<property_traits_value_t<PartitionMap>>;

    auto u = source(e, g);
    auto v = target(e, g);
    if (get(partition_map_, u) == color_traits::white()) {
      put(partition_map_, v, color_traits::black());
    } else {
      put(partition_map_, v, color_traits::white());
    }
  }

  template <typename Edge, typename Graph>
  void back_edge(Edge e, const Graph& g) {
    auto u = source(e, g);
    auto v = target(e, g);
    if (get(partition_map_, u) == get(partition_map_, v)) {
      throw bipartite_visitor_error(u, v);
    }
  }

 private:
  PartitionMap partition_map_;
};

// Find the beginning of a common suffix of two sequences
// `sequence1` Pair of bidirectional iterators defining the first sequence.
// `sequence2` Pair of bidirectional iterators defining the second sequence.
// Returns a pair of iterators pointing to the beginning of the common suffix.
template <std::ranges::bidirectional_range Range1, std::ranges::bidirectional_range Range2>
std::pair<std::ranges::iterator_t<Range1>, std::ranges::iterator_t<Range2>> reverse_mismatch(const Range1& rg1,
                                                                                             const Range2& rg2) {
  if (rg1.empty() || rg2.empty()) {
    return {rg1.begin(), rg2.begin()};
  }

  auto it1 = rg1.end();
  auto it2 = rg2.end();

  while (true) {
    --it1;
    --it2;
    if (*it1 != *it2) {
      ++it1;
      ++it2;
      break;
    }
    if (it1 == rg1.begin()) {
      break;
    }
    if (it2 == rg2.begin()) {
      break;
    }
  }

  return {it1, it2};
}

}  // namespace bipartite_detail

// Checks a given graph for bipartiteness and fills the given color map with
// white and black according to the bipartition. If the graph is not
// bipartite, the contents of the color map are undefined. Runs in linear
// time in the size of the graph, if access to the property maps is in
// constant time.
//
// `index_map` An index map associating vertices with an index.
// `partition_map` A color map to fill with the bipartition.
// Returns true if and only if the given graph is bipartite.
template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> IndexMap,
          concepts::ReadWriteVertexPropertyMap<G> PartitionMap>
bool is_bipartite(const G& g, const IndexMap index_map, PartitionMap partition_map) {
  // General types and variables
  using Vertex = graph_vertex_descriptor_t<G>;

  // Call dfs
  try {
    depth_first_search(g, make_dfs_visitor(bipartition_visitor(partition_map)),
                       make_vector_property_map(num_vertices_or_zero(g), index_map, default_color_type::white_color));
  } catch (const bipartite_detail::bipartite_visitor_error<Vertex>&) {
    return false;
  }

  return true;
}

// Checks a given graph for bipartiteness.
//
// `index_map` An index map associating vertices with an index.
// Returns true if and only if the given graph is bipartite.
template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> IndexMap>
bool is_bipartite(const G& g, const IndexMap index_map) {
  return is_bipartite(g, index_map, make_one_bit_color_map(num_vertices_or_zero(g), index_map));
}

// Checks a given graph for bipartiteness. The graph must
// have an internal vertex_index property. Runs in linear time in the
// size of the graph, if access to the property maps is in constant time.
//
// Returns true if and only if the given graph is bipartite.
template <concepts::VertexListGraph G>
bool is_bipartite(const G& g) {
  return is_bipartite(g, get(vertex_index, g));
}

// Checks a given graph for bipartiteness and fills a given color map with
// white and black according to the bipartition. If the graph is not
// bipartite, a sequence of vertices, producing an odd-cycle, is written to
// the output iterator. The final iterator value is returned. Runs in linear
// time in the size of the graph, if access to the property maps is in
// constant time.
//
// `index_map` An index map associating vertices with an index.
// `partition_map` A color map to fill with the bipartition.
// `result` An iterator to write the odd-cycle vertices to.
// Returns the final iterator value after writing.
template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> IndexMap,
          concepts::ReadWriteVertexPropertyMap<G> PartitionMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
OutputIterator find_odd_cycle(const G& g, IndexMap index_map, PartitionMap partition_map, OutputIterator result) {
  /// General types and variables
  using Vertex = graph_vertex_descriptor_t<G>;

  // Declare predecessor map
  auto predecessor_map = make_vector_property_map(num_vertices(g), index_map, graph_traits<G>::null_vertex());

  // Initialize predecessor map
  for (auto v : vertices(g)) {
    put(predecessor_map, v, v);
  }

  // Call dfs
  try {
    depth_first_search(
        g, make_dfs_visitor(bipartition_visitor(partition_map), predecessor_recorder_on_tree_edge(predecessor_map)),
        make_vector_property_map(num_vertices_or_zero(g), index_map, default_color_type::white_color));
  } catch (const bipartite_detail::bipartite_visitor_error<Vertex>& error) {
    using Path = std::vector<Vertex>;

    Path path1;
    Path path2;

    // First path
    Vertex next = error.witnesses.first;
    Vertex current = next;
    do {
      current = next;
      path1.push_back(current);
      next = predecessor_map[current];
    } while (current != next);

    // Second path
    next = error.witnesses.second;
    do {
      current = next;
      path2.push_back(current);
      next = predecessor_map[current];
    } while (current != next);

    // Find beginning of common suffix
    auto mismatch = bipartite_detail::reverse_mismatch(path1, path2);

    // Copy the odd-length cycle
    result = std::copy(path1.begin(), mismatch.first + 1, result);
    return std::reverse_copy(path2.begin(), mismatch.second, result);
  }

  return result;
}

// Checks a given graph for bipartiteness. If the graph is not bipartite, a
// sequence of vertices, producing an odd-cycle, is written to the output
// iterator. The final iterator value is returned. Runs in linear time in the
// size of the graph, if access to the property maps is in constant time.
//
// `index_map` An index map associating vertices with an index.
// `result` An iterator to write the odd-cycle vertices to.
// Returns the final iterator value after writing.
template <concepts::VertexListGraph G, concepts::ReadableVertexPropertyMap<G> IndexMap,
          std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
OutputIterator find_odd_cycle(const G& g, const IndexMap index_map, OutputIterator result) {
  return find_odd_cycle(g, index_map, make_one_bit_color_map(num_vertices(g), index_map), result);
}

// Checks a given graph for bipartiteness. If the graph is not bipartite, a
// sequence of vertices, producing an odd-cycle, is written to the output
// iterator. The final iterator value is returned. The graph must have an
// internal vertex_index property. Runs in linear time in the size of the
// graph, if access to the property maps is in constant time.
//
// `result` An iterator to write the odd-cycle vertices to.
// Returns the final iterator value after writing.
template <concepts::VertexListGraph G, std::output_iterator<graph_vertex_descriptor_t<G>> OutputIterator>
OutputIterator find_odd_cycle(const G& g, OutputIterator result) {
  return find_odd_cycle(g, get(vertex_index, g), result);
}

}  // namespace bagl

#endif  // BAGL_BAGL_BIPARTITE_H_
