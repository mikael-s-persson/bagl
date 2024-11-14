// Copyright (c) Aaron Windsor 2007
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PLANAR_CANONICAL_ORDERING_H_
#define BAGL_BAGL_PLANAR_CANONICAL_ORDERING_H_

#include <iterator>
#include <list>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"

namespace bagl {

namespace planar_canonical_ordering_state {
constexpr int processed = 0;
constexpr int unprocessed = 1;
constexpr int one_neighbor_processed = 2;
constexpr int ready_to_be_processed = 3;
}  // namespace planar_canonical_ordering_state

template <concepts::VertexListGraph Graph, concepts::LvalueVertexPropertyMap<Graph> PlanarEmbedding,
          std::output_iterator<graph_vertex_descriptor_t<Graph>> OutputIterator,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::AdjacencyGraph<Graph>
void planar_canonical_ordering(const Graph& g, PlanarEmbedding embedding, OutputIterator ordering, VertexIndexMap vm) {
  using vertex_t = graph_vertex_descriptor_t<Graph>;
  using edge_t = graph_edge_descriptor_t<Graph>;

  auto processed_neighbor = vector_property_map(num_vertices(g), vm, graph_traits<Graph>::null_vertex());
  auto status = vector_property_map(num_vertices(g), vm, planar_canonical_ordering_state::unprocessed);

  std::list<vertex_t> ready_to_be_processed;

  vertex_t first_vertex = *vertices(g).begin();
  vertex_t second_vertex = first_vertex;
  for (auto u : adjacent_vertices(first_vertex, g)) {
    if (u != first_vertex) {
      second_vertex = u;
      break;
    }
  }

  ready_to_be_processed.push_back(first_vertex);
  status[first_vertex] = planar_canonical_ordering_state::ready_to_be_processed;
  ready_to_be_processed.push_back(second_vertex);
  status[second_vertex] = planar_canonical_ordering_state::ready_to_be_processed;

  while (!ready_to_be_processed.empty()) {
    vertex_t u = ready_to_be_processed.front();
    ready_to_be_processed.pop_front();

    if (status[u] != planar_canonical_ordering_state::ready_to_be_processed && u != second_vertex) {
      continue;
    }

    auto ei_start = embedding[u].begin();
    auto ei_end = embedding[u].end();
    auto prior_edge_itr = std::prev(ei_end);
    while (source(*prior_edge_itr, g) == target(*prior_edge_itr, g)) {
      prior_edge_itr = std::prev(prior_edge_itr);
    }

    for (auto ei = ei_start; ei != ei_end; ++ei) {
      edge_t e(*ei);  // e = (u,v)
      auto next_edge_itr = std::next(ei) == ei_end ? ei_start : std::next(ei);
      vertex_t v = (source(e, g) == u ? target(e, g) : source(e, g));

      vertex_t prior_vertex =
          (source(*prior_edge_itr, g) == u ? target(*prior_edge_itr, g) : source(*prior_edge_itr, g));
      vertex_t next_vertex = (source(*next_edge_itr, g) == u ? target(*next_edge_itr, g) : source(*next_edge_itr, g));

      // Need prior_vertex, u, v, and next_vertex to all be
      // distinct. This is possible, since the input graph is
      // triangulated. It'll be true all the time in a simple
      // graph, but loops and parallel edges cause some complications.
      if (prior_vertex == v || prior_vertex == u) {
        prior_edge_itr = ei;
        continue;
      }

      // Skip any self-loops
      if (u == v) {
        continue;
      }

      // Move next_edge_itr (and next_vertex) forwards
      // past any loops or parallel edges
      while (next_vertex == v || next_vertex == u) {
        next_edge_itr = (std::next(next_edge_itr) == ei_end ? ei_start : std::next(next_edge_itr));
        next_vertex = (source(*next_edge_itr, g) == u ? target(*next_edge_itr, g) : source(*next_edge_itr, g));
      }

      if (status[v] == planar_canonical_ordering_state::unprocessed) {
        status[v] = planar_canonical_ordering_state::one_neighbor_processed;
        processed_neighbor[v] = u;
      } else if (status[v] == planar_canonical_ordering_state::one_neighbor_processed) {
        vertex_t x = processed_neighbor[v];
        // are edges (v,u) and (v,x) adjacent in the planar
        // embedding? if so, set status[v] = 1. otherwise, set
        // status[v] = 2.

        if ((next_vertex == x && !(first_vertex == u && second_vertex == x)) ||
            (prior_vertex == x && !(first_vertex == x && second_vertex == u))) {
          status[v] = planar_canonical_ordering_state::ready_to_be_processed;
        } else {
          status[v] = planar_canonical_ordering_state::ready_to_be_processed + 1;
        }
      } else if (status[v] > planar_canonical_ordering_state::one_neighbor_processed) {
        // check the two edges before and after (v,u) in the planar
        // embedding, and update status[v] accordingly

        bool processed_before = false;
        if (status[prior_vertex] == planar_canonical_ordering_state::processed) {
          processed_before = true;
        }

        bool processed_after = false;
        if (status[next_vertex] == planar_canonical_ordering_state::processed) {
          processed_after = true;
        }

        if (!processed_before && !processed_after) {
          ++status[v];
        } else if (processed_before && processed_after) {
          --status[v];
        }
      }

      if (status[v] == planar_canonical_ordering_state::ready_to_be_processed) {
        ready_to_be_processed.push_back(v);
      }

      prior_edge_itr = ei;
    }

    status[u] = planar_canonical_ordering_state::processed;
    *ordering = u;
    ++ordering;
  }
}

template <concepts::VertexListGraph Graph, concepts::LvalueVertexPropertyMap<Graph> PlanarEmbedding,
          std::output_iterator<graph_vertex_descriptor_t<Graph>> OutputIterator>
requires concepts::AdjacencyGraph<Graph>
void planar_canonical_ordering(const Graph& g, PlanarEmbedding embedding, OutputIterator ordering) {
  planar_canonical_ordering(g, embedding, ordering, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_PLANAR_CANONICAL_ORDERING_H_
