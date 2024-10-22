// Copyright 2004-5 The Trustees of Indiana University.
// Copyright 2002 Brad King and Douglas Gregor
//  Authors: Douglas Gregor
//           Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_PAGE_RANK_H_
#define BAGL_BAGL_PAGE_RANK_H_

#include <type_traits>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/subgraph.h"

namespace bagl {

namespace page_rank_detail {

struct iterations_limit {
  explicit iterations_limit(std::size_t n) : n_(n) {}

  template <typename RankMap, typename Graph>
  bool operator()(const RankMap& /*unused*/, const Graph& /*unused*/) {
    return n_-- == 0;
  }

 private:
  std::size_t n_;
};

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> RankMap,
          concepts::ReadWriteVertexPropertyMap<Graph> RankMap2>
void page_rank_step(const Graph& g, RankMap from_rank, RankMap2 to_rank, property_traits_value_t<RankMap> damping) {
  using RankValue = property_traits_value_t<RankMap>;
  if constexpr (is_bidirectional_graph_v<Graph>) {
    for (auto v : vertices(g)) {
      auto rank = RankValue{0};
      for (auto e : in_edges(v, g)) {
        rank += get(from_rank, source(e, g)) / out_degree(source(e, g), g);
      }
      put(to_rank, v, (RankValue{1} - damping) + damping * rank);
    }
  } else {
    // Set new rank maps
    for (auto v : vertices(g)) {
      put(to_rank, v, rank_type(1 - damping));
    }

    for (auto u : vertices(g)) {
      auto u_rank_out = damping * get(from_rank, u) / out_degree(u, g);
      for (auto v : adjacent_vertices(u, g)) {
        put(to_rank, v, get(to_rank, v) + u_rank_out);
      }
    }
  }
}
}  // namespace page_rank_detail

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> RankMap, typename Done,
          concepts::ReadWriteVertexPropertyMap<Graph> RankMap2>
void page_rank(const Graph& g, RankMap rank_map, Done done, property_traits_value_t<RankMap> damping, std::size_t n,
               RankMap2 rank_map2) {
  using RankValue = property_traits_value_t<RankMap>;

  auto initial_rank = RankValue{RankValue{1} / n};
  for (auto v : vertices(g)) {
    put(rank_map, v, initial_rank);
  }

  bool to_map_2 = true;
  while ((to_map_2 && !done(rank_map, g)) || (!to_map_2 && !done(rank_map2, g))) {
    if (to_map_2) {
      page_rank_detail::page_rank_step(g, rank_map, rank_map2, damping);
    } else {
      page_rank_detail::page_rank_step(g, rank_map2, rank_map, damping);
    }
    to_map_2 = !to_map_2;
  }

  if (!to_map_2) {
    for (auto v : vertices(g)) {
      put(rank_map, v, get(rank_map2, v));
    }
  }
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> RankMap, typename Done>
void page_rank(const Graph& g, RankMap rank_map, Done done, property_traits_value_t<RankMap> damping, std::size_t n) {
  using RankValue = property_traits_value_t<RankMap>;
  page_rank(g, rank_map, done, damping, n,
            vector_property_map(num_vertices(g), get(vertex_index, g), RankValue{0}).ref());
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> RankMap, typename Done>
void page_rank(const Graph& g, RankMap rank_map, Done done, property_traits_value_t<RankMap> damping = 0.85) {
  page_rank(g, rank_map, done, damping, num_vertices(g));
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> RankMap>
void page_rank(const Graph& g, RankMap rank_map) {
  page_rank(g, rank_map, page_rank_detail::iterations_limit(20));
}

// TBD: this could be _much_ more efficient, using a queue to store
// the vertices that should be reprocessed and keeping track of which
// vertices are in the queue with a property map. Baah, this only
// applies when we have a bidirectional graph.
template <concepts::VertexListGraph G>
requires concepts::MutableGraph<G>
void remove_dangling_links(G& g) {
  bool removed_something = false;
  do {
    removed_something = false;

    auto v_rg = vertices(g);
    if constexpr (std::is_integral_v<graph_vertex_descriptor_t<G>> &&
                  std::ranges::random_access_range<decltype(v_rg)>) {
      // Probably a vector-like storage.
      // Proceed in reverse.
      for (auto v : vertices(g) | std::views::reverse) {
        if (out_degree(v, g) == 0) {
          clear_vertex(v, g);
          remove_vertex(v, g);
          removed_something = true;
          // next `v` is at a lower index than `v`, so we can continue.
        }
      }
    } else if constexpr (std::input_iterator<graph_vertex_descriptor_t<G>> &&
                         !std::ranges::random_access_range<decltype(v_rg)>) {
      // Probably a list/set/map-like storage.
      for (auto it = v_rg.begin(); it != v_rg.end(); /*increment before removal*/) {
        auto v = *it++;
        if (out_degree(v, g) == 0) {
          clear_vertex(v, g);
          remove_vertex(v, g);
          removed_something = true;
          // `it` already moved beyond `v`, so we can continue.
        }
      }
    } else {
      // Something else, like a filtered range.
      // Proceed in reverse in case underlying storage is contiguous.
      for (auto v : vertices(g) | std::views::reverse) {
        if (out_degree(v, g) == 0) {
          clear_vertex(v, g);
          remove_vertex(v, g);
          removed_something = true;
          // We cannot continue to iterate after vertex removal.
          break;
        }
      }
    }
  } while (removed_something);
}

}  // namespace bagl

#endif  // BOOST_GRAPH_PAGE_RANK_HPP
