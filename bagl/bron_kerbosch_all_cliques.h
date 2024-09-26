// (C) Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CLIQUE_H_
#define BAGL_BAGL_CLIQUE_H_

#include <algorithm>
#include <concepts>
#include <deque>
#include <vector>

#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/lookup_edge.h"

namespace bagl {

namespace concepts {
template <typename V, typename Clique, typename G>
concept CliqueVisitor = std::copy_constructible<V> && requires(const V& vis, const Clique& c, G& g) {
  vis.clique(c, g);
};
}  // namespace concepts

// The algorithm implemented in this paper is based on the so-called
// Algorithm 457, published as:
//
//     @article{362367,
//         author = {Coen Bron and Joep Kerbosch},
//         title = {Algorithm 457: finding all cliques of an undirected graph},
//         journal = {Communications of the ACM},
//         volume = {16},
//         number = {9},
//         year = {1973},
//         issn = {0001-0782},
//         pages = {575--577},
//         doi = {http://doi.acm.org/10.1145/362342.362367},
//             publisher = {ACM Press},
//             address = {New York, NY, USA},
//         }
//
// Sort of. This implementation is adapted from the 1st version of the
// algorithm and does not implement the candidate selection optimization
// described as published - it could, it just doesn't yet.
//
// The algorithm is given as proportional to (3.14)^(n/3) power. This is
// not the same as O(...), but based on time measures and approximation.
//
// Unfortunately, this implementation may be less efficient on non-
// AdjacencyMatrix modeled graphs due to the non-constant implementation
// of the edge(u,v,g) functions.
//
// TODO: It might be worthwhile to provide functionality for passing
// a connectivity matrix to improve the efficiency of those lookups
// when needed. This could simply be passed as a BooleanMatrix
// s.t. edge(u,v,B) returns true or false. This could easily be
// abstracted for adjacency matricies.
//
// The following paper is interesting for a number of reasons. First,
// it lists a number of other such algorithms and second, it describes
// a new algorithm (that does not appear to require the edge(u,v,g)
// function and appears fairly efficient. It is probably worth investigating.
//
//      @article{DBLP:journals/tcs/TomitaTT06,
//          author = {Etsuji Tomita and Akira Tanaka and Haruhisa Takahashi},
//          title = {The worst-case time complexity for generating all maximal
//          cliques and computational experiments}, journal = {Theor. Comput.
//          Sci.}, volume = {363}, number = {1}, year = {2006}, pages = {28-42}
//          ee = {https://doi.org/10.1016/j.tcs.2006.06.015}
//      }

// The default clique_visitor supplies an empty visitation function.
struct clique_visitor {
  template <typename VertexSet, typename G>
  void clique(const VertexSet& /*unused*/, G& /*unused*/) {}
};

// The max_clique_visitor records the size of the maximum clique (but not the clique itself).
struct max_clique_visitor {
  explicit max_clique_visitor(std::size_t& max) : maximum(max) {}

  template <typename Clique, typename Graph>
  inline void clique(const Clique& p, const Graph& /*g*/) {
    maximum = std::max(maximum, p.size());
  }
  std::size_t& maximum;
};

inline max_clique_visitor find_max_clique(std::size_t& max) { return max_clique_visitor(max); }

namespace clique_detail {
template <typename G>
bool is_connected_to_clique(const G& g, graph_vertex_descriptor_t<G> u, graph_vertex_descriptor_t<G> v) {
  if constexpr (is_undirected_graph_v<G>) {
    return lookup_edge(u, v, g).second;
  } else {
    // Note that this could alternate between using an || to determine
    // full connectivity. I believe that this should produce strongly
    // connected components. Note that using && instead of || will
    // change the results to a fully connected subgraph (i.e., symmetric
    // edges between all vertices s.t., if a->b, then b->a.
    return lookup_edge(u, v, g).second && lookup_edge(v, u, g).second;
  }
}

template <concepts::Graph G, typename Container>
void filter_unconnected_vertices(const G& g, graph_vertex_descriptor_t<G> v, const Container& in, Container& out) {
  for (auto u : in) {
    if (is_connected_to_clique(g, v, u)) {
      out.push_back(u);
    }
  }
}

template <concepts::Graph G,
          typename Clique,     // compsub type
          typename Container,  // candidates/not type
          concepts::CliqueVisitor<Clique, G> V>
void extend_clique(const G& g, Clique& clique, Container& cands, Container& nots, V vis, std::size_t min) {
  // Is there vertex in nots that is connected to all vertices
  // in the candidate set? If so, no clique can ever be found.
  // This could be broken out into a separate function.
  for (auto n : nots) {
    // n is connected to all cands
    if (std::all_of(cands.begin(), cands.end(), [&](auto c) { return lookup_edge(n, c, g).second; })) {
      return;
    }
  }

  // TODO: the original algorithm 457 describes an alternative
  // (albeit really complicated) mechanism for selecting candidates.
  // The given optimizaiton seeks to bring about the above
  // condition sooner (i.e., there is a vertex in the not set
  // that is connected to all candidates). unfortunately, the
  // method they give for doing this is fairly unclear.

  // basically, for every vertex in nots, we should know how many
  // vertices it is disconnected from in the candidate set. if
  // we fix some vertex in the not set, then we want to keep
  // choosing vertices that are not connected to that fixed vertex.
  // apparently, by selecting fix point with the minimum number
  // of disconnections (i.e., the maximum number of connections
  // within the candidate set), then the previous condition wil
  // be reached sooner.

  // there's some other stuff about using the number of disconnects
  // as a counter, but i'm jot really sure i followed it.

  // TODO: If we min-sized cliques to visit, then theoretically, we
  // should be able to stop recursing if the clique falls below that
  // size - maybe?

  // otherwise, iterate over candidates and and test
  // for maxmimal cliquiness.
  for (auto i = cands.begin(); i != cands.end();) {
    auto candidate = *i;

    // add the candidate to the clique (keeping the iterator!)
    // typename Clique::iterator ci = clique.insert(clique.end(),
    // candidate);
    clique.push_back(candidate);

    // remove it from the candidate set
    i = cands.erase(i);

    // build new candidate and not sets by removing all vertices
    // that are not connected to the current candidate vertex.
    // these actually invert the operation, adding them to the new
    // sets if the vertices are connected. its semantically the same.
    Container new_cands;
    Container new_nots;
    filter_unconnected_vertices(g, candidate, cands, new_cands);
    filter_unconnected_vertices(g, candidate, nots, new_nots);

    if (new_cands.empty() && new_nots.empty()) {
      // our current clique is maximal since there's nothing
      // that's connected that we haven't already visited. If
      // the clique is below our radar, then we won't visit it.
      if (clique.size() >= min) {
        vis.clique(clique, g);
      }
    } else {
      // recurse to explore the new candidates
      extend_clique(g, clique, new_cands, new_nots, vis, min);
    }

    // we're done with this vertex, so we need to move it
    // to the nots, and remove the candidate from the clique.
    nots.push_back(candidate);
    clique.pop_back();
  }
}
}  // namespace clique_detail

template <concepts::IncidenceGraph G, concepts::CliqueVisitor<std::vector<graph_vertex_descriptor_t<G>>, G> V>
requires concepts::VertexListGraph<G> && concepts::AdjacencyMatrix<G>
void bron_kerbosch_all_cliques(const G& g, V vis, std::size_t min) {
  using Vertex = graph_vertex_descriptor_t<G>;
  using VertexSet = std::vector<Vertex>;

  auto v_rg = vertices(g);
  VertexSet cands(v_rg.begin(), v_rg.end());  // start with all vertices as candidates
  VertexSet nots;                             // start with no vertices visited

  VertexSet clique;  // the first clique is an empty vertex set
  clique_detail::extend_clique(g, clique, cands, nots, vis, min);
}

// NOTE: By default the minimum number of vertices per clique is set at 2
// because singleton cliques aren't really very interesting.
template <typename G, typename Visitor>
void bron_kerbosch_all_cliques(const G& g, Visitor vis) {
  bron_kerbosch_all_cliques(g, vis, 2);
}

template <typename G>
std::size_t bron_kerbosch_clique_number(const G& g) {
  std::size_t ret = 0;
  bron_kerbosch_all_cliques(g, find_max_clique(ret));
  return ret;
}

}  // namespace bagl

#endif  // BAGL_BAGL_CLIQUE_H_
