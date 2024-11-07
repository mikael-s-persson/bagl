// Copyright (c) 2005 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MAXIMUM_CARDINALITY_MATCHING_H_
#define BAGL_BAGL_MAXIMUM_CARDINALITY_MATCHING_H_

#include <algorithm>
#include <cassert>
#include <deque>
#include <list>
#include <utility>
#include <vector>

#include "bagl/bfl_d_ary_tree.h"
#include "bagl/depth_first_search.h"
#include "bagl/disjoint_sets.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/property_map.h"
#include "bagl/visitors.h"

namespace bagl {
namespace max_cardinality_detail {
enum class vertex_state { even, odd, unreached };
}  // namespace max_cardinality_detail

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
std::size_t matching_size(const Graph& g, MateMap mate, VertexIndexMap vm) {
  std::size_t size_of_matching = 0;
  for (auto v : vertices(g)) {
    if (get(mate, v) != graph_traits<Graph>::null_vertex() && get(vm, v) < get(vm, get(mate, v))) {
      ++size_of_matching;
    }
  }
  return size_of_matching;
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap>
std::size_t matching_size(const Graph& g, MateMap mate) {
  return matching_size(g, mate, get(vertex_index, g));
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
bool is_a_matching(const Graph& g, MateMap mate, VertexIndexMap) {
  for (auto v : vertices(g)) {
    if (get(mate, v) != graph_traits<Graph>::null_vertex() && v != get(mate, get(mate, v))) {
      return false;
    }
  }
  return true;
}

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap>
bool is_a_matching(const Graph& g, MateMap mate) {
  return is_a_matching(g, mate, get(vertex_index, g));
}

//============== Maximum Cardinality Matching Functors ==============

struct no_augmenting_path_finder {
  bool augment_matching() { return false; }
  template <typename PropertyMap>
  void get_current_matching(PropertyMap /*unused*/) {}
};

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::IncidenceGraph<Graph>
class edmonds_augmenting_path_finder {
  // This implementation of Edmonds' matching algorithm closely
  // follows Tarjan's description of the algorithm in "Data
  // Structures and Network Algorithms."

 public:
  // generates the type of an iterator property map from vertices to type X
  template <typename X>
  using map_vertex_to = vector_property_map<X, VertexIndexMap>;

  using vertex_descriptor_t = graph_vertex_descriptor_t<Graph>;
  using vertex_pair_t = std::pair<vertex_descriptor_t, vertex_descriptor_t>;
  using edge_descriptor_t = graph_edge_descriptor_t<Graph>;
  using vertex_list_t = std::deque<vertex_descriptor_t>;
  using edge_list_t = std::vector<edge_descriptor_t>;
  using vertex_to_vertex_map_t = map_vertex_to<vertex_descriptor_t>;
  using vertex_to_state_map_t = map_vertex_to<max_cardinality_detail::vertex_state>;
  using vertex_to_vertex_pair_map_t = map_vertex_to<vertex_pair_t>;
  using vertex_to_size_map_t = map_vertex_to<std::size_t>;

  edmonds_augmenting_path_finder(const Graph& g, MateMap mate, VertexIndexMap vm)
      : g_(g),
        vm_(vm),
        n_vertices_(num_vertices(g)),

        mate_(n_vertices_, vm_),
        ancestor_of_v_(n_vertices_, vm_),
        ancestor_of_w_(n_vertices_, vm_),
        vertex_state_(n_vertices_, vm_),
        origin_(n_vertices_, vm_),
        pred_(n_vertices_, vm_),
        bridge_(n_vertices_, vm_),
        ds_parent_map_(n_vertices_, vm_),
        ds_rank_map_(n_vertices_, vm_),

        ds_(ds_rank_map_.ref(), ds_parent_map_.ref()) {
    for (auto v : vertices(g_)) {
      mate_[v] = get(mate, v);
    }
  }

  bool augment_matching() {
    // As an optimization, some of these values can be saved from one
    // iteration to the next instead of being re-initialized each
    // iteration, allowing for "lazy blossom expansion." This is not
    // currently implemented.

    std::size_t timestamp = 0;
    even_edges_.clear();

    for (auto u : vertices(g_)) {
      origin_[u] = u;
      pred_[u] = u;
      ancestor_of_v_[u] = 0;
      ancestor_of_w_[u] = 0;
      ds_.make_set(u);

      if (mate_[u] == graph_traits<Graph>::null_vertex()) {
        vertex_state_[u] = max_cardinality_detail::vertex_state::even;
        for (auto e : out_edges(u, g_)) {
          if (target(e, g_) != u) {
            even_edges_.push_back(e);
          }
        }
      } else {
        vertex_state_[u] = max_cardinality_detail::vertex_state::unreached;
      }
    }

    // end initializations

    auto w_free_ancestor = graph_traits<Graph>::null_vertex();
    auto w = w_free_ancestor;
    auto v_free_ancestor = graph_traits<Graph>::null_vertex();
    auto v = v_free_ancestor;
    bool found_alternating_path = false;

    while (!even_edges_.empty() && !found_alternating_path) {
      // since we push even edges onto the back of the list as
      // they're discovered, taking them off the back will search
      // for augmenting paths depth-first.
      auto current_edge = even_edges_.back();
      even_edges_.pop_back();

      v = source(current_edge, g_);
      w = target(current_edge, g_);

      auto v_prime = origin_[ds_.find_set(v)];
      auto w_prime = origin_[ds_.find_set(w)];

      // because of the way we put all of the edges on the queue,
      // v_prime should be labeled V_EVEN; the following is a
      // little paranoid but it could happen...
      if (vertex_state_[v_prime] != max_cardinality_detail::vertex_state::even) {
        std::swap(v_prime, w_prime);
        std::swap(v, w);
      }

      if (vertex_state_[w_prime] == max_cardinality_detail::vertex_state::unreached) {
        vertex_state_[w_prime] = max_cardinality_detail::vertex_state::odd;
        auto w_prime_mate = mate_[w_prime];
        vertex_state_[w_prime_mate] = max_cardinality_detail::vertex_state::even;
        for (auto e : out_edges(w_prime_mate, g_)) {
          if (target(e, g_) != w_prime_mate) {
            even_edges_.push_back(e);
          }
        }
        pred_[w_prime] = v;
      }

      else if (vertex_state_[w_prime] == max_cardinality_detail::vertex_state::even && w_prime != v_prime) {
        // w_prime == v_prime can happen below if we get an edge that has
        // been shrunk into a blossom
        auto w_up = w_prime;
        auto v_up = v_prime;
        auto nearest_common_ancestor = graph_traits<Graph>::null_vertex();
        w_free_ancestor = graph_traits<Graph>::null_vertex();
        v_free_ancestor = graph_traits<Graph>::null_vertex();

        // We now need to distinguish between the case that
        // w_prime and v_prime share an ancestor under the
        // "parent" relation, in which case we've found a
        // blossom and should shrink it, or the case that
        // w_prime and v_prime both have distinct ancestors that
        // are free, in which case we've found an alternating
        // path between those two ancestors.

        ++timestamp;

        while (nearest_common_ancestor == graph_traits<Graph>::null_vertex() &&
               (v_free_ancestor == graph_traits<Graph>::null_vertex() ||
                w_free_ancestor == graph_traits<Graph>::null_vertex())) {
          ancestor_of_w_[w_up] = timestamp;
          ancestor_of_v_[v_up] = timestamp;

          if (w_free_ancestor == graph_traits<Graph>::null_vertex()) {
            w_up = parent(w_up);
          }
          if (v_free_ancestor == graph_traits<Graph>::null_vertex()) {
            v_up = parent(v_up);
          }

          if (mate_[v_up] == graph_traits<Graph>::null_vertex()) {
            v_free_ancestor = v_up;
          }
          if (mate_[w_up] == graph_traits<Graph>::null_vertex()) {
            w_free_ancestor = w_up;
          }

          if (ancestor_of_w_[v_up] == timestamp) {
            nearest_common_ancestor = v_up;
          } else if (ancestor_of_v_[w_up] == timestamp) {
            nearest_common_ancestor = w_up;
          } else if (v_free_ancestor == w_free_ancestor && v_free_ancestor != graph_traits<Graph>::null_vertex()) {
            nearest_common_ancestor = v_up;
          }
        }

        if (nearest_common_ancestor == graph_traits<Graph>::null_vertex()) {
          found_alternating_path = true;  // to break out of the loop
        } else {
          // shrink the blossom
          link_and_set_bridges(w_prime, nearest_common_ancestor, std::pair(w, v));
          link_and_set_bridges(v_prime, nearest_common_ancestor, std::pair(v, w));
        }
      }
    }

    if (!found_alternating_path) {
      return false;
    }

    // retrieve the augmenting path and put it in aug_path
    reversed_retrieve_augmenting_path(v, v_free_ancestor);
    retrieve_augmenting_path(w, w_free_ancestor);

    // augment the matching along aug_path
    while (!aug_path_.empty()) {
      auto a = aug_path_.front();
      aug_path_.pop_front();
      auto b = aug_path_.front();
      aug_path_.pop_front();
      mate_[a] = b;
      mate_[b] = a;
    }

    return true;
  }

  template <typename PropertyMap>
  void get_current_matching(PropertyMap pm) {
    for (auto v : vertices(g_)) {
      put(pm, v, mate_[v]);
    }
  }

  template <typename PropertyMap>
  void get_vertex_state_map(PropertyMap pm) {
    for (auto v : vertices(g_)) {
      put(pm, v, vertex_state_[origin_[ds_.find_set(v)]]);
    }
  }

 private:
  vertex_descriptor_t parent(vertex_descriptor_t x) {
    if (vertex_state_[x] == max_cardinality_detail::vertex_state::even &&
        mate_[x] != graph_traits<Graph>::null_vertex()) {
      return mate_[x];
    }
    if (vertex_state_[x] == max_cardinality_detail::vertex_state::odd) {
      return origin_[ds_.find_set(pred_[x])];
    }
    return x;
  }

  void link_and_set_bridges(vertex_descriptor_t x, vertex_descriptor_t stop_vertex, vertex_pair_t the_bridge) {
    for (vertex_descriptor_t v = x; v != stop_vertex; v = parent(v)) {
      ds_.union_set(v, stop_vertex);
      origin_[ds_.find_set(stop_vertex)] = stop_vertex;

      if (vertex_state_[v] == max_cardinality_detail::vertex_state::odd) {
        bridge_[v] = the_bridge;
        for (auto e : out_edges(v, g_)) {
          if (target(e, g_) != v) {
            even_edges_.push_back(e);
          }
        }
      }
    }
  }

  // Since none of the STL containers support both constant-time
  // concatenation and reversal, the process of expanding an
  // augmenting path once we know one exists is a little more
  // complicated than it has to be. If we know the path is from v to
  // w, then the augmenting path is recursively defined as:
  //
  // path(v,w) = [v], if v = w
  //           = concat([v, mate[v]], path(pred[mate[v]], w),
  //                if v != w and vertex_state[v] == max_cardinality_detail::vertex_state::even
  //           = concat([v], reverse(path(x,mate[v])), path(y,w)),
  //                if v != w, vertex_state[v] == max_cardinality_detail::vertex_state::odd, and
  //                bridge[v] = (x,y)
  //
  // These next two mutually recursive functions implement this definition.

  void retrieve_augmenting_path(vertex_descriptor_t v, vertex_descriptor_t w) {
    if (v == w) {
      aug_path_.push_back(v);
    } else if (vertex_state_[v] == max_cardinality_detail::vertex_state::even) {
      aug_path_.push_back(v);
      aug_path_.push_back(mate_[v]);
      retrieve_augmenting_path(pred_[mate_[v]], w);
    } else  // vertex_state[v] == max_cardinality_detail::vertex_state::odd
    {
      aug_path_.push_back(v);
      reversed_retrieve_augmenting_path(bridge_[v].first, mate_[v]);
      retrieve_augmenting_path(bridge_[v].second, w);
    }
  }

  void reversed_retrieve_augmenting_path(vertex_descriptor_t v, vertex_descriptor_t w) {
    if (v == w) {
      aug_path_.push_back(v);
    } else if (vertex_state_[v] == max_cardinality_detail::vertex_state::even) {
      reversed_retrieve_augmenting_path(pred_[mate_[v]], w);
      aug_path_.push_back(mate_[v]);
      aug_path_.push_back(v);
    } else  // vertex_state[v] == max_cardinality_detail::vertex_state::odd
    {
      reversed_retrieve_augmenting_path(bridge_[v].second, w);
      retrieve_augmenting_path(bridge_[v].first, mate_[v]);
      aug_path_.push_back(v);
    }
  }

  // private data members

  const Graph& g_;
  VertexIndexMap vm_;
  std::size_t n_vertices_;

  // iterator property maps
  vertex_to_vertex_map_t mate_;
  vertex_to_size_map_t ancestor_of_v_;
  vertex_to_size_map_t ancestor_of_w_;
  vertex_to_state_map_t vertex_state_;
  vertex_to_vertex_map_t origin_;
  vertex_to_vertex_map_t pred_;
  vertex_to_vertex_pair_map_t bridge_;
  vertex_to_vertex_map_t ds_parent_map_;
  vertex_to_size_map_t ds_rank_map_;

  vertex_list_t aug_path_;
  edge_list_t even_edges_;
  disjoint_sets<property_map_ref<vertex_to_size_map_t>, property_map_ref<vertex_to_vertex_map_t>> ds_;
};

// ================ Initial Matching Functors ===================

template <concepts::VertexAndEdgeListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
void greedy_matching(const Graph& g, MateMap mate) {
  for (auto v : vertices(g)) {
    put(mate, v, graph_traits<Graph>::null_vertex());
  }

  for (auto e : edges(g)) {
    auto u = source(e, g);
    auto v = target(e, g);
    if (u != v && get(mate, u) == get(mate, v)) {
      // only way equality can hold is if
      //   mate[u] == mate[v] == null_vertex
      put(mate, u, v);
      put(mate, v, u);
    }
  }
}

// The "extra greedy matching" is formed by repeating the
// following procedure as many times as possible: Choose the
// unmatched vertex v of minimum non-zero degree.  Choose the
// neighbor w of v which is unmatched and has minimum degree over
// all of v's neighbors. Add (u,v) to the matching. Ties for
// either choice are broken arbitrarily. This procedure takes time
// O(m log n), where m is the number of edges in the graph and n
// is the number of vertices.
template <concepts::VertexAndEdgeListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
void extra_greedy_matching(const Graph& g, MateMap mate) {
  using vertex_descriptor_t = graph_vertex_descriptor_t<Graph>;
  using directed_edges_vector_t = std::vector<std::pair<vertex_descriptor_t, vertex_descriptor_t>>;

  directed_edges_vector_t edge_list;
  for (auto v : vertices(g)) {
    put(mate, v, graph_traits<Graph>::null_vertex());
  }

  for (auto e : edges(g)) {
    auto u = source(e, g);
    auto v = target(e, g);
    if (u == v) {
      continue;
    }
    edge_list.emplace_back(u, v);
    edge_list.emplace_back(v, u);
  }

  // sort the edges by the degree of the target, then (using a
  // stable sort) by degree of the source
  std::sort(edge_list.begin(), edge_list.end(),
            [&g](const auto& x, const auto& y) { return out_degree(x.second, g) < out_degree(y.second, g); });
  std::stable_sort(edge_list.begin(), edge_list.end(),
                   [&g](const auto& x, const auto& y) { return out_degree(x.first, g) < out_degree(y.first, g); });

  // construct the extra greedy matching
  for (auto [u, v] : edge_list) {
    if (get(mate, u) == get(mate, v)) {
      // only way equality can hold is if mate[itr->first] ==
      // mate[itr->second] == null_vertex
      put(mate, u, v);
      put(mate, v, u);
    }
  }
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
void empty_matching(const Graph& g, MateMap mate) {
  for (auto v : vertices(g)) {
    put(mate, v, graph_traits<Graph>::null_vertex());
  }
}

// =============== Matching Verifiers ================

namespace max_cardinality_detail {

class odd_components_counter : public dfs_visitor<>
// This depth-first search visitor will count the number of connected
// components with an odd number of vertices. It's used by
// maximum_matching_verifier.
{
 public:
  explicit odd_components_counter(std::size_t& c_count) : count_(c_count) {}

  template <class Vertex, class Graph>
  void start_vertex(Vertex /*unused*/, Graph& /*unused*/) {
    parity_ = false;
  }

  template <class Vertex, class Graph>
  void discover_vertex(Vertex /*unused*/, Graph& /*unused*/) {
    parity_ = !parity_;
    parity_ ? ++count_ : --count_;
  }

 private:
  std::size_t& count_;
  bool parity_ = false;
};

}  // namespace max_cardinality_detail

template <concepts::VertexListGraph Graph, concepts::ReadableVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
requires concepts::IncidenceGraph<Graph>
bool maximum_cardinality_matching_verify(const Graph& g, MateMap mate, VertexIndexMap vm) {
  // For any graph G, let o(G) be the number of connected
  // components in G of odd size. For a subset S of G's vertex set
  // V(G), let (G - S) represent the subgraph of G induced by
  // removing all vertices in S from G. Let M(G) be the size of the
  // maximum cardinality matching in G. Then the Tutte-Berge
  // formula guarantees that
  //
  //           2 * M(G) = min ( |V(G)| + |U| + o(G - U) )
  //
  // where the minimum is taken over all subsets U of
  // V(G). Edmonds' algorithm finds a set U that achieves the
  // minimum in the above formula, namely the vertices labeled
  //"ODD." This function runs one iteration of Edmonds' algorithm
  // to find U, then verifies that the size of the matching given
  // by mate satisfies the Tutte-Berge formula.

  // first, make sure it's a valid matching
  if (!is_a_matching(g, mate, vm)) {
    return false;
  }

  // We'll try to augment the matching once. This serves two
  // purposes: first, if we find some augmenting path, the matching
  // is obviously non-maximum. Second, running edmonds' algorithm
  // on a graph with no augmenting path will create the
  // Edmonds-Gallai decomposition that we need as a certificate of
  // maximality - we can get it by looking at the vertex_state map
  // that results.
  edmonds_augmenting_path_finder augmentor(g, mate, vm);
  if (augmentor.augment_matching()) {
    return false;
  }

  auto vertex_state = vector_property_map(num_vertices(g), vm, max_cardinality_detail::vertex_state::unreached);
  augmentor.get_vertex_state_map(vertex_state.ref());

  // count the number of graph::detail::V_ODD vertices
  std::size_t num_odd_vertices = 0;
  for (auto v : vertices(g)) {
    if (vertex_state[v] == max_cardinality_detail::vertex_state::odd) {
      ++num_odd_vertices;
    }
  }

  // count the number of connected components with odd cardinality
  // in the graph without graph::detail::V_ODD vertices
  auto non_odd_vertex = [&vertex_state](const auto& v) {
    return get(vertex_state, v) != max_cardinality_detail::vertex_state::odd;
  };
  filtered_graph fg(g, keep_all(), non_odd_vertex);

  std::size_t num_odd_components = 0;
  max_cardinality_detail::odd_components_counter occ(num_odd_components);
  depth_first_search(fg, occ, two_bit_color_map(num_vertices(g), vm).ref());

  return (2 * matching_size(g, mate, vm) == num_vertices(g) + num_odd_vertices - num_odd_components);
}

template <typename MateMap, typename AugmentingPathFinder>
void matching_loop(MateMap mate, AugmentingPathFinder augmentor) {
  bool not_maximum_yet = true;
  while (not_maximum_yet) {
    not_maximum_yet = augmentor.augment_matching();
  }
  augmentor.get_current_matching(mate);
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
bool checked_edmonds_maximum_cardinality_matching(const Graph& g, MateMap mate, VertexIndexMap vm) {
  extra_greedy_matching(g, mate);
  matching_loop(mate, edmonds_augmenting_path_finder(g, mate, vm));
  return maximum_cardinality_matching_verify(g, mate, vm);
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
bool checked_edmonds_maximum_cardinality_matching(const Graph& g, MateMap mate) {
  return checked_edmonds_maximum_cardinality_matching(g, mate, get(vertex_index, g));
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap,
          concepts::ReadableVertexIndexMap<Graph> VertexIndexMap>
void edmonds_maximum_cardinality_matching(const Graph& g, MateMap mate, VertexIndexMap vm) {
  extra_greedy_matching(g, mate);
  matching_loop(mate, edmonds_augmenting_path_finder(g, mate, vm));
}

template <concepts::VertexListGraph Graph, concepts::ReadWriteVertexPropertyMap<Graph> MateMap>
void edmonds_maximum_cardinality_matching(const Graph& g, MateMap mate) {
  edmonds_maximum_cardinality_matching(g, mate, get(vertex_index, g));
}

}  // namespace bagl

#endif  // BAGL_BAGL_MAXIMUM_CARDINALITY_MATCHING_H_
