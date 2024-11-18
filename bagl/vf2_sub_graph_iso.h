// Copyright (C) 2012 Flavio De Lorenzi (fdlorenzi@gmail.com)
// Copyright (C) 2013 Jakob Lykke Andersen, University of Southern Denmark (jlandersen@imada.sdu.dk)
// Revision History:
//   8 April 2013: Fixed a typo in vf2_print_callback. (Flavio De Lorenzi)
// Copyright 2024 Mikael Persson - Modernized to C++20

// The algorithm implemented here is derived from original ideas by
// Pasquale Foggia and colaborators. For further information see
// e.g. Cordella et al. 2001, 2004.

#ifndef BAGL_BAGL_VF2_SUB_GRAPH_ISO_H_
#define BAGL_BAGL_VF2_SUB_GRAPH_ISO_H_

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <ranges>
#include <tuple>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "bagl/filtered_graph.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/mcgregor_common_subgraphs.h"  // for always_equivalent

namespace bagl {

// Default print_callback
template <typename Graph1, typename Graph2>
struct vf2_print_callback {
  vf2_print_callback(const Graph1& graph1, const Graph2& graph2) : graph1_(graph1), graph2_(graph2) {}

  template <typename CorrespondenceMap1To2, typename CorrespondenceMap2To1>
  bool operator()(CorrespondenceMap1To2 f, CorrespondenceMap2To1 /*unused*/) const {
    // Print (sub)graph isomorphism map
    for (auto v : vertices(graph1_)) {
      std::cout << '(' << get(vertex_index, graph1_, v) << ", " << get(vertex_index, graph2_, get(f, v)) << ") ";
    }

    std::cout << std::endl;

    return true;
  }

 private:
  const Graph1& graph1_;
  const Graph2& graph2_;
};

namespace vf2_subgraph_detail {

// State associated with a single graph (graph_this)
template <typename GraphThis, typename GraphOther, concepts::ReadableVertexIndexMap<GraphThis> IndexMapThis,
          concepts::ReadableVertexIndexMap<GraphOther> IndexMapOther>
class base_state {
  using vertex_this_type = graph_vertex_descriptor_t<GraphThis>;
  using vertex_other_type = graph_vertex_descriptor_t<GraphOther>;

  const GraphThis& graph_this_;
  const GraphOther& graph_other_;

  IndexMapThis index_map_this_;
  IndexMapOther index_map_other_;

  vector_property_map<vertex_other_type, IndexMapThis> core_;

  vector_property_map<std::size_t, IndexMapThis> in_;
  vector_property_map<std::size_t, IndexMapThis> out_;

  std::size_t term_in_count_ = 0;
  std::size_t term_out_count_ = 0;
  std::size_t term_both_count_ = 0;
  std::size_t core_count_ = 0;

 public:
  // Forbidden
  base_state(const base_state&) = delete;
  base_state& operator=(const base_state&) = delete;
  base_state(base_state&&) = delete;
  base_state& operator=(base_state&&) = delete;
  ~base_state() = default;

  base_state(const GraphThis& graph_this, const GraphOther& graph_other, IndexMapThis index_map_this,
             IndexMapOther index_map_other)
      : graph_this_(graph_this),
        graph_other_(graph_other),
        index_map_this_(index_map_this),
        index_map_other_(index_map_other),
        core_(num_vertices(graph_this_), index_map_this_, graph_traits<GraphOther>::null_vertex()),
        in_(num_vertices(graph_this_), index_map_this_, 0),
        out_(num_vertices(graph_this_), index_map_this_, 0) {}

  // Adds a vertex pair to the state of graph graph_this
  void push(const vertex_this_type& v_this, const vertex_other_type& v_other) {
    ++core_count_;

    put(core_, v_this, v_other);

    if (!get(in_, v_this)) {
      put(in_, v_this, core_count_);
      ++term_in_count_;
      if (get(out_, v_this)) {
        ++term_both_count_;
      }
    }

    if (!get(out_, v_this)) {
      put(out_, v_this, core_count_);
      ++term_out_count_;
      if (get(in_, v_this)) {
        ++term_both_count_;
      }
    }

    for (auto e : in_edges(v_this, graph_this_)) {
      vertex_this_type w = source(e, graph_this_);
      if (!get(in_, w)) {
        put(in_, w, core_count_);
        ++term_in_count_;
        if (get(out_, w)) {
          ++term_both_count_;
        }
      }
    }

    for (auto e : out_edges(v_this, graph_this_)) {
      vertex_this_type w = target(e, graph_this_);
      if (!get(out_, w)) {
        put(out_, w, core_count_);
        ++term_out_count_;
        if (get(in_, w)) {
          ++term_both_count_;
        }
      }
    }
  }

  // Removes vertex pair from state of graph_this
  void pop(const vertex_this_type& v_this, const vertex_other_type& /*unused*/) {
    if (core_count_ == 0) {
      return;
    }

    if (get(in_, v_this) == core_count_) {
      put(in_, v_this, 0);
      --term_in_count_;
      if (get(out_, v_this)) {
        --term_both_count_;
      }
    }

    for (auto e : in_edges(v_this, graph_this_)) {
      vertex_this_type w = source(e, graph_this_);
      if (get(in_, w) == core_count_) {
        put(in_, w, 0);
        --term_in_count_;
        if (get(out_, w)) {
          --term_both_count_;
        }
      }
    }

    if (get(out_, v_this) == core_count_) {
      put(out_, v_this, 0);
      --term_out_count_;
      if (get(in_, v_this)) {
        --term_both_count_;
      }
    }

    for (auto e : out_edges(v_this, graph_this_)) {
      vertex_this_type w = target(e, graph_this_);
      if (get(out_, w) == core_count_) {
        put(out_, w, 0);
        --term_out_count_;
        if (get(in_, w)) {
          --term_both_count_;
        }
      }
    }
    put(core_, v_this, graph_traits<GraphOther>::null_vertex());

    --core_count_;
  }

  // Returns true if the in-terminal set is not empty
  [[nodiscard]] bool term_in() const { return core_count_ < term_in_count_; }

  // Returns true if vertex belongs to the in-terminal set
  [[nodiscard]] bool term_in(const vertex_this_type& v) {
    return (get(in_, v) > 0) && (get(core_, v) == graph_traits<GraphOther>::null_vertex());
  }

  // Returns true if the out-terminal set is not empty
  [[nodiscard]] bool term_out() const { return core_count_ < term_out_count_; }

  // Returns true if vertex belongs to the out-terminal set
  [[nodiscard]] bool term_out(const vertex_this_type& v) {
    return (get(out_, v) > 0) && (get(core_, v) == graph_traits<GraphOther>::null_vertex());
  }

  // Returns true of both (in- and out-terminal) sets are not empty
  [[nodiscard]] bool term_both() const { return core_count_ < term_both_count_; }

  // Returns true if vertex belongs to both (in- and out-terminal) sets
  [[nodiscard]] bool term_both(const vertex_this_type& v) {
    return (get(in_, v) > 0) && (get(out_, v) > 0) && (get(core_, v) == graph_traits<GraphOther>::null_vertex());
  }

  // Returns true if vertex belongs to the core map, i.e. it is in the
  // present mapping
  [[nodiscard]] bool in_core(const vertex_this_type& v) {
    return get(core_, v) != graph_traits<GraphOther>::null_vertex();
  }

  // Returns the number of vertices in the mapping
  [[nodiscard]] std::size_t count() const { return core_count_; }

  // Returns the image (in graph_other) of vertex v (in graph_this)
  [[nodiscard]] vertex_other_type core(const vertex_this_type& v) { return get(core_, v); }

  // Returns the mapping
  [[nodiscard]] auto get_map() { return core_.ref(); }

  // Returns the "time" (or depth) when vertex was added to the
  // in-terminal set
  [[nodiscard]] std::size_t in_depth(const vertex_this_type& v) { return get(in_, v); }

  // Returns the "time" (or depth) when vertex was added to the
  // out-terminal set
  [[nodiscard]] std::size_t out_depth(const vertex_this_type& v) { return get(out_, v); }

  // Returns the terminal set counts
  [[nodiscard]] std::tuple<std::size_t, std::size_t, std::size_t> term_set() const {
    return {term_in_count_, term_out_count_, term_both_count_};
  }
};

// Function object that checks whether a valid edge
// exists. For multi-graphs matched edges are excluded
template <typename Graph, typename Enable = void>
struct equivalent_edge_exists {
  template <typename EdgePredicate>
  bool operator()(graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> t, EdgePredicate is_valid_edge,
                  const Graph& g) {
    for (auto e : out_edges(s, g)) {
      if ((target(e, g) == t) && is_valid_edge(e) && (matched_edges_.find(e) == matched_edges_.end())) {
        matched_edges_.insert(e);
        return true;
      }
    }

    return false;
  }

 private:
  std::unordered_set<graph_edge_descriptor_t<Graph>> matched_edges_;
};

template <typename Graph>
struct equivalent_edge_exists<Graph, std::enable_if_t<!is_multigraph_v<Graph>>> {
  template <typename EdgePredicate>
  bool operator()(graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> t, EdgePredicate is_valid_edge,
                  const Graph& g) {
    auto [e, found] = edge(s, t, g);
    if (!found) {
      return false;
    }
    return is_valid_edge(e);
  }
};

enum class problem_selector { subgraph_mono, subgraph_iso, isomorphism };

// The actual state associated with both graphs
template <typename Graph1, typename Graph2, typename IndexMap1, typename IndexMap2, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate, typename SubGraphIsoMapCallback, problem_selector problem_selection>
class state {
  using vertex1_type = graph_vertex_descriptor_t<Graph1>;
  using vertex2_type = graph_vertex_descriptor_t<Graph2>;

  const Graph1& graph1_;
  const Graph2& graph2_;

  IndexMap1 index_map1_;

  EdgeEquivalencePredicate edge_comp_;
  VertexEquivalencePredicate vertex_comp_;

  base_state<Graph1, Graph2, IndexMap1, IndexMap2> state1_;
  base_state<Graph2, Graph1, IndexMap2, IndexMap1> state2_;

  // Three helper functions used in Feasibility and Valid functions to
  // test terminal set counts when testing for:
  [[nodiscard]] bool comp_term_sets(std::size_t a, std::size_t b) const {
    if constexpr (problem_selection == problem_selector::isomorphism) {
      return a == b;
    } else {
      return a <= b;
    }
  }

 public:
  // Forbidden
  state(const state&) = delete;
  state& operator=(const state&) = delete;
  state(state&&) = delete;
  state& operator=(state&&) = delete;
  ~state() = default;

  state(const Graph1& graph1, const Graph2& graph2, IndexMap1 index_map1, IndexMap2 index_map2,
        EdgeEquivalencePredicate edge_comp, VertexEquivalencePredicate vertex_comp)
      : graph1_(graph1),
        graph2_(graph2),
        index_map1_(index_map1),
        edge_comp_(edge_comp),
        vertex_comp_(vertex_comp),
        state1_(graph1, graph2, index_map1, index_map2),
        state2_(graph2, graph1, index_map2, index_map1) {}

  // Add vertex pair to the state
  void push(const vertex1_type& v, const vertex2_type& w) {
    state1_.push(v, w);
    state2_.push(w, v);
  }

  // Remove vertex pair from state
  void pop(const vertex1_type& v, const vertex2_type& /*unused*/) {
    vertex2_type w = state1_.core(v);
    state1_.pop(v, w);
    state2_.pop(w, v);
  }

  // Checks the feasibility of a new vertex pair
  bool feasible(const vertex1_type& v_new, const vertex2_type& w_new) {
    if (!vertex_comp_(v_new, w_new)) {
      return false;
    }

    // graph1
    std::size_t term_in1_count = 0;
    std::size_t term_out1_count = 0;
    std::size_t rest1_count = 0;

    {
      equivalent_edge_exists<Graph2> edge2_exists;

      for (auto e1 : in_edges(v_new, graph1_)) {
        vertex1_type v = source(e1, graph1_);

        if (state1_.in_core(v) || (v == v_new)) {
          vertex2_type w = w_new;
          if (v != v_new) {
            w = state1_.core(v);
          }
          if (!edge2_exists(
                  w, w_new, [this, e1](auto e2) { return edge_comp_(e1, e2); }, graph2_)) {
            return false;
          }
        } else {
          if (0 < state1_.in_depth(v)) {
            ++term_in1_count;
          }
          if (0 < state1_.out_depth(v)) {
            ++term_out1_count;
          }
          if ((state1_.in_depth(v) == 0) && (state1_.out_depth(v) == 0)) {
            ++rest1_count;
          }
        }
      }
    }

    {
      equivalent_edge_exists<Graph2> edge2_exists;

      for (auto e1 : out_edges(v_new, graph1_)) {
        vertex1_type v = target(e1, graph1_);
        if (state1_.in_core(v) || (v == v_new)) {
          vertex2_type w = w_new;
          if (v != v_new) {
            w = state1_.core(v);
          }

          if (!edge2_exists(
                  w_new, w, [this, e1](auto e2) { return edge_comp_(e1, e2); }, graph2_)) {
            return false;
          }
        } else {
          if (0 < state1_.in_depth(v)) {
            ++term_in1_count;
          }
          if (0 < state1_.out_depth(v)) {
            ++term_out1_count;
          }
          if ((state1_.in_depth(v) == 0) && (state1_.out_depth(v) == 0)) {
            ++rest1_count;
          }
        }
      }
    }

    // graph2
    std::size_t term_out2_count = 0;
    std::size_t term_in2_count = 0;
    std::size_t rest2_count = 0;

    {
      equivalent_edge_exists<Graph1> edge1_exists;

      for (auto e2 : in_edges(w_new, graph2_)) {
        vertex2_type w = source(e2, graph2_);
        if (state2_.in_core(w) || (w == w_new)) {
          if constexpr (problem_selection != problem_selector::subgraph_mono) {
            vertex1_type v = v_new;
            if (w != w_new) {
              v = state2_.core(w);
            }

            if (!edge1_exists(
                    v, v_new, [this, e2](auto e1) { return edge_comp_(e1, e2); }, graph1_)) {
              return false;
            }
          }
        } else {
          if (0 < state2_.in_depth(w)) {
            ++term_in2_count;
          }
          if (0 < state2_.out_depth(w)) {
            ++term_out2_count;
          }
          if ((state2_.in_depth(w) == 0) && (state2_.out_depth(w) == 0)) {
            ++rest2_count;
          }
        }
      }
    }

    {
      equivalent_edge_exists<Graph1> edge1_exists;

      for (auto e2 : out_edges(w_new, graph2_)) {
        vertex2_type w = target(e2, graph2_);
        if (state2_.in_core(w) || (w == w_new)) {
          if constexpr (problem_selection != problem_selector::subgraph_mono) {
            vertex1_type v = v_new;
            if (w != w_new) {
              v = state2_.core(w);
            }

            if (!edge1_exists(
                    v_new, v, [this, e2](auto e1) { return edge_comp_(e1, e2); }, graph1_)) {
              return false;
            }
          }
        } else {
          if (0 < state2_.in_depth(w)) {
            ++term_in2_count;
          }
          if (0 < state2_.out_depth(w)) {
            ++term_out2_count;
          }
          if ((state2_.in_depth(w) == 0) && (state2_.out_depth(w) == 0)) {
            ++rest2_count;
          }
        }
      }
    }

    if constexpr (problem_selection != problem_selector::subgraph_mono) {
      // subgraph_iso and isomorphism
      return comp_term_sets(term_in1_count, term_in2_count) && comp_term_sets(term_out1_count, term_out2_count) &&
             comp_term_sets(rest1_count, rest2_count);
    } else {
      // subgraph_mono
      return comp_term_sets(term_in1_count, term_in2_count) && comp_term_sets(term_out1_count, term_out2_count) &&
             comp_term_sets(term_in1_count + term_out1_count + rest1_count,
                            term_in2_count + term_out2_count + rest2_count);
    }
  }

  // Returns true if vertex v in graph1 is a possible candidate to
  // be added to the current state
  [[nodiscard]] bool possible_candidate1(const vertex1_type& v) {
    if (state1_.term_both() && state2_.term_both()) {
      return state1_.term_both(v);
    }
    if (state1_.term_out() && state2_.term_out()) {
      return state1_.term_out(v);
    }
    if (state1_.term_in() && state2_.term_in()) {
      return state1_.term_in(v);
    }
    return !state1_.in_core(v);
  }

  // Returns true if vertex w in graph2 is a possible candidate to
  // be added to the current state
  [[nodiscard]] bool possible_candidate2(const vertex2_type& w) {
    if (state1_.term_both() && state2_.term_both()) {
      return state2_.term_both(w);
    }
    if (state1_.term_out() && state2_.term_out()) {
      return state2_.term_out(w);
    }
    if (state1_.term_in() && state2_.term_in()) {
      return state2_.term_in(w);
    }
    return !state2_.in_core(w);
  }

  // Returns true if a mapping was found
  [[nodiscard]] bool success() const { return state1_.count() == num_vertices(graph1_); }

  // Returns true if a state is valid
  [[nodiscard]] bool valid() const {
    auto term1 = state1_.term_set();
    auto term2 = state2_.term_set();

    return comp_term_sets(std::get<0>(term1), std::get<0>(term2)) &&
           comp_term_sets(std::get<1>(term1), std::get<1>(term2)) &&
           comp_term_sets(std::get<2>(term1), std::get<2>(term2));
  }

  // Calls the user_callback with a graph (sub)graph mapping
  bool call_back(SubGraphIsoMapCallback user_callback) { return user_callback(state1_.get_map(), state2_.get_map()); }
};

// Non-recursive method that explores state space using a depth-first
// search strategy.  At each depth possible pairs candidate are compute
// and tested for feasibility to extend the mapping. If a complete
// mapping is found, the mapping is output to user_callback in the form
// of a correspondence map (graph1 to graph2). Returning false from the
// user_callback will terminate the search. Function match will return
// true if the entire search space was explored.
template <typename Graph1, typename Graph2, typename IndexMap1, typename IndexMap2, typename VertexOrder1,
          typename EdgeEquivalencePredicate, typename VertexEquivalencePredicate, typename SubGraphIsoMapCallback,
          problem_selector problem_selection>
bool match(const Graph1& graph1, const Graph2& graph2, SubGraphIsoMapCallback user_callback,
           const VertexOrder1& vertex_order1,
           state<Graph1, Graph2, IndexMap1, IndexMap2, EdgeEquivalencePredicate, VertexEquivalencePredicate,
                 SubGraphIsoMapCallback, problem_selection>& s) {
  // K&R C style variable declation due to the mess of "goto" below.

  auto graph1_verts_iter = vertex_order1.begin();

  auto graph2_verts_rg = vertices(graph2);
  auto graph2_verts_iter = graph2_verts_rg.begin();
  auto graph2_verts_iter_end = graph2_verts_rg.end();

  using match_continuation_type = std::pair<decltype(graph1_verts_iter), decltype(graph2_verts_iter)>;
  std::vector<match_continuation_type> k;
  bool found_match = false;

recur:
  if (s.success()) {
    if (!s.call_back(user_callback)) {
      return true;
    }
    found_match = true;

    goto back_track;
  }

  if (!s.valid()) {
    goto back_track;
  }

  graph1_verts_iter = vertex_order1.begin();
  while (graph1_verts_iter != vertex_order1.end() && !s.possible_candidate1(*graph1_verts_iter)) {
    ++graph1_verts_iter;
  }

  graph2_verts_iter = graph2_verts_rg.begin();
  while (graph2_verts_iter != graph2_verts_iter_end) {
    if (s.possible_candidate2(*graph2_verts_iter)) {
      if (s.feasible(*graph1_verts_iter, *graph2_verts_iter)) {
        k.emplace_back(graph1_verts_iter, graph2_verts_iter);
        s.push(*graph1_verts_iter, *graph2_verts_iter);
        goto recur;
      }
    }
  graph2_loop:
    ++graph2_verts_iter;
  }

back_track:
  if (k.empty()) {
    return found_match;
  }

  std::tie(graph1_verts_iter, graph2_verts_iter) = k.back();
  k.pop_back();
  s.pop(*graph1_verts_iter, *graph2_verts_iter);

  goto graph2_loop;
}

// Sorts vertices of a graph by multiplicity of in/out degrees
template <typename Graph, typename IndexMap, typename VertexOrder>
void sort_vertices(const Graph& graph, IndexMap index_map, VertexOrder& order) {
  std::ranges::sort(order, [&graph](const auto& v, const auto& w) {
    // lexicographical comparison
    return std::pair{in_degree(v, graph), out_degree(v, graph)} < std::pair{in_degree(w, graph), out_degree(w, graph)};
  });

  auto freq = vector_property_map(num_vertices(graph), index_map, std::size_t{0});

  for (auto order_iter = order.begin(); order_iter != order.end();) {
    std::size_t count = 0;
    for (auto count_iter = order_iter;
         (count_iter != order.end()) && (in_degree(*order_iter, graph) == in_degree(*count_iter, graph)) &&
         (out_degree(*order_iter, graph) == out_degree(*count_iter, graph));
         ++count_iter) {
      ++count;
    }

    for (std::size_t i = 0; i < count; ++i) {
      freq[*order_iter] = count;
      ++order_iter;
    }
  }

  std::ranges::sort(order, [&graph, &freq](const auto& v, const auto& w) {
    // lexicographical comparison
    return std::pair{freq[v], in_degree(v, graph) + out_degree(v, graph)} <
           std::pair{freq[w], in_degree(w, graph) + out_degree(w, graph)};
  });
}

// Enumerates all graph sub-graph mono-/iso-morphism mappings between graphs
// graph_small and graph_large. Continues until user_callback returns true
// or the search space has been fully explored.
template <
    problem_selector problem_selection, concepts::BidirectionalGraph GraphSmall,
    concepts::BidirectionalGraph GraphLarge, concepts::ReadableVertexIndexMap<GraphSmall> IndexMapSmall,
    concepts::ReadableVertexIndexMap<GraphLarge> IndexMapLarge,
    std::ranges::output_range<graph_vertex_descriptor_t<GraphSmall>> VertexOrderSmall,
    std::invocable<graph_edge_descriptor_t<GraphSmall>, graph_edge_descriptor_t<GraphLarge>> EdgeEquivalencePredicate,
    std::invocable<graph_vertex_descriptor_t<GraphSmall>, graph_vertex_descriptor_t<GraphLarge>>
        VertexEquivalencePredicate,
    typename SubGraphIsoMapCallback>
requires concepts::VertexAndEdgeListGraph<GraphSmall> && concepts::AdjacencyMatrix<GraphSmall> &&
    concepts::VertexAndEdgeListGraph<GraphLarge> && concepts::AdjacencyMatrix<GraphLarge>
bool vf2_subgraph_morphism(const GraphSmall& graph_small, const GraphLarge& graph_large,
                           SubGraphIsoMapCallback user_callback, IndexMapSmall index_map_small,
                           IndexMapLarge index_map_large, const VertexOrderSmall& vertex_order_small,
                           EdgeEquivalencePredicate edge_comp, VertexEquivalencePredicate vertex_comp) {
  // Vertex order requirements
  assert(num_vertices(graph_small) == vertex_order_small.size());

  if (num_vertices(graph_small) > num_vertices(graph_large)) {
    return false;
  }

  std::size_t num_edges_small = num_edges(graph_small);
  std::size_t num_edges_large = num_edges(graph_large);

  // Double the number of edges for undirected graphs: each edge counts as
  // in-edge and out-edge
  if constexpr (is_undirected_graph_v<GraphSmall>) {
    num_edges_small *= 2;
  }
  if constexpr (is_undirected_graph_v<GraphLarge>) {
    num_edges_large *= 2;
  }
  if (num_edges_small > num_edges_large) {
    return false;
  }

  state<GraphSmall, GraphLarge, IndexMapSmall, IndexMapLarge, EdgeEquivalencePredicate, VertexEquivalencePredicate,
        SubGraphIsoMapCallback, problem_selection>
      s(graph_small, graph_large, index_map_small, index_map_large, edge_comp, vertex_comp);

  return match(graph_small, graph_large, user_callback, vertex_order_small, s);
}

}  // namespace vf2_subgraph_detail

// Returns vertex order (vertices sorted by multiplicity of in/out degrees)
template <typename Graph>
auto vertex_order_by_mult(const Graph& graph) {
  std::vector<graph_vertex_descriptor_t<Graph>> vertex_order;
  std::ranges::copy(vertices(graph), std::back_inserter(vertex_order));

  vf2_subgraph_detail::sort_vertices(graph, get(vertex_index, graph), vertex_order);
  return vertex_order;
}

// Enumerates all graph sub-graph monomorphism mappings between graphs
// graph_small and graph_large. Continues until user_callback returns true or
// the search space has been fully explored.
template <typename GraphSmall, typename GraphLarge, typename IndexMapSmall, typename IndexMapLarge,
          typename VertexOrderSmall, typename EdgeEquivalencePredicate, typename VertexEquivalencePredicate,
          typename SubGraphIsoMapCallback>
bool vf2_subgraph_mono(const GraphSmall& graph_small, const GraphLarge& graph_large,
                       SubGraphIsoMapCallback user_callback, IndexMapSmall index_map_small,
                       IndexMapLarge index_map_large, const VertexOrderSmall& vertex_order_small,
                       EdgeEquivalencePredicate edge_comp, VertexEquivalencePredicate vertex_comp) {
  return vf2_subgraph_detail::vf2_subgraph_morphism<vf2_subgraph_detail::problem_selector::subgraph_mono>(
      graph_small, graph_large, user_callback, index_map_small, index_map_large, vertex_order_small, edge_comp,
      vertex_comp);
}

// All default interface for vf2_subgraph_iso
template <typename GraphSmall, typename GraphLarge, typename SubGraphIsoMapCallback>
bool vf2_subgraph_mono(const GraphSmall& graph_small, const GraphLarge& graph_large,
                       SubGraphIsoMapCallback user_callback) {
  return vf2_subgraph_mono(graph_small, graph_large, user_callback, get(vertex_index, graph_small),
                           get(vertex_index, graph_large), vertex_order_by_mult(graph_small), always_equivalent(),
                           always_equivalent());
}

// Enumerates all graph sub-graph isomorphism mappings between graphs
// graph_small and graph_large. Continues until user_callback returns true or
// the search space has been fully explored.
template <typename GraphSmall, typename GraphLarge, typename IndexMapSmall, typename IndexMapLarge,
          typename VertexOrderSmall, typename EdgeEquivalencePredicate, typename VertexEquivalencePredicate,
          typename SubGraphIsoMapCallback>
bool vf2_subgraph_iso(const GraphSmall& graph_small, const GraphLarge& graph_large,
                      SubGraphIsoMapCallback user_callback, IndexMapSmall index_map_small,
                      IndexMapLarge index_map_large, const VertexOrderSmall& vertex_order_small,
                      EdgeEquivalencePredicate edge_comp, VertexEquivalencePredicate vertex_comp) {
  return vf2_subgraph_detail::vf2_subgraph_morphism<vf2_subgraph_detail::problem_selector::subgraph_iso>(
      graph_small, graph_large, user_callback, index_map_small, index_map_large, vertex_order_small, edge_comp,
      vertex_comp);
}

// All default interface for vf2_subgraph_iso
template <typename GraphSmall, typename GraphLarge, typename SubGraphIsoMapCallback>
bool vf2_subgraph_iso(const GraphSmall& graph_small, const GraphLarge& graph_large,
                      SubGraphIsoMapCallback user_callback) {
  return vf2_subgraph_iso(graph_small, graph_large, user_callback, get(vertex_index, graph_small),
                          get(vertex_index, graph_large), vertex_order_by_mult(graph_small), always_equivalent(),
                          always_equivalent());
}

// Enumerates all isomorphism mappings between graphs graph1_ and graph2_.
// Continues until user_callback returns true or the search space has been
// fully explored.
template <
    concepts::BidirectionalGraph Graph1, concepts::BidirectionalGraph Graph2, typename GraphIsoMapCallback,
    concepts::ReadableVertexIndexMap<Graph1> IndexMap1, concepts::ReadableVertexIndexMap<Graph2> IndexMap2,
    std::ranges::output_range<graph_vertex_descriptor_t<Graph1>> VertexOrder1,
    std::invocable<graph_edge_descriptor_t<Graph1>, graph_edge_descriptor_t<Graph2>> EdgeEquivalencePredicate,
    std::invocable<graph_vertex_descriptor_t<Graph1>, graph_vertex_descriptor_t<Graph2>> VertexEquivalencePredicate>
requires concepts::VertexAndEdgeListGraph<Graph1> && concepts::AdjacencyMatrix<Graph1> &&
    concepts::VertexAndEdgeListGraph<Graph2> && concepts::AdjacencyMatrix<Graph2>
bool vf2_graph_iso(const Graph1& graph1, const Graph2& graph2, GraphIsoMapCallback user_callback, IndexMap1 index_map1,
                   IndexMap2 index_map2, const VertexOrder1& vertex_order1, EdgeEquivalencePredicate edge_comp,
                   VertexEquivalencePredicate vertex_comp) {
  assert(num_vertices(graph1) == vertex_order1.size());
  if (num_vertices(graph1) != num_vertices(graph2)) {
    return false;
  }

  std::size_t num_edges1 = num_edges(graph1);
  std::size_t num_edges2 = num_edges(graph2);

  // Double the number of edges for undirected graphs: each edge counts as
  // in-edge and out-edge
  if constexpr (is_undirected_graph_v<Graph1>) {
    num_edges1 *= 2;
  }
  if constexpr (is_undirected_graph_v<Graph2>) {
    num_edges2 *= 2;
  }
  if (num_edges1 != num_edges2) {
    return false;
  }

  vf2_subgraph_detail::state<Graph1, Graph2, IndexMap1, IndexMap2, EdgeEquivalencePredicate, VertexEquivalencePredicate,
                             GraphIsoMapCallback, vf2_subgraph_detail::problem_selector::isomorphism>
      s(graph1, graph2, index_map1, index_map2, edge_comp, vertex_comp);

  return vf2_subgraph_detail::match(graph1, graph2, user_callback, vertex_order1, s);
}

// All default interface for vf2_graph_iso
template <typename Graph1, typename Graph2, typename GraphIsoMapCallback>
bool vf2_graph_iso(const Graph1& graph1, const Graph2& graph2, GraphIsoMapCallback user_callback) {
  return vf2_graph_iso(graph1, graph2, user_callback, get(vertex_index, graph1), get(vertex_index, graph2),
                       vertex_order_by_mult(graph1), always_equivalent(), always_equivalent());
}

// Verifies a graph (sub)graph isomorphism map
template <concepts::EdgeListGraph Graph1, concepts::AdjacencyMatrix Graph2,
          concepts::ReadableVertexPropertyMap<Graph1> CorresponenceMap1To2, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate>
bool verify_vf2_subgraph_iso(const Graph1& graph1, const Graph2& graph2, const CorresponenceMap1To2 f,
                             EdgeEquivalencePredicate edge_comp, VertexEquivalencePredicate vertex_comp) {
  vf2_subgraph_detail::equivalent_edge_exists<Graph2> edge2_exists;

  for (auto e1 : edges(graph1)) {
    auto s1 = source(e1, graph1);
    auto t1 = target(e1, graph1);
    auto s2 = get(f, s1);
    auto t2 = get(f, t1);

    if (!vertex_comp(s1, s2) || !vertex_comp(t1, t2)) {
      return false;
    }

    if (!edge2_exists(
            s2, t2, [&edge_comp, e1](auto e2) { return edge_comp(e1, e2); }, graph2)) {
      return false;
    }
  }

  return true;
}

// Variant of verify_subgraph_iso with all default parameters
template <typename Graph1, typename Graph2, typename CorresponenceMap1To2>
bool verify_vf2_subgraph_iso(const Graph1& graph1, const Graph2& graph2, const CorresponenceMap1To2 f) {
  return verify_vf2_subgraph_iso(graph1, graph2, f, always_equivalent(), always_equivalent());
}

}  // namespace bagl

#endif  // BAGL_BAGL_VF2_SUB_GRAPH_ISO_H_
