// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen, Andrew Lumsdaine
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_MCGREGOR_COMMON_SUBGRAPHS_H_
#define BAGL_BAGL_MCGREGOR_COMMON_SUBGRAPHS_H_

#include <algorithm>
#include <memory>
#include <stack>
#include <type_traits>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"

namespace bagl {

// ==========================================================================

// Binary function object that returns true if the values for item1
// in property_map1 and item2 in property_map2 are equivalent.
template <typename PropertyMapFirst, typename PropertyMapSecond>
struct property_map_equivalent {
  property_map_equivalent(const PropertyMapFirst property_map1, const PropertyMapSecond property_map2)
      : property_map1_(property_map1), property_map2_(property_map2) {}

  template <typename ItemFirst, typename ItemSecond>
  bool operator()(const ItemFirst item1, const ItemSecond item2) {
    return (get(property_map1_, item1) == get(property_map2_, item2));
  }

 private:
  const PropertyMapFirst property_map1_;
  const PropertyMapSecond property_map2_;
};

// Binary function object that always returns true.  Used when
// vertices or edges are always equivalent (i.e. have no labels).
struct always_equivalent {
  template <typename ItemFirst, typename ItemSecond>
  bool operator()(const ItemFirst& /*unused*/, const ItemSecond& /*unused*/) {
    return true;
  }
};

// ==========================================================================

namespace mcgregor_detail {

template <concepts::AdjacencyMatrix Graph>
auto find_edge(graph_vertex_descriptor_t<Graph> u, graph_vertex_descriptor_t<Graph> v, const Graph& g) {
  return edge(u, v, g);
}
// Fallback to incidence graph.
template <typename Graph>
std::pair<graph_edge_descriptor_t<Graph>, bool> find_edge(graph_vertex_descriptor_t<Graph> u,
                                                          graph_vertex_descriptor_t<Graph> v, const Graph& g) {
  for (auto e : out_edges(u, g)) {
    if (target(e, g) == v) {
      return {e, true};
    }
  }
  return {{}, false};
}

// Return true if new_vertex1 and new_vertex2 can extend the
// subgraph represented by correspondence_map_1_to_2 and
// correspondence_map_2_to_1.  The vertices_equivalent and
// edges_equivalent predicates are used to test vertex and edge
// equivalency between the two graphs.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexPropertyMap<GraphFirst> CorrespondenceMapFirstToSecond,
          concepts::ReadableVertexPropertyMap<GraphSecond> CorrespondenceMapSecondToFirst,
          typename EdgeEquivalencePredicate, typename VertexEquivalencePredicate>
bool can_extend_graph(const GraphFirst& graph1, const GraphSecond& graph2,
                      CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                      CorrespondenceMapSecondToFirst /*correspondence_map_2_to_1*/, std::size_t subgraph_size,
                      graph_vertex_descriptor_t<GraphFirst> new_vertex1,
                      graph_vertex_descriptor_t<GraphSecond> new_vertex2, EdgeEquivalencePredicate edges_equivalent,
                      VertexEquivalencePredicate vertices_equivalent, bool only_connected_subgraphs) {
  // Check vertex equality
  if (!vertices_equivalent(new_vertex1, new_vertex2)) {
    return (false);
  }

  // Vertices match and graph is empty, so we can extend the subgraph
  if (subgraph_size == 0) {
    return (true);
  }

  bool has_one_edge = false;

  // Verify edges with existing sub-graph
  for (auto existing_vertex1 : vertices(graph1)) {
    auto existing_vertex2 = get(correspondence_map_1_to_2, existing_vertex1);

    // Skip unassociated vertices
    if (existing_vertex2 == graph_traits<GraphSecond>::null_vertex()) {
      continue;
    }

    // NOTE: This will not work with parallel edges, since the
    // first matching edge is always chosen.

    // Search for edge from existing to new vertex (graph1)
    auto [edge_to_new1, edge_to_new_exists1] = find_edge(existing_vertex1, new_vertex1, graph1);

    // Search for edge from existing to new vertex (graph2)
    auto [edge_to_new2, edge_to_new_exists2] = find_edge(existing_vertex2, new_vertex2, graph2);

    // Make sure edges from existing to new vertices are equivalent
    if ((edge_to_new_exists1 != edge_to_new_exists2) ||
        ((edge_to_new_exists1 && edge_to_new_exists2) && !edges_equivalent(edge_to_new1, edge_to_new2))) {
      return false;
    }

    constexpr bool is_undirected1 = is_undirected_graph_v<GraphFirst>;
    constexpr bool is_undirected2 = is_undirected_graph_v<GraphSecond>;

    if constexpr (is_undirected1 && is_undirected2) {
      // Edge in both graphs exists and both graphs are undirected
      if (edge_to_new_exists1 && edge_to_new_exists2) {
        has_one_edge = true;
      }

      continue;
    } else {
      // Search for edge from new to existing vertex (graph1)
      auto [edge_from_new1, edge_from_new_exists1] =
          (!is_undirected1 ? find_edge(new_vertex1, existing_vertex1, graph1)
                           : std::pair{graph_edge_descriptor_t<GraphFirst>{}, false});

      // Search for edge from new to existing vertex (graph2)
      auto [edge_from_new2, edge_from_new_exists2] =
          (!is_undirected2 ? find_edge(new_vertex2, existing_vertex2, graph2)
                           : std::pair{graph_edge_descriptor_t<GraphSecond>{}, false});

      // Make sure edges from new to existing vertices are equivalent
      if ((edge_from_new_exists1 != edge_from_new_exists2) ||
          ((edge_from_new_exists1 && edge_from_new_exists2) && !edges_equivalent(edge_from_new1, edge_from_new2))) {
        return false;
      }

      if ((edge_from_new_exists1 && edge_from_new_exists2) || (edge_to_new_exists1 && edge_to_new_exists2)) {
        has_one_edge = true;
      }
    }
  }

  // Make sure new vertices are connected to the existing subgraph
  return !only_connected_subgraphs || has_one_edge;
}

// Recursive method that does a depth-first search in the space of
// potential subgraphs.  At each level, every new vertex pair from
// both graphs is tested to see if it can extend the current
// subgraph.  If so, the subgraph is output to subgraph_callback
// in the form of two correspondence maps (one for each graph).
// Returning false from subgraph_callback will terminate the
// search.  Function returns true if the entire search space was
// explored.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond,
          concepts::ReadableVertexPropertyMap<GraphFirst> CorrespondenceMapFirstToSecond,
          concepts::ReadableVertexPropertyMap<GraphSecond> CorrespondenceMapSecondToFirst, typename VertexStackFirst,
          typename EdgeEquivalencePredicate, typename VertexEquivalencePredicate, typename SubGraphInternalCallback>
bool mcgregor_common_subgraphs_internal(const GraphFirst& graph1, const GraphSecond& graph2,
                                        const VertexIndexMapFirst& vindex_map1, const VertexIndexMapSecond& vindex_map2,
                                        CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                                        CorrespondenceMapSecondToFirst correspondence_map_2_to_1,
                                        VertexStackFirst& vertex_stack1, EdgeEquivalencePredicate edges_equivalent,
                                        VertexEquivalencePredicate vertices_equivalent, bool only_connected_subgraphs,
                                        SubGraphInternalCallback subgraph_callback) {
  // Iterate until all vertices have been visited
  for (auto new_vertex1 : vertices(graph1)) {
    auto existing_vertex2 = get(correspondence_map_1_to_2, new_vertex1);

    // Skip already matched vertices in first graph
    if (existing_vertex2 != graph_traits<GraphSecond>::null_vertex()) {
      continue;
    }

    for (auto new_vertex2 : vertices(graph2)) {
      auto existing_vertex1 = get(correspondence_map_2_to_1, new_vertex2);

      // Skip already matched vertices in second graph
      if (existing_vertex1 != graph_traits<GraphFirst>::null_vertex()) {
        continue;
      }

      // Check if current sub-graph can be extended with the matched
      // vertex pair
      if (!can_extend_graph(graph1, graph2, correspondence_map_1_to_2, correspondence_map_2_to_1, vertex_stack1.size(),
                            new_vertex1, new_vertex2, edges_equivalent, vertices_equivalent,
                            only_connected_subgraphs)) {
        continue;
      }

      // Keep track of old graph size for restoring later
      auto old_graph_size = vertex_stack1.size();
      auto new_graph_size = old_graph_size + 1;

      // Extend subgraph
      put(correspondence_map_1_to_2, new_vertex1, new_vertex2);
      put(correspondence_map_2_to_1, new_vertex2, new_vertex1);
      vertex_stack1.push(new_vertex1);

      // Returning false from the callback will cancel iteration
      if (!subgraph_callback(correspondence_map_1_to_2, correspondence_map_2_to_1, new_graph_size)) {
        return false;
      }

      // Depth-first search into the state space of possible
      // sub-graphs
      bool continue_iteration = mcgregor_common_subgraphs_internal(
          graph1, graph2, vindex_map1, vindex_map2, correspondence_map_1_to_2, correspondence_map_2_to_1, vertex_stack1,
          edges_equivalent, vertices_equivalent, only_connected_subgraphs, subgraph_callback);

      if (!continue_iteration) {
        return false;
      }

      // Restore previous state
      if (vertex_stack1.size() > old_graph_size) {
        auto stack_vertex1 = vertex_stack1.top();
        auto stack_vertex2 = get(correspondence_map_1_to_2, stack_vertex1);

        // Contract subgraph
        put(correspondence_map_1_to_2, stack_vertex1, graph_traits<GraphSecond>::null_vertex());

        put(correspondence_map_2_to_1, stack_vertex2, graph_traits<GraphFirst>::null_vertex());

        vertex_stack1.pop();
      }
    }
  }

  return true;
}

// Internal method that initializes blank correspondence maps and
// a vertex stack for use in mcgregor_common_subgraphs_internal.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate, typename SubGraphInternalCallback>
void mcgregor_common_subgraphs_internal_init(const GraphFirst& graph1, const GraphSecond& graph2,
                                             const VertexIndexMapFirst vindex_map1,
                                             const VertexIndexMapSecond vindex_map2,
                                             EdgeEquivalencePredicate edges_equivalent,
                                             VertexEquivalencePredicate vertices_equivalent,
                                             bool only_connected_subgraphs,
                                             SubGraphInternalCallback subgraph_callback) {
  auto correspondence_map_1_to_2 =
      vector_property_map(num_vertices(graph1), vindex_map1, graph_traits<GraphSecond>::null_vertex());
  for (auto vertex1 : vertices(graph1)) {
    put(correspondence_map_1_to_2, vertex1, graph_traits<GraphSecond>::null_vertex());
  }

  auto correspondence_map_2_to_1 =
      vector_property_map(num_vertices(graph2), vindex_map2, graph_traits<GraphFirst>::null_vertex());
  for (auto vertex2 : vertices(graph2)) {
    put(correspondence_map_2_to_1, vertex2, graph_traits<GraphFirst>::null_vertex());
  }

  using VertexFirst = graph_vertex_descriptor_t<GraphFirst>;

  std::stack<VertexFirst> vertex_stack1;

  mcgregor_common_subgraphs_internal(graph1, graph2, vindex_map1, vindex_map2, correspondence_map_1_to_2,
                                     correspondence_map_2_to_1, vertex_stack1, edges_equivalent, vertices_equivalent,
                                     only_connected_subgraphs, subgraph_callback);
}

}  // namespace mcgregor_detail

// ==========================================================================

// Enumerates all common subgraphs present in graph1 and graph2.
// Continues until the search space has been fully explored or false
// is returned from user_callback.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate, typename SubGraphCallback>
void mcgregor_common_subgraphs(const GraphFirst& graph1, const GraphSecond& graph2,
                               const VertexIndexMapFirst vindex_map1,           // get(vertex_index, graph1)
                               const VertexIndexMapSecond vindex_map2,          // get(vertex_index, graph1)
                               EdgeEquivalencePredicate edges_equivalent,       // always_equivalent()
                               VertexEquivalencePredicate vertices_equivalent,  // always_equivalent()
                               bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_detail::mcgregor_common_subgraphs_internal_init(graph1, graph2, vindex_map1, vindex_map2, edges_equivalent,
                                                           vertices_equivalent, only_connected_subgraphs,
                                                           user_callback);
}

// Variant of mcgregor_common_subgraphs with all default parameters
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond, typename SubGraphCallback>
void mcgregor_common_subgraphs(const GraphFirst& graph1, const GraphSecond& graph2, bool only_connected_subgraphs,
                               SubGraphCallback user_callback) {
  mcgregor_detail::mcgregor_common_subgraphs_internal_init(
      graph1, graph2, get(vertex_index, graph1), get(vertex_index, graph2), always_equivalent(), always_equivalent(),
      only_connected_subgraphs, user_callback);
}

// ==========================================================================

namespace mcgregor_detail {

template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond,
          concepts::ReadableVertexPropertyMap<GraphFirst> CorrespondenceMapFirstToSecond,
          concepts::ReadableVertexPropertyMap<GraphSecond> CorrespondenceMapSecondToFirst>
auto copy_subgraph_maps(std::size_t subgraph_size, const GraphFirst& graph1, const GraphSecond& graph2,
                        VertexIndexMapFirst vindex_map1, VertexIndexMapSecond vindex_map2,
                        CorrespondenceMapFirstToSecond corr_1_to_2, CorrespondenceMapSecondToFirst corr_2_to_1) {
  // Subgraph is unique, so make a cached copy
  auto new_subgraph_1_to_2 =
      vector_property_map(num_vertices(graph1), vindex_map1, graph_traits<GraphSecond>::null_vertex());
  for (auto vertex1 : vertices(graph1)) {
    put(new_subgraph_1_to_2, vertex1, get(corr_1_to_2, vertex1));
  }

  auto new_subgraph_2_to_1 =
      vector_property_map(num_vertices(graph2), vindex_map2, graph_traits<GraphFirst>::null_vertex());
  for (auto vertex2 : vertices(graph2)) {
    put(new_subgraph_2_to_1, vertex2, get(corr_2_to_1, vertex2));
  }

  return std::tuple{subgraph_size, new_subgraph_1_to_2, new_subgraph_2_to_1};
}

// Binary function object that intercepts subgraphs from
// mcgregor_common_subgraphs_internal and maintains a cache of
// unique subgraphs.  The user callback is invoked for each unique
// subgraph.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename SubGraphCallback>
struct unique_subgraph_interceptor {
  using CachedCorrespondenceMapFirstToSecond =
      vector_property_map<graph_vertex_descriptor_t<GraphSecond>, VertexIndexMapFirst>;
  using CachedCorrespondenceMapSecondToFirst =
      vector_property_map<graph_vertex_descriptor_t<GraphFirst>, VertexIndexMapSecond>;

  using SubGraph = std::tuple<std::size_t, CachedCorrespondenceMapFirstToSecond, CachedCorrespondenceMapSecondToFirst>;
  using SubGraphList = std::vector<SubGraph>;

  unique_subgraph_interceptor(const GraphFirst& graph1, const GraphSecond& graph2,
                              const VertexIndexMapFirst vindex_map1, const VertexIndexMapSecond vindex_map2,
                              SubGraphCallback user_callback)
      : graph1_(graph1),
        graph2_(graph2),
        vindex_map1_(vindex_map1),
        vindex_map2_(vindex_map2),
        subgraphs_(std::make_shared<SubGraphList>()),
        user_callback_(user_callback) {}

  template <concepts::ReadableVertexPropertyMap<GraphFirst> CorrespondenceMapFirstToSecond,
            concepts::ReadableVertexPropertyMap<GraphSecond> CorrespondenceMapSecondToFirst>
  bool operator()(CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                  CorrespondenceMapSecondToFirst correspondence_map_2_to_1, std::size_t subgraph_size) {
    for (auto [cached_size, cached_1_to_2, cached_2_to_1] : *subgraphs_) {
      // Compare subgraph sizes
      if (subgraph_size != cached_size) {
        continue;
      }

      if (!are_property_maps_different(correspondence_map_1_to_2, cached_1_to_2, graph1_)) {
        // New subgraph is a duplicate
        return true;
      }
    }

    // Subgraph is unique, so make a cached copy
    subgraphs_->emplace_back(copy_subgraph_maps(subgraph_size, graph1_, graph2_, vindex_map1_, vindex_map2_,
                                                correspondence_map_1_to_2, correspondence_map_2_to_1));

    return user_callback_(correspondence_map_1_to_2, correspondence_map_2_to_1, subgraph_size);
  }

 protected:
  const GraphFirst& graph1_;
  const GraphFirst& graph2_;
  const VertexIndexMapFirst vindex_map1_;
  const VertexIndexMapSecond vindex_map2_;
  std::shared_ptr<SubGraphList> subgraphs_;
  SubGraphCallback user_callback_;
};

}  // namespace mcgregor_detail

// Enumerates all unique common subgraphs between graph1 and graph2.
// The user callback is invoked for each unique subgraph as they are
// discovered.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate, typename SubGraphCallback>
void mcgregor_common_subgraphs_unique(const GraphFirst& graph1, const GraphSecond& graph2,
                                      const VertexIndexMapFirst vindex_map1,           // get(vertex_index, graph1)
                                      const VertexIndexMapSecond vindex_map2,          // get(vertex_index, graph1)
                                      EdgeEquivalencePredicate edges_equivalent,       // always_equivalent()
                                      VertexEquivalencePredicate vertices_equivalent,  // always_equivalent()
                                      bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_detail::unique_subgraph_interceptor unique_callback(graph1, graph2, vindex_map1, vindex_map2, user_callback);

  mcgregor_detail::mcgregor_common_subgraphs_internal_init(graph1, graph2, vindex_map1, vindex_map2, edges_equivalent,
                                                           vertices_equivalent, only_connected_subgraphs,
                                                           unique_callback);
}

// Variant of mcgregor_common_subgraphs_unique with all default
// parameters.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond, typename SubGraphCallback>
void mcgregor_common_subgraphs_unique(const GraphFirst& graph1, const GraphSecond& graph2,
                                      bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_common_subgraphs_unique(graph1, graph2, get(vertex_index, graph1), get(vertex_index, graph2),
                                   always_equivalent(), always_equivalent(), only_connected_subgraphs, user_callback);
}

// ==========================================================================

namespace mcgregor_detail {

// Binary function object that intercepts subgraphs from
// mcgregor_common_subgraphs_internal and maintains a cache of the
// largest subgraphs.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename SubGraphCallback,
          typename CheckUniqueness>
struct maximum_subgraph_interceptor : unique_subgraph_interceptor<GraphFirst, GraphSecond, VertexIndexMapFirst,
                                                                  VertexIndexMapSecond, SubGraphCallback> {
  using Base =
      unique_subgraph_interceptor<GraphFirst, GraphSecond, VertexIndexMapFirst, VertexIndexMapSecond, SubGraphCallback>;

  maximum_subgraph_interceptor(const GraphFirst& graph1, const GraphSecond& graph2,
                               const VertexIndexMapFirst vindex_map1, const VertexIndexMapSecond vindex_map2,
                               SubGraphCallback user_callback, CheckUniqueness /*unused*/)
      : Base(graph1, graph2, vindex_map1, vindex_map2, user_callback) {}

  template <concepts::ReadableVertexPropertyMap<GraphFirst> CorrespondenceMapFirstToSecond,
            concepts::ReadableVertexPropertyMap<GraphSecond> CorrespondenceMapSecondToFirst>
  bool operator()(CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                  CorrespondenceMapSecondToFirst correspondence_map_2_to_1, std::size_t subgraph_size) {
    if (subgraph_size > largest_size_so_far_) {
      this->subgraphs_->clear();
      largest_size_so_far_ = subgraph_size;
    }

    if (subgraph_size == largest_size_so_far_) {
      if constexpr (CheckUniqueness::value) {
        // Check if subgraph is unique
        for (auto [cached_size, cached_1_to_2, cached_2_to_1] : *this->subgraphs_) {
          if (!are_property_maps_different(correspondence_map_1_to_2, cached_1_to_2, this->graph1_)) {
            // New subgraph is a duplicate
            return true;
          }
        }
      }

      // Make a cached copy
      this->subgraphs_->emplace_back(copy_subgraph_maps(subgraph_size, this->graph1_, this->graph2_, this->vindex_map1_,
                                                        this->vindex_map2_, correspondence_map_1_to_2,
                                                        correspondence_map_2_to_1));
    }

    return true;
  }

  void output_subgraphs() {
    for (auto [cached_size, cached_1_to_2, cached_2_to_1] : *this->subgraphs_) {
      user_callback_(cached_1_to_2, cached_2_to_1, cached_size);
    }
  }

 private:
  std::size_t largest_size_so_far_ = 0;
};

}  // namespace mcgregor_detail

// Enumerates the largest common subgraphs found between graph1
// and graph2.  Note that the ENTIRE search space is explored before
// user_callback is actually invoked.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate, typename SubGraphCallback>
void mcgregor_common_subgraphs_maximum(const GraphFirst& graph1, const GraphSecond& graph2,
                                       const VertexIndexMapFirst vindex_map1,           // get(vertex_index, graph1)
                                       const VertexIndexMapSecond vindex_map2,          // get(vertex_index, graph1)
                                       EdgeEquivalencePredicate edges_equivalent,       // always_equivalent()
                                       VertexEquivalencePredicate vertices_equivalent,  // always_equivalent()
                                       bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_detail::maximum_subgraph_interceptor max_interceptor(graph1, graph2, vindex_map1, vindex_map2, user_callback,
                                                                /*CheckUniqueness=*/std::false_type{});

  mcgregor_detail::mcgregor_common_subgraphs_internal_init(graph1, graph2, vindex_map1, vindex_map2, edges_equivalent,
                                                           vertices_equivalent, only_connected_subgraphs,
                                                           max_interceptor);

  // Only output the largest subgraphs
  max_interceptor.output_subgraphs();
}

// Variant of mcgregor_common_subgraphs_maximum with all default
// parameters.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond, typename SubGraphCallback>
void mcgregor_common_subgraphs_maximum(const GraphFirst& graph1, const GraphSecond& graph2,
                                       bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_common_subgraphs_maximum(graph1, graph2, get(vertex_index, graph1), get(vertex_index, graph2),
                                    always_equivalent(), always_equivalent(), only_connected_subgraphs, user_callback);
}

// ==========================================================================

// Enumerates the largest, unique common subgraphs found between
// graph1 and graph2.  Note that the ENTIRE search space is explored
// before user_callback is actually invoked.
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond,
          concepts::ReadableVertexIndexMap<GraphFirst> VertexIndexMapFirst,
          concepts::ReadableVertexIndexMap<GraphSecond> VertexIndexMapSecond, typename EdgeEquivalencePredicate,
          typename VertexEquivalencePredicate, typename SubGraphCallback>
void mcgregor_common_subgraphs_maximum_unique(const GraphFirst& graph1, const GraphSecond& graph2,
                                              const VertexIndexMapFirst vindex_map1,      // get(vertex_index, graph1)
                                              const VertexIndexMapSecond vindex_map2,     // get(vertex_index, graph1)
                                              EdgeEquivalencePredicate edges_equivalent,  // always_equivalent()
                                              VertexEquivalencePredicate vertices_equivalent,  // always_equivalent()
                                              bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_detail::maximum_subgraph_interceptor unique_max_interceptor(
      graph1, graph2, vindex_map1, vindex_map2, user_callback, /*CheckUniqueness=*/std::true_type{});

  mcgregor_detail::mcgregor_common_subgraphs_internal_init(graph1, graph2, vindex_map1, vindex_map2, edges_equivalent,
                                                           vertices_equivalent, only_connected_subgraphs,
                                                           unique_max_interceptor);

  // Only output the largest, unique subgraphs
  unique_max_interceptor.output_subgraphs();
}

// Variant of mcgregor_common_subgraphs_maximum_unique with all default
// parameters
template <concepts::VertexListGraph GraphFirst, concepts::VertexListGraph GraphSecond, typename SubGraphCallback>
void mcgregor_common_subgraphs_maximum_unique(const GraphFirst& graph1, const GraphSecond& graph2,
                                              bool only_connected_subgraphs, SubGraphCallback user_callback) {
  mcgregor_common_subgraphs_maximum_unique(graph1, graph2, get(vertex_index, graph1), get(vertex_index, graph2),
                                           always_equivalent(), always_equivalent(), only_connected_subgraphs,
                                           user_callback);
}

// ==========================================================================

// Fills a membership map (vertex -> bool) using the information
// present in correspondence_map_1_to_2. Every vertex in a
// membership map will have a true value only if it is not
// associated with a null vertex in the correspondence map.
template <concepts::Graph GraphSecond, concepts::VertexListGraph GraphFirst,
          concepts::ReadableVertexPropertyMap<GraphFirst> CorrespondenceMapFirstToSecond,
          concepts::WritableVertexPropertyMap<GraphFirst> MembershipMapFirst>
void fill_membership_map(const GraphFirst& graph1, const CorrespondenceMapFirstToSecond correspondence_map_1_to_2,
                         MembershipMapFirst membership_map1) {
  for (auto vertex1 : vertices(graph1)) {
    put(membership_map1, vertex1, get(correspondence_map_1_to_2, vertex1) != graph_traits<GraphSecond>::null_vertex());
  }
}

// Returns a filtered sub-graph of graph whose edge and vertex
// inclusion is dictated by membership_map.
template <concepts::Graph Graph, concepts::ReadableVertexPropertyMap<Graph> MembershipMap>
auto make_membership_filtered_graph(const Graph& graph, MembershipMap& membership_map) {
  return filtered_graph{graph, keep_all(), property_map_filter{membership_map}};
}

}  // namespace bagl

#endif  // BAGL_BAGL_MCGREGOR_COMMON_SUBGRAPHS_H_
