// Copyright 2009 Trustees of Indiana University.
// Authors: Michael Hansen
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_utility.h"
#include "bagl/isomorphism.h"
#include "bagl/mcgregor_common_subgraphs.h"
#include "bagl/random.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

// Callback that compares incoming graphs to the supplied common
// subgraph.
template <typename Graph>
struct test_callback {
  test_callback(Graph& common_subgraph, const Graph& graph1, const Graph& graph2, bool& was_common_subgraph_found,
                std::vector<std::string>& all_subgraphs_found)
      : graph1_(graph1),
        graph2_(graph2),
        common_subgraph_(common_subgraph),
        was_common_subgraph_found_(&was_common_subgraph_found),
        all_subgraphs_found_(&all_subgraphs_found) {}

  template <typename CorrFirstToSecond, typename CorrSecondToFirst>
  bool operator()(CorrFirstToSecond corr_1_to_2, CorrSecondToFirst corr_2_to_1, std::size_t subgraph_size) {
    using Vertex = graph_vertex_descriptor_t<Graph>;

    if (subgraph_size != num_vertices(common_subgraph_)) {
      return true;
    }

    // Fill membership maps for both graphs
    auto membership_map1 = vector_property_map(num_vertices(graph1_), get(vertex_index, graph1_), false);
    auto membership_map2 = vector_property_map(num_vertices(graph2_), get(vertex_index, graph2_), false);

    fill_membership_map<Graph>(graph1_, corr_1_to_2, membership_map1.ref());
    fill_membership_map<Graph>(graph2_, corr_2_to_1, membership_map2.ref());

    // Generate filtered graphs using membership maps
    auto subgraph1 = make_membership_filtered_graph(graph1_, membership_map1.ref());
    auto subgraph2 = make_membership_filtered_graph(graph2_, membership_map2.ref());

    // For debug print in case we fail.
    all_subgraphs_found_->emplace_back(
        print_graph_to_string(subgraph1, get(vertex_name, graph1_), get(edge_name, graph1_)));

    auto vindex_map1 = get(vertex_index, subgraph1);
    auto vindex_map2 = get(vertex_index, subgraph2);

    auto vname_map_common = get(vertex_name, common_subgraph_);
    auto vname_map1 = get(vertex_name, subgraph1);
    auto vname_map2 = get(vertex_name, subgraph2);

    auto ename_map_common = get(edge_name, common_subgraph_);
    auto ename_map1 = get(edge_name, subgraph1);
    auto ename_map2 = get(edge_name, subgraph2);

    // Verify that subgraph1 matches the supplied common subgraph
    for (auto vertex1 : vertices(subgraph1)) {
      Vertex vertex_common = *std::next(vertices(common_subgraph_).begin(), get(vindex_map1, vertex1));

      // Match vertex names
      if (get(vname_map_common, vertex_common) != get(vname_map1, vertex1)) {
        // Keep looking
        return true;
      }

      for (auto vertex1_2 : vertices(subgraph1)) {
        Vertex vertex_common2 = *std::next(vertices(common_subgraph_).begin(), get(vindex_map1, vertex1_2));
        auto [edge_common, edge_common_found] = edge(vertex_common, vertex_common2, common_subgraph_);
        auto [edge1, edge1_found] = edge(vertex1, vertex1_2, subgraph1);

        if ((edge_common_found != edge1_found) ||
            ((edge_common_found && edge1_found) && (get(ename_map_common, edge_common) != get(ename_map1, edge1)))) {
          // Keep looking
          return true;
        }
      }
    }

    // Verify that subgraph2 matches the supplied common subgraph
    for (auto vertex2 : vertices(subgraph2)) {
      Vertex vertex_common = *std::next(vertices(common_subgraph_).begin(), get(vindex_map2, vertex2));

      // Match vertex names
      if (get(vname_map_common, vertex_common) != get(vname_map2, vertex2)) {
        // Keep looking
        return true;
      }

      for (auto vertex2_2 : vertices(subgraph2)) {
        Vertex vertex_common2 = *std::next(vertices(common_subgraph_).begin(), get(vindex_map2, vertex2_2));
        auto [edge_common, edge_common_found] = edge(vertex_common, vertex_common2, common_subgraph_);
        auto [edge2, edge2_found] = edge(vertex2, vertex2_2, subgraph2);

        if ((edge_common_found != edge2_found) ||
            ((edge_common_found && edge2_found) && (get(ename_map_common, edge_common) != get(ename_map2, edge2)))) {
          // Keep looking
          return true;
        }
      }
    }

    // Check isomorphism just to be thorough
    if (verify_isomorphism(subgraph1, subgraph2, corr_1_to_2)) {
      *was_common_subgraph_found_ = true;
      // Stop iterating
      return false;
    }

    // Keep looking
    return true;
  }

 private:
  const Graph& graph1_;
  const Graph& graph2_;
  Graph& common_subgraph_;

  bool* was_common_subgraph_found_;
  std::vector<std::string>* all_subgraphs_found_;
};

template <typename Graph, typename RandomNumberGenerator, typename VertexNameMap, typename EdgeNameMap>
void add_random_vertices(Graph& graph, RandomNumberGenerator& generator, int vertices_to_create,
                         int max_edges_per_vertex, VertexNameMap vname_map, EdgeNameMap ename_map) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  for (int v_index = 0; v_index < vertices_to_create; ++v_index) {
    Vertex new_vertex = add_vertex(graph);
    put(vname_map, new_vertex, generator());
  }

  // Add edges for every new vertex. Care is taken to avoid parallel
  // edges.
  for (auto source_vertex : vertices(graph)) {
    int edges_for_vertex =
        std::min(static_cast<int>(generator() % max_edges_per_vertex) + 1, static_cast<int>(num_vertices(graph)));

    while (edges_for_vertex > 0) {
      Vertex target_vertex = random_vertex(graph, generator);
      if (source_vertex == target_vertex || std::ranges::any_of(out_edges(source_vertex, graph), [&](auto e) {
            return target(e, graph) == target_vertex;
          })) {
        continue;
      }
      auto [new_e, new_e_created] = add_edge(source_vertex, target_vertex, graph);
      EXPECT_TRUE(new_e_created);
      put(ename_map, new_e, generator());
      edges_for_vertex--;
    }
  }
}

TEST(McGregorSubgraphsTest, RandomGraph) {
  int vertices_to_create = 10;
  int max_edges_per_vertex = 2;
  std::size_t random_seed = 42;

  std::mt19937 generator(random_seed);

  // Using a vec_s graph here so that we don't have to mess around with
  // a vertex index map; it will be implicit.
  using Graph = adjacency_list<list_s, vec_s, directed_s, property<vertex_name_t, unsigned int>,
                               property<edge_name_t, unsigned int> >;

  // Generate a random common subgraph and then add random vertices
  // and edges to the two parent graphs.
  Graph common_subgraph;
  Graph graph1;
  Graph graph2;

  auto vname_map_common = get(vertex_name, common_subgraph);
  auto vname_map1 = get(vertex_name, graph1);
  auto vname_map2 = get(vertex_name, graph2);

  auto ename_map_common = get(edge_name, common_subgraph);
  auto ename_map1 = get(edge_name, graph1);
  auto ename_map2 = get(edge_name, graph2);

  for (int vindex = 0; vindex < vertices_to_create; ++vindex) {
    put(vname_map_common, add_vertex(common_subgraph), generator());
  }

  for (auto source_vertex : vertices(common_subgraph)) {
    for (auto target_vertex : vertices(common_subgraph)) {
      if (source_vertex == target_vertex) {
        continue;
      }
      auto [new_e, new_e_created] = add_edge(source_vertex, target_vertex, common_subgraph);
      EXPECT_TRUE(new_e_created);
      put(ename_map_common, new_e, generator());
    }
  }

  randomize_property<vertex_name_t>(common_subgraph, generator);
  randomize_property<edge_name_t>(common_subgraph, generator);

  copy_graph(common_subgraph, graph1);
  copy_graph(common_subgraph, graph2);

  // Randomly add vertices and edges to graph1 and graph2.
  add_random_vertices(graph1, generator, vertices_to_create, max_edges_per_vertex, vname_map1, ename_map1);

  add_random_vertices(graph2, generator, vertices_to_create, max_edges_per_vertex, vname_map2, ename_map2);

  bool was_common_subgraph_found = false;
  std::vector<std::string> all_subgraphs_found;
  test_callback user_callback(common_subgraph, graph1, graph2, was_common_subgraph_found, all_subgraphs_found);

  mcgregor_common_subgraphs(graph1, graph2, get(vertex_index, graph1), get(vertex_index, graph2),
                            property_map_equivalent(ename_map1, ename_map2),
                            property_map_equivalent(vname_map1, vname_map2), true, user_callback);

  EXPECT_TRUE(was_common_subgraph_found);
  if (::testing::Test::HasFailure()) {
    std::cerr << "Graph1:\n"
              << print_graph_to_string(graph1, vname_map1, ename_map1) << "\nGraph2:\n"
              << print_graph_to_string(graph2, vname_map2, ename_map2) << "\nExpected Common Graph:\n"
              << print_graph_to_string(common_subgraph, vname_map_common, ename_map_common);
    for (auto& s : all_subgraphs_found) {
      std::cerr << "\n\nGot subgraph:\n" << s;
    }
  }
}

class McGregorSubgraphsSimpleTest : public ::testing::Test {
 public:
  using Graph = adjacency_list<list_s, vec_s, directed_s, property<vertex_name_t, unsigned int>,
                               property<edge_name_t, unsigned int> >;

  McGregorSubgraphsSimpleTest() {
    auto vname_map_simple1 = get(vertex_name, graph_simple1);
    auto vname_map_simple2 = get(vertex_name, graph_simple2);

    put(vname_map_simple1, add_vertex(graph_simple1), 1);
    put(vname_map_simple1, add_vertex(graph_simple1), 2);
    put(vname_map_simple1, add_vertex(graph_simple1), 3);

    add_edge(0, 1, graph_simple1);
    add_edge(0, 2, graph_simple1);
    add_edge(1, 2, graph_simple1);

    put(vname_map_simple2, add_vertex(graph_simple2), 1);
    put(vname_map_simple2, add_vertex(graph_simple2), 2);
    put(vname_map_simple2, add_vertex(graph_simple2), 3);
    put(vname_map_simple2, add_vertex(graph_simple2), 4);

    add_edge(0, 1, graph_simple2);
    add_edge(0, 2, graph_simple2);
    add_edge(1, 2, graph_simple2);
    add_edge(1, 3, graph_simple2);
  }

  auto vertex_name_equivalence() const {
    return property_map_equivalent(get(vertex_name, graph_simple1), get(vertex_name, graph_simple2));
  }

  template <typename CorrFirstToSecond>
  std::string make_subgraph_str(CorrFirstToSecond corr_1_to_2) const {
    std::stringstream subgraph_string;
    for (auto vertex1 : vertices(graph_simple1)) {
      auto vertex2 = get(corr_1_to_2, vertex1);
      if (vertex2 != graph_traits<Graph>::null_vertex()) {
        subgraph_string << vertex1 << "," << vertex2 << " ";
      }
    }
    return subgraph_string.str();
  }

 protected:
  Graph graph_simple1;
  Graph graph_simple2;
};

TEST_F(McGregorSubgraphsSimpleTest, UniqueSubGraph) {
  std::vector<std::string> simple_subgraph_list;
  auto user_callback = [this, &simple_subgraph_list](auto corr_1_to_2, const auto& /*unused*/, const auto& /*unused*/) {
    simple_subgraph_list.emplace_back(make_subgraph_str(corr_1_to_2));
    return true;
  };
  mcgregor_common_subgraphs_unique(graph_simple1, graph_simple2, get(vertex_index, graph_simple1),
                                   get(vertex_index, graph_simple2), always_equivalent(), vertex_name_equivalence(),
                                   true, user_callback);
  EXPECT_THAT(simple_subgraph_list, ::testing::IsSupersetOf({"0,0 1,1 ", "0,0 1,1 2,2 ", "0,0 2,2 ", "1,1 2,2 "}));
}

TEST_F(McGregorSubgraphsSimpleTest, MaximumSubGraph) {
  std::vector<std::string> simple_subgraph_list;
  auto user_callback = [this, &simple_subgraph_list](auto corr_1_to_2, const auto& /*unused*/, const auto& /*unused*/) {
    simple_subgraph_list.emplace_back(make_subgraph_str(corr_1_to_2));
    return true;
  };
  mcgregor_common_subgraphs_maximum(graph_simple1, graph_simple2, get(vertex_index, graph_simple1),
                                    get(vertex_index, graph_simple2), always_equivalent(), vertex_name_equivalence(),
                                    true, user_callback);

  EXPECT_THAT(simple_subgraph_list, ::testing::Contains("0,0 1,1 2,2 "));
}

TEST_F(McGregorSubgraphsSimpleTest, MaximumUniqueSubGraph) {
  std::vector<std::string> simple_subgraph_list;
  auto user_callback = [this, &simple_subgraph_list](auto corr_1_to_2, const auto& /*unused*/, const auto& /*unused*/) {
    simple_subgraph_list.emplace_back(make_subgraph_str(corr_1_to_2));
    return true;
  };
  mcgregor_common_subgraphs_maximum_unique(graph_simple1, graph_simple2, get(vertex_index, graph_simple1),
                                           get(vertex_index, graph_simple2), always_equivalent(),
                                           vertex_name_equivalence(), true, user_callback);

  EXPECT_THAT(simple_subgraph_list, ::testing::ElementsAre("0,0 1,1 2,2 "));
}

}  // namespace
}  // namespace bagl
