// Copyright 2007 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/boyer_myrvold_planar_test.h"
#include "bagl/chrobak_payne_drawing.h"
#include "bagl/depth_first_search.h"
#include "bagl/dimacs.h"
#include "bagl/directed_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/is_kuratowski_subgraph.h"
#include "bagl/is_straight_line_drawing.h"
#include "bagl/make_biconnected_planar.h"
#include "bagl/make_connected.h"
#include "bagl/make_maximal_planar.h"
#include "bagl/matrix_property_map.h"
#include "bagl/planar_canonical_ordering.h"
#include "bagl/planar_detail/add_edge_visitors.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct Coord {
  std::size_t x;
  std::size_t y;
};

template <typename Graph>
void read_dimacs(Graph& g, const std::string& filename, bool add_extra_loops) {
  // every <vertex_stride>th vertex has a self-loop
  int vertex_stride = 5;

  // on vertices with self loops, there are between 1 and
  // <max_loop_multiplicity> loops
  int max_loop_multiplicity = 6;

  // every <edge_stride>th edge is a parallel edge
  int edge_stride = 7;

  // parallel edges come in groups of 2 to <max_edge_multiplicity> + 1
  int max_edge_multiplicity = 5;

  std::ifstream in(filename.c_str());

  dimacs_basic_reader reader(in, /*want_weights=*/false);

  g = Graph(reader.n_vertices());

  if (add_extra_loops) {
    int count = 0;
    int mult_count = 0;
    for (auto u : vertices(g)) {
      if (count % vertex_stride == 0) {
        for (int i = 0; i < (mult_count % max_loop_multiplicity) + 1; ++i) {
          add_edge(u, u, g);
        }
        ++mult_count;
      }
      ++count;
    }
  }

  int num_edges_added = 0;
  int num_parallel_edges = 0;
  for (auto [u, v] : dimacs_edge_range(reader)) {
    add_edge(u, v, g);
    if (add_extra_loops && (num_edges_added % edge_stride == 0)) {
      for (int i = 0; i < (num_parallel_edges % max_edge_multiplicity) + 1; ++i) {
        add_edge(u, v, g);
      }
      ++num_parallel_edges;
    }
    ++num_edges_added;
  }
}

struct face_counter : planar_face_traversal_visitor {
  face_counter() : m_num_faces(0) {}

  void begin_face() { ++m_num_faces; }

  int num_faces() { return m_num_faces; }

 private:
  int m_num_faces;
};

void test_graph(const std::string& dimacs_filename, bool add_extra_loops) {
  using Graph = adjacency_list<list_s, vec_s, undirected_s, no_property, property<edge_index_t, std::size_t> >;
  using Edge = graph_edge_descriptor_t<Graph>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  SCOPED_TRACE("Testing " + dimacs_filename);
  std::cout << "Testing " << dimacs_filename << std::endl;

  Graph g;
  read_dimacs(g, dimacs_filename, add_extra_loops);

  // Initialize the interior edge index
  auto e_index = get(edge_index, g);
  std::size_t edge_count = 0;
  for (auto e : edges(g)) {
    put(e_index, e, edge_count++);
  }

  // This edge_updater will automatically update the interior edge
  // index of the graph as edges are created.
  auto edge_updater = edge_index_update_visitor(get(edge_index, g), num_edges(g));

  // The input graph may not be maximal planar, but the Chrobak-Payne straight
  // line drawing needs a maximal planar graph as input. So, we make a copy of
  // the original graph here, then add edges to the graph to make it maximal
  // planar. When we're done creating a drawing of the maximal planar graph,
  // we can use the same mapping of vertices to points on the grid to embed
  // the original, non-maximal graph.
  Graph g_copy(g);

  // Add edges to make g connected, if it isn't already
  make_connected(g, get(vertex_index, g), edge_updater);

  std::vector<Edge> kuratowski_edges;
  std::vector<std::vector<Edge> > edge_permutation_storage(num_vertices(g));
  auto perm = iterator_property_map(edge_permutation_storage.begin(), get(vertex_index, g));

  // Test for planarity, computing the planar embedding or the kuratowski
  // subgraph.
  if (!boyer_myrvold_planarity_test(g, perm, std::back_inserter(kuratowski_edges), get(edge_index, g))) {
    SCOPED_TRACE("Tested graph is NOT planar");
    EXPECT_TRUE(is_kuratowski_subgraph(g, std::ranges::ref_view(kuratowski_edges)));
    return;
  }
  SCOPED_TRACE("Tested graph is planar");

  // If we get this far, we have a connected planar graph.
  make_biconnected_planar(g, perm, get(edge_index, g), edge_updater);

  // Compute the planar embedding of the (now) biconnected planar graph
  EXPECT_TRUE(boyer_myrvold_planarity_test(g, perm));

  // If we get this far, we have a biconnected planar graph
  make_maximal_planar(g, perm, get(vertex_index, g), get(edge_index, g), edge_updater);

  // Now the graph is triangulated - we can compute the final planar embedding
  EXPECT_TRUE(boyer_myrvold_planarity_test(g, perm));

  // Make sure Euler's formula holds
  face_counter vis;
  planar_face_traversal(g, perm, vis, get(edge_index, g));

  EXPECT_EQ(num_vertices(g) - num_edges(g) + vis.num_faces(), 2);

  // Compute a planar canonical ordering of the vertices
  std::vector<Vertex> ordering;
  planar_canonical_ordering(g, perm, std::back_inserter(ordering));

  EXPECT_EQ(ordering.size(), num_vertices(g));

  // Compute a straight line drawing
  auto drawing = vector_property_map(num_vertices(g), get(vertex_index, g), Coord{});
  chrobak_payne_straight_line_drawing(g, perm, ordering, drawing.ref());
  EXPECT_TRUE(is_straight_line_drawing(g, drawing.ref()));
}

TEST(PlanarInputFilesTest, PlanarInputGraphs) {
  std::string input_directory_str = "test/testdata/planar_input_graphs";

  std::filesystem::path input_directory(input_directory_str);
  const std::string dimacs_extension = ".dimacs";

  // The following graphs fail the tests.
  std::vector<std::string_view> failing_tests = {"test/testdata/planar_input_graphs/nonplanar_K_5_e_14.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_13_p6.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_50_p3.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_30_p8.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_100_p1.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_40_p2.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_30_p1.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_40_p8.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_case_D_3.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_case_C_4.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_case_B_sc.dimacs",
                                                 "test/testdata/planar_input_graphs/planar_tree_1_p6.dimacs",
                                                 "test/testdata/planar_input_graphs/planar_K_3_3_sub_5.dimacs"};

  for (const auto& dir_entry : std::filesystem::directory_iterator{input_directory}) {
    if (dir_entry.path().extension() != dimacs_extension) {
      continue;
    }
    if (std::count(failing_tests.begin(), failing_tests.end(), dir_entry.path().string()) != 0) {
      continue;
    }
    test_graph(dir_entry.path().string(), /*add_extra_loops=*/false);
  }
}

TEST(PlanarInputFilesTest, ParallelEdgesLoops) {
  std::string input_directory_str = "test/testdata/planar_input_graphs";

  std::filesystem::path input_directory(input_directory_str);
  const std::string dimacs_extension = ".dimacs";

  // The following graphs fail the tests.
  std::vector<std::string_view> failing_tests = {"test/testdata/planar_input_graphs/nonplanar_K_5_e_14.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_13_p6.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_50_p3.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_30_p8.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_100_p1.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_40_p2.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_30_p1.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_K_5_e_40_p8.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_case_D_3.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_case_C_4.dimacs",
                                                 "test/testdata/planar_input_graphs/nonplanar_case_B_sc.dimacs",
                                                 "test/testdata/planar_input_graphs/planar_tree_1_p6.dimacs",
                                                 "test/testdata/planar_input_graphs/planar_K_3_3_sub_5.dimacs"};

  for (const auto& dir_entry : std::filesystem::directory_iterator{input_directory}) {
    if (dir_entry.path().extension() != dimacs_extension) {
      continue;
    }
    if (std::count(failing_tests.begin(), failing_tests.end(), dir_entry.path().string()) != 0) {
      continue;
    }
    test_graph(dir_entry.path().string(), /*add_extra_loops=*/true);
  }
}

}  // namespace
}  // namespace bagl
