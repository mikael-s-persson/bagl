// Copyright 2002 Indiana University.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

/* Expected Output
         0    1    2    3    4    5    6    7    8    9
    0    0   99  111  123  103  107  125  105  111  123
    1   99    0   12   24    4    8   26    6   12   24
    2  111   12    0   12   16   20   24   18   24   26
    3  123   24   12    0   28   30   12   30   26   14
    4  103    4   16   28    0    4   22    2    8   20
    5  107    8   20   30    4    0   18    6    4   16
    6  125   26   24   12   22   18    0   24   14    2
    7  105    6   18   30    2    6   24    0   10   22
    8  111   12   24   26    8    4   14   10    0   12
    9  123   24   26   14   20   16    2   22   12    0
*/

#include <iomanip>
#include <iostream>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/johnson_all_pairs_shortest.h"
#include "bagl/property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(JohnsonAllPairsTest, Basic) {
  using Graph = adjacency_list<vec_s, vec_s, undirected_s, no_property,
                               property<edge_weight_t, int, property<edge_weight2_t, int> > >;
  using Edge = std::pair<std::size_t, std::size_t>;

  std::array edge_array = {Edge(0, 1), Edge(1, 2), Edge(2, 3), Edge(1, 4), Edge(2, 5), Edge(3, 6), Edge(4, 5),
                           Edge(5, 6), Edge(4, 7), Edge(5, 8), Edge(6, 9), Edge(7, 8), Edge(8, 9)};
  const std::size_t V = 10;

  Graph g(V, edge_array);

  auto w = get(edge_weight, g);
  std::array weights = {99, 12, 12, 4, 99, 12, 4, 99, 2, 4, 2, 99, 12};

  for (auto [e, value] : zip_range(edges(g), weights)) {
    w[e] = value;
  }

  int D[V][V];
  EXPECT_TRUE(johnson_all_pairs_shortest_paths(g, D));

  // Check symmetry
  for (int i = 0; i < V; ++i) {
    EXPECT_EQ(D[i][i], 0);
    for (int j = i + 1; j < V; ++j) {
      EXPECT_EQ(D[i][j], D[j][i]);
    }
  }

  if (::testing::Test::HasFailure()) {
    std::cout << "Got following distance matrix:\n";
    std::cout << std::setw(5) << " ";
    for (int k = 0; k < 10; ++k) {
      std::cout << std::setw(5) << k;
    }
    std::cout << std::endl;
    for (int i = 0; i < 10; ++i) {
      std::cout << std::setw(5) << i;
      for (int j = 0; j < 10; ++j) {
        std::cout << std::setw(5) << D[i][j];
      }
      std::cout << std::endl;
    }
  }
}

}  // namespace
}  // namespace bagl
