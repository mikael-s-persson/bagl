// Copyright 2008 Aaron Windsor
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/graph_traits.h"
#include "bagl/is_straight_line_drawing.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/vector_property_map.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct Coord {
  std::size_t x;
  std::size_t y;
};
using Graph = adjacency_list<vec_s, vec_s, undirected_s>;

TEST(IsStraightLineDrawTest, Simple) {
  Graph g(4);
  add_edge(0, 1, g);
  add_edge(2, 3, g);

  auto drawing = vector_property_map(num_vertices(g), get(vertex_index, g), Coord{});

  // two perpendicular lines that intersect at (1,1)
  drawing[0].x = 1;
  drawing[0].y = 0;
  drawing[1].x = 1;
  drawing[1].y = 2;
  drawing[2].x = 0;
  drawing[2].y = 1;
  drawing[3].x = 2;
  drawing[3].y = 1;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));

  // two parallel horizontal lines
  drawing[0].x = 0;
  drawing[0].y = 0;
  drawing[1].x = 2;
  drawing[1].y = 0;

  EXPECT_TRUE(is_straight_line_drawing(g, drawing.ref()));

  // two parallel vertical lines
  drawing[0].x = 0;
  drawing[0].y = 0;
  drawing[1].x = 0;
  drawing[1].y = 2;
  drawing[2].x = 1;
  drawing[2].y = 0;
  drawing[3].x = 1;
  drawing[3].y = 2;

  EXPECT_TRUE(is_straight_line_drawing(g, drawing.ref()));

  // two lines that intersect at (1,1)
  drawing[0].x = 0;
  drawing[0].y = 0;
  drawing[1].x = 2;
  drawing[1].y = 2;
  drawing[2].x = 0;
  drawing[2].y = 2;
  drawing[3].x = 2;
  drawing[3].y = 0;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));
}

TEST(IsStraightLineDrawTest, N4) {
  // K_4 arranged in a diamond pattern, so that edges intersect
  Graph g(4);
  add_edge(0, 1, g);
  add_edge(0, 2, g);
  add_edge(0, 3, g);
  add_edge(1, 2, g);
  add_edge(1, 3, g);
  add_edge(2, 3, g);

  auto drawing = vector_property_map(num_vertices(g), get(vertex_index, g), Coord{});

  drawing[0].x = 1;
  drawing[0].y = 2;
  drawing[1].x = 2;
  drawing[1].y = 1;
  drawing[2].x = 1;
  drawing[2].y = 0;
  drawing[3].x = 0;
  drawing[3].y = 1;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));

  // K_4 arranged so that no edges intersect
  drawing[0].x = 0;
  drawing[0].y = 0;
  drawing[1].x = 1;
  drawing[1].y = 1;
  drawing[2].x = 1;
  drawing[2].y = 2;
  drawing[3].x = 2;
  drawing[3].y = 0;

  EXPECT_TRUE(is_straight_line_drawing(g, drawing.ref()));
}

TEST(IsStraightLineDrawTest, N8) {
  // a slightly more complicated example - edges (0,1) and (4,5)
  // intersect
  Graph g(8);
  add_edge(0, 1, g);
  add_edge(2, 3, g);
  add_edge(4, 5, g);
  add_edge(6, 7, g);

  auto drawing = vector_property_map(num_vertices(g), get(vertex_index, g), Coord{});

  drawing[0].x = 1;
  drawing[0].y = 1;
  drawing[1].x = 5;
  drawing[1].y = 4;
  drawing[2].x = 2;
  drawing[2].y = 5;
  drawing[3].x = 4;
  drawing[3].y = 4;
  drawing[4].x = 3;
  drawing[4].y = 4;
  drawing[5].x = 3;
  drawing[5].y = 2;
  drawing[6].x = 4;
  drawing[6].y = 2;
  drawing[7].x = 1;
  drawing[7].y = 1;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));
}

TEST(IsStraightLineDrawTest, N22) {
  // form a graph consisting of a bunch of parallel vertical edges,
  // then place an edge at various positions to intersect edges
  Graph g(22);
  for (int i = 0; i < 11; ++i) {
    add_edge(2 * i, 2 * i + 1, g);
  }

  auto drawing = vector_property_map(num_vertices(g), get(vertex_index, g), Coord{});

  for (int i = 0; i < 10; ++i) {
    drawing[2 * i].x = i;
    drawing[2 * i].y = 0;
    drawing[2 * i + 1].x = i;
    drawing[2 * i + 1].y = 10;
  }

  // put the final edge as a horizontal edge intersecting one other edge
  drawing[20].x = 5;
  drawing[20].y = 5;
  drawing[21].x = 7;
  drawing[21].y = 5;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));

  // make the final edge a diagonal intersecting multiple edges
  drawing[20].x = 2;
  drawing[20].y = 4;
  drawing[21].x = 9;
  drawing[21].y = 7;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));

  // reverse the slope
  drawing[20].x = 2;
  drawing[20].y = 7;
  drawing[21].x = 9;
  drawing[21].y = 4;

  EXPECT_FALSE(is_straight_line_drawing(g, drawing.ref()));
}

}  // namespace
}  // namespace bagl
