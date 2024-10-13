// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_matrix.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(AdjacencyMatrix, ConceptChecks) {
  // Check adjacency_matrix without properties
  {
    typedef adjacency_matrix<directed_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableGraph<Graph>);
    static_assert(concepts::AdjacencyMatrix<Graph>);
  }
  {
    typedef adjacency_matrix<undirected_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableGraph<Graph>);
    static_assert(concepts::AdjacencyMatrix<Graph>);
  }
  // Check adjacency_matrix with properties
  {
    typedef adjacency_matrix<directed_s, property<vertex_color_t, int>, property<edge_weight_t, float> > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::AdjacencyMatrix<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
    static_assert(concepts::PropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::PropertyGraph<Graph, Edge, edge_weight_t>);
  }
  {
    typedef adjacency_matrix<undirected_s, property<vertex_color_t, int>, property<edge_weight_t, float> > Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::AdjacencyMatrix<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
    static_assert(concepts::PropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::PropertyGraph<Graph, Edge, edge_weight_t>);
  }
}

}  // namespace
}  // namespace bagl