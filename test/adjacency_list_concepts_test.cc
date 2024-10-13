// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/adjacency_list.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_traits.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

TEST(AdjacencyList, ConceptChecks) {
  // Check adjacency_list with properties
  {
    using Graph =
        adjacency_list<vec_s, vec_s, directed_s, property<vertex_color_t, int>, property<edge_weight_t, int> >;
    using Vertex = graph_vertex_descriptor_t<Graph>;
    using Edge = graph_edge_descriptor_t<Graph>;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  {
    typedef adjacency_list<vec_s, vec_s, bidirectional_s, property<vertex_color_t, int>, property<edge_weight_t, int> >
        Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  {
    typedef adjacency_list<list_s, list_s, directed_s, property<vertex_color_t, int>, property<edge_weight_t, int> >
        Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  {
    typedef adjacency_list<list_s, list_s, undirected_s, property<vertex_color_t, int>, property<edge_weight_t, int> >
        Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  // Checking adjacency_list with EdgeList=set_s
  {
    typedef adjacency_list<set_s, vec_s, bidirectional_s, property<vertex_color_t, int>, property<edge_weight_t, int> >
        Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  {
    typedef adjacency_list<set_s, list_s, directed_s, property<vertex_color_t, int>, property<edge_weight_t, int> >
        Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  {
    typedef adjacency_list<set_s, list_s, undirected_s, property<vertex_color_t, int>, property<edge_weight_t, int> >
        Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    typedef graph_traits<Graph>::edge_descriptor Edge;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Vertex, vertex_color_t>);
    static_assert(concepts::LvaluePropertyGraph<Graph, Edge, edge_weight_t>);
  }
  // Check adjacency_list without any properties
  {
    typedef adjacency_list<vec_s, vec_s, directed_s> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::VertexMutablePropertyGraph<Graph>);
    static_assert(concepts::EdgeMutablePropertyGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
  }
  {
    typedef adjacency_list<vec_s, vec_s, bidirectional_s> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
  }
  {
    typedef adjacency_list<list_s, list_s, directed_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<list_s, list_s, bidirectional_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<list_s, list_s, undirected_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  // Checking EdgeList=set_s with no properties
  {
    typedef adjacency_list<set_s, vec_s, bidirectional_s> Graph;
    typedef graph_traits<Graph>::vertex_descriptor Vertex;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
    static_assert(concepts::ReadablePropertyGraph<Graph, Vertex, vertex_index_t>);
  }
  {
    typedef adjacency_list<set_s, list_s, directed_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<set_s, list_s, bidirectional_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<set_s, list_s, undirected_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<vec_s, pool_s, directed_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<vec_s, pool_s, bidirectional_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<vec_s, pool_s, undirected_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<pool_s, pool_s, directed_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableIncidenceGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<pool_s, pool_s, bidirectional_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::BidirectionalGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
  {
    typedef adjacency_list<pool_s, pool_s, undirected_s> Graph;
    static_assert(concepts::VertexListGraph<Graph>);
    static_assert(concepts::EdgeListGraph<Graph>);
    static_assert(concepts::IncidenceGraph<Graph>);
    static_assert(concepts::AdjacencyGraph<Graph>);
    static_assert(concepts::MutableBidirectionalGraph<Graph>);
    static_assert(concepts::MutableEdgeListGraph<Graph>);
  }
}

}  // namespace
}  // namespace bagl
