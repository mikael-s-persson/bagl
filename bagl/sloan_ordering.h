// Copyright 2002 Marc Wintermantel (wintermantel@even-ag.ch)
// ETH Zurich, Center of Structure Technologies
// (https://web.archive.org/web/20050307090307/http://www.structures.ethz.ch/)
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_SLOAN_ORDERING_H_
#define BAGL_BAGL_SLOAN_ORDERING_H_

#include <algorithm>
#include <limits>
#include <queue>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/breadth_first_search.h"
#include "bagl/cuthill_mckee_ordering.h"
#include "bagl/graph_traits.h"
#include "bagl/properties.h"
#include "bagl/property_map.h"
#include "bagl/visitors.h"

//==========================================================
// Sloan-Algorithm for graph reordering
//(optimzes profile and wavefront, not primiraly bandwidth
//==========================================================

namespace bagl {

constexpr int sloan_default_weight1 = 1;  // default weight for the distance in the Sloan algorithm
constexpr int sloan_default_weight2 = 2;  // default weight for the degree in the Sloan algorithm

// Function that returns the maximum depth of
// a rooted level strucutre (rls)
template <class Distance>
auto rls_depth(Distance& d) {
  using Degree = typename Distance::value_type;

  Degree h_s = 0;
  for (auto i : d) {
    if (i > h_s) {
      h_s = i;
    }
  }
  return h_s;
}

// Function that returns the width of the largest level of
// a rooted level strucutre (rls)
template <class Distance>
auto rls_max_width(Distance& d, std::size_t depth) {
  using Degree = typename Distance::value_type;

  // Searching for the maximum width of a level
  std::vector<Degree> dummy_width(depth + 1, 0);
  for (auto w : d) {
    dummy_width[w]++;
  }

  Degree w_max = 0;
  for (auto w : dummy_width) {
    if (w > w_max) {
      w_max = w;
    }
  }
  return w_max;
}

// Function for finding a good starting node for Sloan algorithm
//
// This is to find a good starting node. "good" is in the sense
// of the ordering generated.
template <class Graph, class ColorMap, class DegreeMap>
auto sloan_start_end_vertices(Graph& g, graph_vertex_descriptor_t<Graph>& s, ColorMap color, DegreeMap degree) {
  using Degree = property_traits_value_t<DegreeMap>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Creating a property_map for the indices of a vertex
  auto index_map = get(vertex_index, g);

  // Creating a std-vector for storing the distance from the start vertex in
  // dist
  std::vector<std::size_t> dist(num_vertices(g), 0);
  auto dist_pmap = iterator_property_map(dist.begin(), index_map);

  // Creating a priority queue
  auto comp = [degree](Vertex lhs, Vertex rhs) { return get(degree, lhs) > get(degree, rhs); };
  std::priority_queue<Vertex, std::vector<Vertex>, decltype(comp)> degree_queue(comp);

  // step 1
  // Scan for the vertex with the smallest degree and the maximum degree
  s = *vertices(g).begin();
  Degree my_degree = get(degree, s);
  Degree maximum_degree = 0;
  for (auto u : vertices(g)) {
    Degree dummy = get(degree, u);

    if (dummy < my_degree) {
      my_degree = dummy;
      s = u;
    }

    if (dummy > maximum_degree) {
      maximum_degree = dummy;
    }
  }
  // end 1

  bool new_start = true;
  Vertex e = s;
  do {
    new_start = false;  // Setting the loop repetition status to false

    // step 2
    // initialize the the disance std-vector with 0
    std::ranges::fill(dist, 0);

    // generating the RLS (rooted level structure)
    breadth_first_search(g, std::ranges::single_view{s}, make_bfs_visitor(distance_recorder_on_tree_edge(dist_pmap)));
    // end 2

    // step 3
    // calculating the depth of the RLS
    Degree h_s = rls_depth(dist);

    // step 4
    // pushing one node of each degree in an ascending manner into
    // degree_queue
    std::vector<bool> shrink_trace(maximum_degree, false);
    for (auto u : vertices(g)) {
      Degree dummy = get(degree, u);
      if ((dist[index_map[u]] == h_s) && (!shrink_trace[dummy])) {
        degree_queue.push(u);
        shrink_trace[dummy] = true;
      }
    }
    // end 3 & 4

    // step 5
    // Initializing w
    Degree w_e = std::numeric_limits<Degree>::max();
    // end 5

    // step 6
    // Testing for termination
    while (!degree_queue.empty()) {
      Vertex i = degree_queue.top();  // getting the node with the lowest degree
                                      // from the degree queue
      degree_queue.pop();             // ereasing the node with the lowest degree from
                                      // the degree queue

      // generating a RLS
      std::ranges::fill(dist, 0);
      breadth_first_search(g, std::ranges::single_view{i}, make_bfs_visitor(distance_recorder_on_tree_edge(dist_pmap)));

      // Calculating depth and width of the rooted level
      Degree h_i = rls_depth(dist);
      Degree w_i = rls_max_width(dist, h_i);

      // Testing for termination
      if ((h_i > h_s) && (w_i < w_e)) {
        h_s = h_i;
        s = i;
        while (!degree_queue.empty()) {
          degree_queue.pop();
        }
        new_start = true;
      } else if (w_i < w_e) {
        w_e = w_i;
        e = i;
      }
    }
    // end 6

  } while (new_start);

  return e;
}

// Sloan algorithm with a given starting Vertex.
//
// This algorithm requires user to provide a starting vertex to
// compute Sloan ordering.
template <class Graph, class OutputIterator, class ColorMap, class DegreeMap, class PriorityMap, class Weight>
OutputIterator sloan_ordering(Graph& g, graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> e,
                              OutputIterator permutation, ColorMap color, DegreeMap degree, PriorityMap priority,
                              Weight w1, Weight w2) {
  using Degree = property_traits_value_t<PriorityMap>;
  using ColorValue = property_traits_value_t<ColorMap>;
  using Color = color_traits<ColorValue>;
  using Vertex = graph_vertex_descriptor_t<Graph>;

  // Creating a property_map for the indices of a vertex
  auto index_map = get(vertex_index, g);

  // Creating a std-vector for storing the distance from the end vertex in it
  std::vector<std::size_t> dist(num_vertices(g), 0);
  auto dist_pmap = iterator_property_map(dist.begin(), index_map);

  breadth_first_search(g, std::ranges::single_view{e}, make_bfs_visitor(distance_recorder_on_tree_edge(dist_pmap)));

  // Sets the color and priority to their initial status
  for (auto u : vertices(g)) {
    put(color, u, Color::white());
    Degree cdeg = get(degree, u) + 1;
    put(priority, u, w1 * dist[index_map[u]] - w2 * cdeg);
  }

  // Priority list
  auto comp = [priority](Vertex lhs, Vertex rhs) { return get(priority, lhs) > get(priority, rhs); };
  std::list<Vertex> priority_list;

  // Sets the color of the starting vertex to gray
  put(color, s, Color::gray());
  priority_list.push_front(s);  // Puts s into the priority_list

  while (!priority_list.empty()) {
    // Orders the elements in the priority list in an ascending manner
    priority_list.sort(comp);

    Vertex u = priority_list.front();
    priority_list.pop_front();

    if (get(color, u) == Color::gray()) {
      for (auto oe : out_edges(u, g)) {
        Vertex v = target(oe, g);

        put(priority, v, get(priority, v) + w2);

        // test if the vertex is inactive
        if (get(color, v) == Color::white()) {
          put(color, v, Color::gray());
          priority_list.push_front(v);
        }
      }
    }

    // Here starts step 8
    *permutation++ = u;
    // Gives u an inactive status
    put(color, u, Color::black());

    for (auto oe : out_edges(u, g)) {
      Vertex v = target(oe, g);

      if (get(color, v) != Color::gray()) {
        continue;
      }

      put(color, v,
          Color::red());                        // giving the vertex an active status
      put(priority, v, get(priority, v) + w2);  // updates the priority

      // for loop over alll adjacent vertices of v
      for (auto e2 : out_edges(v, g)) {
        Vertex w = target(e2, g);

        if (get(color, w) == Color::black()) {
          continue;
        }

        put(priority, w, get(priority, w) + w2);  // updates the priority

        if (get(color, w) != Color::white()) {
          continue;
        }

        put(color, w, Color::gray());
        priority_list.push_front(w);
      }
    }
  }

  return permutation;
}

// Same algorithm as before, but without the weights given (taking default
// weights
template <class Graph, class OutputIterator, class ColorMap, class DegreeMap, class PriorityMap>
OutputIterator sloan_ordering(Graph& g, graph_vertex_descriptor_t<Graph> s, graph_vertex_descriptor_t<Graph> e,
                              OutputIterator permutation, ColorMap color, DegreeMap degree, PriorityMap priority) {
  return sloan_ordering(g, s, e, permutation, color, degree, priority, sloan_default_weight1, sloan_default_weight2);
}

// Sloan algorithm without a given starting Vertex.
//
// This algorithm finds a good starting vertex itself to
// compute Sloan-ordering.

template <class Graph, class OutputIterator, class Color, class Degree, class Priority, class Weight>
OutputIterator sloan_ordering(Graph& G, OutputIterator permutation, Color color, Degree degree, Priority priority,
                              Weight w1, Weight w2) {
  using Vertex = graph_vertex_descriptor_t<Graph>;

  Vertex s = graph_traits<Graph>::null_vertex();
  Vertex e = sloan_start_end_vertices(G, s, color, degree);

  return sloan_ordering(G, s, e, permutation, color, degree, priority, w1, w2);
}

/////////////////////////////////////////////////////////////////////////////////////////
// Same as before, but without given weights (default weights are taken instead)
template <class Graph, class OutputIterator, class Color, class Degree, class Priority>
OutputIterator sloan_ordering(Graph& G, OutputIterator permutation, Color color, Degree degree, Priority priority) {
  return sloan_ordering(G, permutation, color, degree, priority, sloan_default_weight1, sloan_default_weight2);
}

}  // namespace bagl

#endif  // BAGL_BAGL_SLOAN_ORDERING_H_
