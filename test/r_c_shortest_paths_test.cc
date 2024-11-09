// Copyright Michael Drexl 2005, 2006.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include "bagl/r_c_shortest_paths.h"

#include <iostream>

#include "bagl/adjacency_list.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct VertProp {
  VertProp(int n = 0, int e = 0, int l = 0) : num(n), eat(e), lat(l) {}
  int num;
  // earliest arrival time
  int eat;
  // latest arrival time
  int lat;
};

struct ArcProp {
  ArcProp(int n = 0, int c = 0, int t = 0) : num(n), cost(c), time(t) {}
  int num;
  // traversal cost
  int cost;
  // traversal time
  int time;
};

using Graph = adjacency_list<vec_s, vec_s, directed_s, VertProp, ArcProp>;

// data structures for spp without resource constraints:
// ResourceContainer model
struct spp_no_rc_res_cont {
  spp_no_rc_res_cont(int c = 0) : cost(c){};
  int cost;
};

// ResourceExtensionFunction model
struct ref_no_res_cont {
  bool operator()(const Graph& g, spp_no_rc_res_cont& new_cont, const spp_no_rc_res_cont& old_cont,
                  graph_edge_descriptor_t<Graph> ed) const {
    new_cont.cost = old_cont.cost + g[ed].cost;
    return true;
  }
};

// DominanceFunction model
struct dominance_no_res_cont {
  bool operator()(const spp_no_rc_res_cont& res_cont_1, const spp_no_rc_res_cont& res_cont_2) const {
    // must be "<=" here!!!
    // must NOT be "<"!!!
    return res_cont_1.cost <= res_cont_2.cost;
    // this is not a contradiction to the documentation
    // the documentation says:
    // "A label $l_1$ dominates a label $l_2$ if and only if both are
    // resident at the same vertex, and if, for each resource, the resource
    // consumption of $l_1$ is less than or equal to the resource
    // consumption of $l_2$, and if there is at least one resource where
    // $l_1$ has a lower resource consumption than $l_2$." one can think of
    // a new label with a resource consumption equal to that of an old label
    // as being dominated by that old label, because the new one will have a
    // higher number and is created at a later point in time, so one can
    // implicitly use the number or the creation time as a resource for
    // tie-breaking
  }
};
// end data structures for spp without resource constraints:

// data structures for shortest path problem with time windows (spptw)
// ResourceContainer model
struct spp_spptw_res_cont {
  spp_spptw_res_cont(int c = 0, int t = 0) : cost(c), time(t) {}
  int cost;
  int time;
};

bool operator==(const spp_spptw_res_cont& res_cont_1, const spp_spptw_res_cont& res_cont_2) {
  return (res_cont_1.cost == res_cont_2.cost && res_cont_1.time == res_cont_2.time);
}

bool operator<(const spp_spptw_res_cont& res_cont_1, const spp_spptw_res_cont& res_cont_2) {
  if (res_cont_1.cost > res_cont_2.cost) {
    return false;
  }
  if (res_cont_1.cost == res_cont_2.cost) {
    return res_cont_1.time < res_cont_2.time;
  }
  return true;
}

// ResourceExtensionFunction model
struct ref_spptw {
  bool operator()(const Graph& g, spp_spptw_res_cont& new_cont, const spp_spptw_res_cont& old_cont,
                  graph_edge_descriptor_t<Graph> ed) const {
    const ArcProp& arc_prop = g[ed];
    const VertProp& vert_prop = g[target(ed, g)];
    new_cont.cost = old_cont.cost + arc_prop.cost;
    int& i_time = new_cont.time;
    i_time = old_cont.time + arc_prop.time;
    i_time < vert_prop.eat ? i_time = vert_prop.eat : 0;
    return i_time <= vert_prop.lat ? true : false;
  }
};

// DominanceFunction model
struct dominance_spptw {
  bool operator()(const spp_spptw_res_cont& res_cont_1, const spp_spptw_res_cont& res_cont_2) const {
    // must be "<=" here!!!
    // must NOT be "<"!!!
    return res_cont_1.cost <= res_cont_2.cost && res_cont_1.time <= res_cont_2.time;
    // this is not a contradiction to the documentation
    // the documentation says:
    // "A label $l_1$ dominates a label $l_2$ if and only if both are
    // resident at the same vertex, and if, for each resource, the resource
    // consumption of $l_1$ is less than or equal to the resource
    // consumption of $l_2$, and if there is at least one resource where
    // $l_1$ has a lower resource consumption than $l_2$." one can think of
    // a new label with a resource consumption equal to that of an old label
    // as being dominated by that old label, because the new one will have a
    // higher number and is created at a later point in time, so one can
    // implicitly use the number or the creation time as a resource for
    // tie-breaking
  }
};
// end data structures for shortest path problem with time windows (spptw)

struct spp_spptw_marked_res_cont {
  spp_spptw_marked_res_cont(Graph::vertex_descriptor v, int c = 0, int t = 0) : cost(c), time(t), marked() {
    marked.insert(v);
  }
  int cost;
  int time;
  std::set<Graph::vertex_descriptor> marked;
};

struct ref_spptw_marked {
  bool operator()(const Graph& g, spp_spptw_marked_res_cont& new_cont, const spp_spptw_marked_res_cont& old_cont,
                  graph_edge_descriptor_t<Graph> ed) const {
    const graph_vertex_descriptor_t<Graph> dest = target(ed, g);

    if (old_cont.marked.find(dest) != old_cont.marked.end()) {
      return false;
    }

    const ArcProp& arc_prop = g[ed];
    const VertProp& vert_prop = g[dest];
    new_cont.cost = old_cont.cost + arc_prop.cost;
    new_cont.marked = old_cont.marked;
    new_cont.marked.insert(dest);
    int& i_time = new_cont.time;
    i_time = old_cont.time + arc_prop.time;
    i_time < vert_prop.eat ? i_time = vert_prop.eat : 0;
    return i_time <= vert_prop.lat;
  }
};

struct dominance_spptw_marked {
  bool operator()(const spp_spptw_marked_res_cont& res_cont_1, const spp_spptw_marked_res_cont& res_cont_2) const {
    return res_cont_1.time <= res_cont_2.time && res_cont_1.cost <= res_cont_2.cost &&
           std::includes(res_cont_1.marked.begin(), res_cont_1.marked.end(), res_cont_2.marked.begin(),
                         res_cont_2.marked.end());
  }
};

TEST(RCShortestPathsTest, Large) {
  Graph g;
  add_vertex(g, 0, 0, 1000000000);
  add_vertex(g, 1, 56, 142);
  add_vertex(g, 2, 0, 1000000000);
  add_vertex(g, 3, 89, 178);
  add_vertex(g, 4, 0, 1000000000);
  add_vertex(g, 5, 49, 76);
  add_vertex(g, 6, 0, 1000000000);
  add_vertex(g, 7, 98, 160);
  add_vertex(g, 8, 0, 1000000000);
  add_vertex(g, 9, 90, 158);
  add_edge(0, 7, g, 6, 33, 2);
  add_edge(0, 6, g, 5, 31, 6);
  add_edge(0, 4, g, 3, 14, 4);
  add_edge(0, 1, g, 0, 43, 8);
  add_edge(0, 4, g, 4, 28, 10);
  add_edge(0, 3, g, 1, 31, 10);
  add_edge(0, 3, g, 2, 1, 7);
  add_edge(0, 9, g, 7, 25, 9);
  add_edge(1, 0, g, 8, 37, 4);
  add_edge(1, 6, g, 9, 7, 3);
  add_edge(2, 6, g, 12, 6, 7);
  add_edge(2, 3, g, 10, 13, 7);
  add_edge(2, 3, g, 11, 49, 9);
  add_edge(2, 8, g, 13, 47, 5);
  add_edge(3, 4, g, 17, 5, 10);
  add_edge(3, 1, g, 15, 47, 1);
  add_edge(3, 2, g, 16, 26, 9);
  add_edge(3, 9, g, 21, 24, 10);
  add_edge(3, 7, g, 20, 50, 10);
  add_edge(3, 0, g, 14, 41, 4);
  add_edge(3, 6, g, 19, 6, 1);
  add_edge(3, 4, g, 18, 8, 1);
  add_edge(4, 5, g, 26, 38, 4);
  add_edge(4, 9, g, 27, 32, 10);
  add_edge(4, 3, g, 24, 40, 3);
  add_edge(4, 0, g, 22, 7, 3);
  add_edge(4, 3, g, 25, 28, 9);
  add_edge(4, 2, g, 23, 39, 6);
  add_edge(5, 8, g, 32, 6, 2);
  add_edge(5, 2, g, 30, 26, 10);
  add_edge(5, 0, g, 28, 38, 9);
  add_edge(5, 2, g, 31, 48, 10);
  add_edge(5, 9, g, 33, 49, 2);
  add_edge(5, 1, g, 29, 22, 7);
  add_edge(6, 1, g, 34, 15, 7);
  add_edge(6, 7, g, 35, 20, 3);
  add_edge(7, 9, g, 40, 1, 3);
  add_edge(7, 0, g, 36, 23, 5);
  add_edge(7, 6, g, 38, 36, 2);
  add_edge(7, 6, g, 39, 18, 10);
  add_edge(7, 2, g, 37, 2, 1);
  add_edge(8, 5, g, 46, 36, 5);
  add_edge(8, 1, g, 42, 13, 10);
  add_edge(8, 0, g, 41, 40, 5);
  add_edge(8, 1, g, 43, 32, 8);
  add_edge(8, 6, g, 47, 25, 1);
  add_edge(8, 2, g, 44, 44, 3);
  add_edge(8, 3, g, 45, 11, 9);
  add_edge(9, 0, g, 48, 41, 5);
  add_edge(9, 1, g, 49, 44, 7);

  // spp without resource constraints

  std::vector<std::vector<graph_edge_descriptor_t<Graph> > > opt_solutions;
  std::vector<spp_no_rc_res_cont> pareto_opt_rcs_no_rc;
  std::vector<int> i_vec_opt_solutions_spp_no_rc;
  for (int s = 0; s < 10; ++s) {
    for (int t = 0; t < 10; ++t) {
      r_c_shortest_paths(g, get(&VertProp::num, g), get(&ArcProp::num, g), s, t, opt_solutions, pareto_opt_rcs_no_rc,
                         spp_no_rc_res_cont(0), ref_no_res_cont(), dominance_no_res_cont(),
                         default_r_c_shortest_paths_visitor());
      i_vec_opt_solutions_spp_no_rc.push_back(pareto_opt_rcs_no_rc[0].cost);
    }
  }

  // std::vector<graph_traits<Graph>::vertex_descriptor>
  //  p( num_vertices( g ) );
  // std::vector<int> d( num_vertices( g ) );
  // std::vector<int> i_vec_dijkstra_distances;
  // std::cout << "Dijkstra:" << std::endl;
  // for( int s = 0; s < 10; ++s )
  //{
  //  dijkstra_shortest_paths( g,
  //                           s,
  //                           &p[0],
  //                           &d[0],
  //                           get( &ArcProp::cost, g ),
  //                           get( &VertProp::num, g ),
  //                           std::less<int>(),
  //                           closed_plus<int>(),
  //                           (std::numeric_limits<int>::max)(),
  //                           0,
  //                           default_dijkstra_visitor() );
  //  for( int t = 0; t < 10; ++t )
  //  {
  //    i_vec_dijkstra_distances.push_back( d[t] );
  //    std::cout << "From " << s << " to " << t << ": " << d[t] << std::endl;
  //  }
  //}

  std::vector<int> i_vec_correct_solutions = {
      0,  22, 27, 1,  6,  44, 7,  27, 50, 25, 37, 0,  29, 38, 43, 81, 7,  27, 76, 28, 25, 21, 0,  13, 18,
      56, 6,  26, 47, 27, 12, 21, 26, 0,  5,  43, 6,  26, 49, 24, 7,  29, 34, 8,  0,  38, 14, 34, 44, 32,
      29, 19, 26, 17, 22, 0,  23, 43, 6,  41, 43, 15, 22, 35, 40, 78, 0,  20, 69, 21, 23, 23, 2,  15, 20,
      58, 8,  0,  49, 1,  23, 13, 37, 11, 16, 36, 17, 37, 0,  35, 41, 44, 68, 42, 47, 85, 48, 68, 91, 0};
  EXPECT_EQ(i_vec_opt_solutions_spp_no_rc.size(), i_vec_correct_solutions.size());
  for (int i = 0; i < i_vec_correct_solutions.size(); ++i) {
    EXPECT_EQ(i_vec_opt_solutions_spp_no_rc[i], i_vec_correct_solutions[i]);
  }

  // spptw
  std::vector<std::vector<graph_edge_descriptor_t<Graph> > > opt_solutions_spptw;
  std::vector<spp_spptw_res_cont> pareto_opt_rcs_spptw;
  std::vector<std::vector<std::vector<std::vector<graph_edge_descriptor_t<Graph> > > > >
      vec_vec_vec_vec_opt_solutions_spptw(10);

  for (int s = 0; s < 10; ++s) {
    for (int t = 0; t < 10; ++t) {
      r_c_shortest_paths(
          g, get(&VertProp::num, g), get(&ArcProp::num, g), s, t, opt_solutions_spptw, pareto_opt_rcs_spptw,
          // be careful, do not simply take 0 as initial value for time
          spp_spptw_res_cont(0, g[s].eat), ref_spptw(), dominance_spptw(), default_r_c_shortest_paths_visitor());
      vec_vec_vec_vec_opt_solutions_spptw[s].push_back(opt_solutions_spptw);
      if (!opt_solutions_spptw.empty()) {
        bool b_is_a_path_at_all = false;
        bool b_feasible = false;
        bool b_correctly_extended = false;
        spp_spptw_res_cont actual_final_resource_levels(0, 0);
        graph_edge_descriptor_t<Graph> ed_last_extended_arc;
        check_r_c_path(g, opt_solutions_spptw[0], spp_spptw_res_cont(0, g[s].eat), true, pareto_opt_rcs_spptw[0],
                       actual_final_resource_levels, ref_spptw(), b_is_a_path_at_all, b_feasible, b_correctly_extended,
                       ed_last_extended_arc);
        EXPECT_TRUE(b_is_a_path_at_all);
        EXPECT_TRUE(b_feasible);
        EXPECT_TRUE(b_correctly_extended);
        b_is_a_path_at_all = false;
        b_feasible = false;
        b_correctly_extended = false;
        spp_spptw_res_cont actual_final_resource_levels2(0, 0);
        graph_edge_descriptor_t<Graph> ed_last_extended_arc2;
        check_r_c_path(g, opt_solutions_spptw[0], spp_spptw_res_cont(0, g[s].eat), false, pareto_opt_rcs_spptw[0],
                       actual_final_resource_levels2, ref_spptw(), b_is_a_path_at_all, b_feasible, b_correctly_extended,
                       ed_last_extended_arc2);
        EXPECT_TRUE(b_is_a_path_at_all);
        EXPECT_TRUE(b_feasible);
        EXPECT_TRUE(b_correctly_extended);
      }
    }
  }

  std::vector<int> i_vec_correct_num_solutions_spptw = {
      1, 2, 3, 1, 3, 1, 2, 1, 2, 1, 1, 1, 4, 1, 3, 1, 1, 1, 2, 2, 4, 1, 1, 1, 4, 1, 2, 1, 1, 3, 2, 2, 2, 1,
      2, 0, 1, 1, 2, 1, 1, 2, 2, 1, 1, 1, 2, 1, 2, 1, 4, 2, 2, 1, 4, 1, 4, 1, 1, 2, 2, 1, 4, 2, 5, 1, 1, 1,
      2, 2, 1, 3, 1, 1, 2, 0, 2, 1, 1, 1, 3, 1, 2, 1, 3, 1, 3, 1, 1, 2, 1, 1, 4, 1, 3, 0, 2, 3, 4, 1};

  for (int s = 0; s < 10; ++s) {
    for (int t = 0; t < 10; ++t) {
      EXPECT_EQ(static_cast<int>(vec_vec_vec_vec_opt_solutions_spptw[s][t].size()),
                i_vec_correct_num_solutions_spptw[10 * s + t]);
    }
  }
}

TEST(RCShortestPathsTest, Small) {
  // one pareto-optimal solution
  Graph g;
  add_vertex(g, 0, 0, 1000000000);
  add_vertex(g, 1, 0, 1000000000);
  add_vertex(g, 2, 0, 1000000000);
  add_vertex(g, 3, 0, 1000000000);
  add_edge(0, 1, g, 0, 1, 1);
  add_edge(0, 2, g, 1, 2, 1);
  add_edge(1, 3, g, 2, 3, 1);
  add_edge(2, 3, g, 3, 1, 1);
  std::vector<graph_edge_descriptor_t<Graph> > opt_solution;
  spp_spptw_res_cont pareto_opt_rc;

  auto g_source = *vertices(g).begin();
  auto g_target = *std::next(vertices(g).begin(), 3);
  r_c_shortest_paths(g, get(&VertProp::num, g), get(&ArcProp::num, g), g_source, g_target, opt_solution, pareto_opt_rc,
                     spp_spptw_res_cont(0, 0), ref_spptw(), dominance_spptw(), default_r_c_shortest_paths_visitor());

  EXPECT_EQ(pareto_opt_rc.cost, 3);
}

TEST(RCShortestPathsTest, Medium) {
  Graph g;
  add_vertex(g, 0, 0, 1000);
  add_vertex(g, 1, 0, 1000);
  add_vertex(g, 2, 0, 974);
  add_vertex(g, 3, 0, 972);
  add_vertex(g, 4, 0, 967);
  add_vertex(g, 5, 678, 801);
  add_edge(0, 2, g, 0, 0, 16);
  add_edge(0, 3, g, 1, 0, 18);
  add_edge(0, 4, g, 2, 0, 23);
  add_edge(0, 5, g, 3, 0, 25);
  add_edge(2, 3, g, 4, 0, 33);
  add_edge(2, 4, g, 5, 0, 15);
  add_edge(2, 5, g, 6, 0, 33);
  add_edge(2, 1, g, 7, 0, 16);
  add_edge(3, 2, g, 8, 0, 33);
  add_edge(3, 4, g, 9, 0, 35);
  add_edge(3, 5, g, 10, 0, 21);
  add_edge(3, 1, g, 11, 0, 18);
  add_edge(4, 2, g, 12, 0, 15);
  add_edge(4, 3, g, 13, 0, 35);
  add_edge(4, 5, g, 14, 0, 25);
  add_edge(4, 1, g, 15, 0, 23);
  add_edge(5, 2, g, 16, 0, 33);
  add_edge(5, 3, g, 17, 0, 21);
  add_edge(5, 4, g, 18, 0, 25);
  add_edge(5, 1, g, 19, 0, 25);

  std::vector<std::vector<graph_edge_descriptor_t<Graph> > > pareto_opt_marked_solutions;
  std::vector<spp_spptw_marked_res_cont> pareto_opt_marked_resource_containers;

  auto g_source = *vertices(g).begin();
  auto g_target = *std::next(vertices(g).begin());
  r_c_shortest_paths(g, get(&VertProp::num, g), get(&ArcProp::num, g), g_source, g_target, pareto_opt_marked_solutions,
                     pareto_opt_marked_resource_containers, spp_spptw_marked_res_cont(0, 0, 0), ref_spptw_marked(),
                     dominance_spptw_marked(), default_r_c_shortest_paths_visitor());

  EXPECT_FALSE(pareto_opt_marked_solutions.empty());
  for (const auto& path : pareto_opt_marked_solutions) {
    EXPECT_FALSE(path.empty());
    const auto front = path.front();
    EXPECT_EQ(target(front, g), g_target);

    auto prev_edge = front;
    for (auto edge_it = std::next(path.begin()), edge_it_end = path.end(); edge_it != edge_it_end; ++edge_it) {
      auto edge = *edge_it;
      auto prev_end = source(prev_edge, g);
      auto current_end = target(edge, g);
      EXPECT_EQ(prev_end, current_end);
      prev_edge = edge;
    }

    const auto back = path.back();
    EXPECT_EQ(source(back, g), g_source);
  }
}

}  // namespace
}  // namespace bagl
