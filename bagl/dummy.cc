
#include "bagl/adjacency_list.h"
#include "bagl/adjacency_range.h"
#include "bagl/astar_search.h"
#include "bagl/bandwidth.h"
#include "bagl/bc_clustering.h"
#include "bagl/bellman_ford_shortest_paths.h"
#include "bagl/betweenness_centrality.h"
#include "bagl/bfl_d_ary_tree.h"
#include "bagl/biconnected_components.h"
#include "bagl/bipartite.h"
#include "bagl/boyer_myrvold_planar_test.h"
#include "bagl/boykov_kolmogorov_max_flow.h"
#include "bagl/breadth_first_search.h"
#include "bagl/bron_kerbosch_all_cliques.h"
#include "bagl/bst_inorder_iterator.h"
#include "bagl/buffer_concepts.h"
#include "bagl/chrobak_payne_drawing.h"
#include "bagl/circle_layout.h"
#include "bagl/closeness_centrality.h"
#include "bagl/clustering_coefficient.h"
#include "bagl/compressed_sparse_row_graph.h"
#include "bagl/connected_components.h"
#include "bagl/constant_property_map.h"
#include "bagl/container_property_map.h"
#include "bagl/container_traits.h"
#include "bagl/copy.h"
#include "bagl/core_numbers.h"
#include "bagl/create_condensation_graph.h"
#include "bagl/cuthill_mckee_ordering.h"
#include "bagl/cycle_canceling.h"
#include "bagl/d_ary_heap.h"
#include "bagl/dag_shortest_paths.h"
#include "bagl/degree_centrality.h"
#include "bagl/depth_first_search.h"
#include "bagl/detail/adjlist_containers.h"
#include "bagl/detail/adjlist_ranges.h"
#include "bagl/detail/array_binary_tree.h"
#include "bagl/detail/bfl_tree_ranges.h"
#include "bagl/detail/compressed_sparse_row_struct.h"
#include "bagl/detail/container_generators.h"
#include "bagl/detail/disjoint_sets.h"
#include "bagl/detail/edge.h"
#include "bagl/detail/histogram_sort.h"
#include "bagl/detail/incremental_components.h"
#include "bagl/detail/index.h"
#include "bagl/detail/indexed_properties.h"
#include "bagl/detail/is_distributed_selector.h"
#include "bagl/detail/ltree_containers.h"
#include "bagl/detail/set_adaptor.h"
#include "bagl/detail/vebl_tree_ranges.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/dijkstra_shortest_paths_no_color_map.h"
#include "bagl/disjoint_sets.h"
#include "bagl/edges_from_out_edges.h"
#include "bagl/edmonds_karp_max_flow.h"
#include "bagl/eccentricity.h"
#include "bagl/edge_coloring.h"
#include "bagl/edge_connectivity.h"
#include "bagl/edge_list.h"
#include "bagl/erdos_renyi_generator.h"
#include "bagl/exception.h"
#include "bagl/exterior_property.h"
#include "bagl/filtered_graph.h"
#include "bagl/find_flow_cost.h"
#include "bagl/floyd_warshall_shortest.h"
#include "bagl/fruchterman_reingold.h"
#include "bagl/function_property_map.h"
#include "bagl/geodesic_distance.h"
#include "bagl/geodesic.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_as_tree.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_mutability_traits.h"
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/gursoy_atun_layout.h"
#include "bagl/has_trait_member.h"
#include "bagl/hawick_circuits.h"
#include "bagl/howard_cycle_ratio.h"
#include "bagl/incremental_components.h"
#include "bagl/indirect_cmp.h"
#include "bagl/is_heap.h"
#include "bagl/is_kuratowski_subgraph.h"
#include "bagl/is_straight_line_drawing.h"
#include "bagl/isomorphism.h"
#include "bagl/johnson_all_pairs_shortest.h"
#include "bagl/kamada_kawai_spring_layout.h"
#include "bagl/king_ordering.h"
#include "bagl/kruskal_min_spanning_tree.h"
#include "bagl/labeled_graph.h"
#include "bagl/linked_tree.h"
#include "bagl/lookup_edge.h"
#include "bagl/loop_erased_random_walk.h"
#include "bagl/make_biconnected_planar.h"
#include "bagl/make_connected.h"
#include "bagl/make_maximal_planar.h"
#include "bagl/matrix_property_map.h"
#include "bagl/maximum_adjacency_search.h"
#include "bagl/more_property_maps.h"
#include "bagl/mutable_heap.h"
#include "bagl/mutable_queue.h"
#include "bagl/named_graph.h"
#include "bagl/neighbor_bfs.h"
#include "bagl/null_property_map.h"
#include "bagl/numeric_values.h"
#include "bagl/one_bit_color_map.h"
#include "bagl/partial_range.h"
#include "bagl/planar_face_traversal.h"
#include "bagl/prim_minimum_spanning_tree.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"
#include "bagl/random.h"
#include "bagl/relax.h"
#include "bagl/reverse_graph.h"
#include "bagl/tiernan_all_cycles.h"
#include "bagl/topological_sort.h"
#include "bagl/topology.h"
#include "bagl/transform_iterator.h"
#include "bagl/transform_value_property_map.h"
#include "bagl/tree_adaptor.h"
#include "bagl/tree_concepts.h"
#include "bagl/tree_traits.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/undirected_dfs.h"
#include "bagl/vebl_d_ary_tree.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"
#include "bagl/wavefront.h"
#include "bagl/zip_range.h"

int dummy() {
    #if 0
  auto rg = bagl::partial_view(std::ranges::iota_view(0, 10));
  auto rg2 = rg | std::views::transform([](int i) { return i * 2; });
  int u = 5;
  using small_graph_t = bagl::adjacency_list<bagl::vec_s, bagl::vec_s, bagl::undirected_s>;
  //auto eo_rg = std::ranges::iota_view(0, 10) | std::views::transform([u](const auto& e_id) { return small_graph_t::bidir_edge_descriptor(u, e_id); });
  small_graph_t g;
  auto v_rg = vertices(g);
  auto eo_rg = g.m_pack.out_edges(0);
  //auto ueo_rg = bagl::adjlist_detail::adjlist_undir_ioerange<small_graph_t::edge_descriptor, decltype(g.m_pack.in_edges(0)), decltype(g.m_pack.out_edges(0))>(true, g.m_pack.in_edges(0), g.m_pack.out_edges(0));
  auto ueo_rg = bagl::adjlist_detail::make_adjlist_undir_ioerange<small_graph_t::edge_descriptor>(true, g.m_pack.in_edges(0), g.m_pack.out_edges(0));
  for (auto e : ueo_rg) {
    (void)e;
  }
  using EdgesFromOutEdges = decltype(g.m_pack.edges());
//  EdgesFromOutEdges e_rg{g.m_pack};
//  auto b = e_rg.begin();
//  auto en = std::ranges::end(e_rg);
  static_assert(std::ranges::input_range<decltype(g.m_pack.edges())>);
//  auto e_rg = std::ranges::transform_view(g.m_pack.edges(), [](const auto& e) { return bagl::container_detail::undir_edge_desc(e); });
//  auto e_rg = bagl::adjlist_detail::adjlist_undir_eiter_range(g.m_pack.edges());
  auto e_rg = bagl::edges_from_out_edges<small_graph_t>(g);
  static_assert(std::ranges::input_range<decltype(edges(g))>);
  #endif
  return 0;
}
