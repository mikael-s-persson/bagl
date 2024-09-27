
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
#include "bagl/boykov_kolmogorov_max_flow.h"
#include "bagl/breadth_first_search.h"
#include "bagl/bron_kerbosch_all_cliques.h"
#include "bagl/bst_inorder_iterator.h"
#include "bagl/buffer_concepts.h"
#include "bagl/chrobak_payne_drawing.h"
#include "bagl/constant_property_map.h"
#include "bagl/container_traits.h"
#include "bagl/core_numbers.h"
#include "bagl/d_ary_heap.h"
#include "bagl/depth_first_search.h"
#include "bagl/detail/adjlist_containers.h"
#include "bagl/detail/adjlist_ranges.h"
#include "bagl/detail/array_binary_tree.h"
#include "bagl/detail/bfl_tree_ranges.h"
#include "bagl/detail/container_generators.h"
#include "bagl/detail/edge.h"
#include "bagl/detail/ltree_containers.h"
#include "bagl/detail/vebl_tree_ranges.h"
#include "bagl/dijkstra_shortest_paths.h"
#include "bagl/exception.h"
#include "bagl/function_property_map.h"
#include "bagl/graph_archetypes.h"
#include "bagl/graph_as_tree.h"
#include "bagl/graph_concepts.h"
#include "bagl/graph_mutability_traits.h"
#include "bagl/graph_selectors.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/has_trait_member.h"
#include "bagl/indirect_cmp.h"
#include "bagl/is_heap.h"
#include "bagl/linked_tree.h"
#include "bagl/lookup_edge.h"
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
#include "bagl/prim_minimum_spanning_tree.h"
#include "bagl/properties.h"
#include "bagl/property.h"
#include "bagl/property_map.h"
#include "bagl/relax.h"
#include "bagl/tiernan_all_cycles.h"
#include "bagl/topology.h"
#include "bagl/transform_iterator.h"
#include "bagl/tree_adaptor.h"
#include "bagl/tree_concepts.h"
#include "bagl/tree_traits.h"
#include "bagl/two_bit_color_map.h"
#include "bagl/undirected_dfs.h"
#include "bagl/vebl_d_ary_tree.h"
#include "bagl/vector_property_map.h"
#include "bagl/visitors.h"

int dummy() { return 0; }
