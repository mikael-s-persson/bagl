// Copyright Douglas Gregor 2004.
// Copyright 2024 Mikael Persson - Modernized to C++20

#include <algorithm>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include "bagl/adjacency_list.h"
#include "bagl/adjacency_matrix.h"
#include "bagl/filtered_graph.h"
#include "bagl/graph_traits.h"
#include "bagl/graph_utility.h"
#include "bagl/subgraph.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"

namespace bagl {
namespace {

struct City {
  City() {}
  City(const std::string& a_name, int a_pop, int a_zipcode) : name(a_name), population(a_pop) {
    zipcodes.push_back(a_zipcode);
  }

  std::string name;
  int population;
  std::vector<int> zipcodes;
};

std::ostream& operator<<(std::ostream& out, const City& city) {
  out << city.name << ' ' << city.population << ' ';
  std::copy(city.zipcodes.begin(), city.zipcodes.end(), std::ostream_iterator<int>(out, " "));
  out << -1;
  return out;
}

bool operator==(const City& c1, const City& c2) {
  return (c1.name == c2.name && c1.population == c2.population && c1.zipcodes == c2.zipcodes);
}

struct Highway {
  Highway() {}
  Highway(const std::string& a_name, double a_miles, int a_speed_limit = 65, int a_lanes = 4, bool a_divided = true)
      : name(a_name), miles(a_miles), speed_limit(a_speed_limit), lanes(a_lanes), divided(a_divided) {}

  std::string name;
  double miles;
  int speed_limit;
  int lanes;
  bool divided;
};

std::ostream& operator<<(std::ostream& out, const Highway& highway) {
  return out << highway.name << ' ' << highway.miles << ' ' << highway.miles << ' ' << highway.speed_limit << ' '
             << highway.lanes << ' ' << highway.divided;
}

bool operator==(const Highway& h1, const Highway& h2) {
  return (h1.name == h2.name && h1.miles == h2.miles && h1.speed_limit == h2.speed_limit && h1.lanes == h2.lanes &&
          h1.divided == h2.divided);
}

template <typename T>
class BundledPropertiesTest : public ::testing::Test {};
TYPED_TEST_SUITE_P(BundledPropertiesTest);

TYPED_TEST_P(BundledPropertiesTest, Basic) {
  using Graph = TypeParam;
  using Vertex = graph_vertex_descriptor_t<Graph>;
  constexpr bool CanAddVertex = concepts::VertexMutableGraph<Graph>;

  Graph map(CanAddVertex ? 2 : 3);

  auto v_rg = vertices(map);
  auto vi = v_rg.begin();
  Vertex v = *vi++;
  map[v].name = "Troy";
  map[v].population = 49170;
  map[v].zipcodes.push_back(12180);

  Vertex u = *vi++;
  map[u].name = "Albany";
  map[u].population = 95658;
  map[u].zipcodes.push_back(12201);

  // Try adding a vertex with a property value
  Vertex bloomington = u;
  if constexpr (CanAddVertex) {
    bloomington = add_vertex(map, City("Bloomington", 39000, 47401));
  } else {
    bloomington = *vi++;
    get(vertex_bundle, map)[bloomington] = City("Bloomington", 39000, 47401);
  }
  EXPECT_EQ(get(vertex_bundle, map)[bloomington].zipcodes[0], 47401);

  auto [e, e_added] = add_edge(v, u, map);
  ASSERT_TRUE(e_added);
  map[e].name = "I-87";
  map[e].miles = 10;
  map[e].speed_limit = 65;
  map[e].lanes = 4;
  map[e].divided = true;

  auto [our_trip, out_trip_added] = add_edge(v, bloomington, map, Highway("Long", 1000));
  ASSERT_TRUE(out_trip_added);
  EXPECT_EQ(get(edge_bundle, map, our_trip).miles, 1000);

  EXPECT_EQ(get(get(&City::name, map), v), "Troy");
  EXPECT_EQ(get(get(&Highway::name, map), e), "I-87");
  EXPECT_EQ(get(&City::name, map, u), "Albany");
  EXPECT_EQ(get(&Highway::name, map, e), "I-87");
  put(&City::population, map, v, 49168);
  EXPECT_EQ(get(&City::population, map)[v], 49168);

  filtered_graph<Graph, keep_all> fmap(map, keep_all());
  EXPECT_EQ(get(edge_bundle, map, our_trip).miles, 1000);

  EXPECT_EQ(get(get(&City::name, fmap), v), "Troy");
  EXPECT_EQ(get(get(&Highway::name, fmap), e), "I-87");
  EXPECT_EQ(get(&City::name, fmap, u), "Albany");
  EXPECT_EQ(get(&Highway::name, fmap, e), "I-87");
  put(&City::population, fmap, v, 49169);
  EXPECT_EQ(get(&City::population, fmap)[v], 49169);

  std::stringstream out;
  print_graph(map, get(vertex_bundle, map), get(edge_bundle, map), out);
  EXPECT_THAT(
      out.str(),
      ::testing::ContainsRegex(
          R"(Troy\s*49169\s*12180\s*-1\s*-->\s*Albany\s*95658\s*12201\s*-1\s*Bloomington\s*39000\s*47401\s*-1)"));
  EXPECT_THAT(
      out.str(),
      ::testing::ContainsRegex(
          R"(I-87\s*10\s*10\s*65\s*4\s*1\s*\(\s*Troy\s*49169\s*12180\s*-1\s*-->\s*Albany\s*95658\s*12201\s*-1\s*\))"));
}

REGISTER_TYPED_TEST_SUITE_P(BundledPropertiesTest, Basic);

using TestAdjList = adjacency_list<list_s, vec_s, bidirectional_s, City, Highway>;
using TestAdjMatrix = adjacency_matrix<directed_s, City, Highway>;
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyListBundledProperties, BundledPropertiesTest, TestAdjList);
INSTANTIATE_TYPED_TEST_SUITE_P(AdjacencyMatrixBundledProperties, BundledPropertiesTest, TestAdjMatrix);

TEST(SubgraphBundledPropertiesTest, Basic) {
  typedef subgraph<adjacency_list<vec_s, vec_s, bidirectional_s, City, property<edge_index_t, int, Highway> > > SubMap;
  typedef graph_vertex_descriptor_t<SubMap> Vertex;

  SubMap map(3);
  auto v_rg = vertices(map);
  auto vi = v_rg.begin();
  Vertex troy = *vi++;
  map[troy].name = "Troy";
  map[*vi++].name = "Bloomington";
  map[*vi++].name = "Endicott";

  SubMap& g1 = map.create_subgraph();
  Vertex troy1 = add_vertex(troy, g1);
  EXPECT_EQ(map[troy].name, g1[troy1].name);
}

}  // namespace
}  // namespace bagl
