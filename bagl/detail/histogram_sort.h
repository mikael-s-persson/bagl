// Copyright 2009 The Trustees of Indiana University.
//  Authors: Jeremiah Willcock, Andrew Lumsdaine

#ifndef BAGL_BAGL_DETAIL_HISTOGRAM_SORT_H_
#define BAGL_BAGL_DETAIL_HISTOGRAM_SORT_H_

#include <cassert>
#include <iterator>
#include <ranges>
#include <vector>

#include "bagl/zip_range.h"

namespace bagl::histogram_detail {

template <std::ranges::range Range>
std::size_t reserve_count_for_single_pass(Range rg) {
  if constexpr (std::ranges::sized_range<Range>) {
    return rg.size();
  } else {
    return 0;
  }
}

template <std::ranges::input_range KeyRange, std::ranges::random_access_range RowstartRange, typename KeyFilter,
          typename KeyTransform>
requires std::ranges::output_range<RowstartRange, std::ranges::range_value_t<RowstartRange>>
void count_starts(KeyRange keys, RowstartRange starts, std::size_t numkeys, KeyFilter key_filter,
                  KeyTransform key_transform) {
  assert(starts.size() >= numkeys + 1);
  // Put the degree of each vertex v into m_rowstart[v + 1]
  for (auto k : keys) {
    if (key_filter(k)) {
      assert(key_transform(k) < numkeys);
      ++starts[key_transform(k) + 1];
    }
  }

  // Compute the partial sum of the degrees to get the actual values
  // of m_rowstart
  std::ranges::range_value_t<RowstartRange> start_of_this_row = 0;
  starts[0] = start_of_this_row;
  for (std::size_t i = 1; i < numkeys + 1; ++i) {
    start_of_this_row += starts[i];
    starts[i] = start_of_this_row;
  }
}

template <std::ranges::input_range KeyRange, std::ranges::random_access_range RowstartRange,
          std::ranges::input_range InputRange, std::ranges::random_access_range OutputRange, typename KeyFilter,
          typename KeyTransform>
requires std::ranges::output_range<OutputRange, std::ranges::range_value_t<InputRange>>
void histogram_sort(KeyRange keys, RowstartRange rowstarts, std::size_t numkeys, InputRange values_in,
                    OutputRange values_out, KeyFilter key_filter, KeyTransform key_transform) {
  assert(rowstarts.size() >= numkeys + 1);
  using EdgeIndex = std::ranges::range_value_t<RowstartRange>;

  // Histogram sort the edges by their source vertices, putting the
  // targets into m_column.  The index current_insert_positions[v]
  // contains the next location to insert out edges for vertex v.
  std::vector<EdgeIndex> current_insert_positions(rowstarts.begin(), rowstarts.begin() + numkeys);
  auto vit = values_in.begin();
  for (auto k : keys) {
    if (key_filter(k)) {
      std::size_t source = key_transform(k);
      assert(source < numkeys);
      EdgeIndex insert_pos = current_insert_positions[source];
      ++current_insert_positions[source];
      values_out[insert_pos] = *vit;
    }
    ++vit;
  }
}

template <std::ranges::input_range KeyRange, std::ranges::random_access_range RowstartRange,
          std::ranges::input_range InputRange1, std::ranges::random_access_range OutputRange1,
          std::ranges::input_range InputRange2, std::ranges::random_access_range OutputRange2, typename KeyFilter,
          typename KeyTransform>
requires std::ranges::output_range<OutputRange1, std::ranges::range_value_t<InputRange1>> &&
    std::ranges::output_range<OutputRange2, std::ranges::range_value_t<InputRange2>>
void histogram_sort(KeyRange keys, RowstartRange rowstarts, std::size_t numkeys, InputRange1 values1_in,
                    OutputRange1 values1_out, InputRange2 values2_in, OutputRange2 values2_out, KeyFilter key_filter,
                    KeyTransform key_transform) {
  histogram_sort(keys, rowstarts, numkeys, zip_range(values1_in, values2_in), zip_range(values1_out, values2_out),
                 key_filter, key_transform);
}

template <std::ranges::random_access_range KeyRange, std::ranges::random_access_range RowstartRange,
          std::ranges::random_access_range OutputRange, typename KeyTransform>
requires std::ranges::output_range<OutputRange, std::ranges::range_value_t<OutputRange>> &&
    std::ranges::output_range<KeyRange, std::ranges::range_value_t<KeyRange>>
void histogram_sort_inplace(KeyRange keys, RowstartRange rowstarts, std::size_t numkeys, OutputRange values,
                            KeyTransform key_transform) {
  assert(rowstarts.size() >= numkeys + 1);
  using EdgeIndex = std::ranges::range_value_t<RowstartRange>;

  // 1. Copy m_rowstart (except last element) to get insert positions
  std::vector<EdgeIndex> insert_positions(rowstarts.begin(), rowstarts.begin() + numkeys);
  // 2. Swap the sources and targets into place
  for (std::size_t i = 0; i < rowstarts[numkeys]; ++i) {
    assert(key_transform(keys[i]) < numkeys);
    // While edge i is not in the right bucket:
    while (!(i >= rowstarts[key_transform(keys[i])] && i < insert_positions[key_transform(keys[i])])) {
      // Add a slot in the right bucket
      std::size_t target_pos = insert_positions[key_transform(keys[i])]++;
      assert(target_pos < rowstarts[key_transform(keys[i]) + 1]);
      if (target_pos == i) {
        continue;
      }
      // Swap this edge into place
      using std::swap;
      swap(keys[i], keys[target_pos]);
      swap(values[i], values[target_pos]);
    }
  }
}

template <std::ranges::random_access_range KeyRange, std::ranges::random_access_range RowstartRange,
          std::ranges::random_access_range OutputRange1, std::ranges::random_access_range OutputRange2,
          typename KeyTransform>
requires std::ranges::output_range<OutputRange1, std::ranges::range_value_t<OutputRange1>> &&
    std::ranges::output_range<OutputRange2, std::ranges::range_value_t<OutputRange2>> &&
    std::ranges::output_range<KeyRange, std::ranges::range_value_t<KeyRange>>
void histogram_sort_inplace(KeyRange keys, RowstartRange rowstarts, std::size_t numkeys, OutputRange1 values1,
                            OutputRange2 values2, KeyTransform key_transform) {
  histogram_sort_inplace(keys, rowstarts, numkeys, zip_range(values1, values2), key_transform);
}

template <std::ranges::input_range InputRange, typename VerticesSize>
void split_into_separate_coords(InputRange r, std::vector<VerticesSize>& firsts, std::vector<VerticesSize>& seconds) {
  firsts.clear();
  seconds.clear();
  std::size_t reserve_size = reserve_count_for_single_pass(r);
  firsts.reserve(reserve_size);
  seconds.reserve(reserve_size);
  for (const auto& [f, s] : r) {
    firsts.push_back(f);
    seconds.push_back(s);
  }
}

template <std::ranges::input_range InputRange, typename VerticesSize, typename SourceFilter>
void split_into_separate_coords_filtered(InputRange r, std::vector<VerticesSize>& firsts,
                                         std::vector<VerticesSize>& seconds, SourceFilter&& filter) {
  firsts.clear();
  seconds.clear();
  for (const auto& [f, s] : r) {
    if (filter(f)) {
      firsts.push_back(f);
      seconds.push_back(s);
    }
  }
}

template <std::ranges::input_range InputRange, std::ranges::input_range PropInputRange, typename VerticesSize,
          typename PropType, typename SourceFilter>
void split_into_separate_coords_filtered(InputRange r, PropInputRange props, std::vector<VerticesSize>& firsts,
                                         std::vector<VerticesSize>& seconds, std::vector<PropType>& props_out,
                                         const SourceFilter& filter) {
  firsts.clear();
  seconds.clear();
  props_out.clear();
  auto props_it = props.begin();
  for (const auto& [f, s] : r) {
    if (filter(f)) {
      firsts.push_back(f);
      seconds.push_back(s);
      props_out.push_back(*props_it);
    }
    ++props_it;
  }
}

}  // namespace bagl::histogram_detail

#endif  // BAGL_BAGL_DETAIL_HISTOGRAM_SORT_H_
