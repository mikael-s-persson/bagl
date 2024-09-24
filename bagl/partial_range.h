// Copyright 2013-2024 Mikael Persson.

#ifndef BAGL_BAGL_PARTIAL_RANGE_H_
#define BAGL_BAGL_PARTIAL_RANGE_H_

#include <iterator>
#include <ranges>
#include <utility>
#include <memory>

#include "bagl/has_trait_member.h"

namespace bagl {
namespace partial_range_detail {
// Working around missing std::ranges::sized_range.
BAGL_GRAPH_HAS_MEMBER_FUNCTION(size);
} // partial_range_detail

// This is a partial view over an underlying range which ensures that an iterator moving
// within the range continues to be correctly tied to its parent range, i.e., not dangling.
template <typename BaseRange>
class partial_view : public std::ranges::view_interface<partial_view<BaseRange>> {
 public:
  static constexpr bool is_base_borrowable = std::ranges::borrowed_range<BaseRange>;
  using iterator = std::ranges::iterator_t<BaseRange>;
  using sentinel = std::ranges::sentinel_t<BaseRange>;

  template <typename... Args>
  explicit partial_view(Args&&... args) : base_range_(make_base(std::forward<Args>(args)...)), current_begin_(get_base().begin()) {}

  explicit partial_view(BaseRange&& base_range) : base_range_(make_base(std::move(base_range))), current_begin_(get_base().begin()) {}
  explicit partial_view(const BaseRange& base_range) : base_range_(make_base(base_range)), current_begin_(get_base().begin()) {}

  iterator begin() const { return current_begin_; }
  sentinel end() const { return get_base().end(); }

  // Move the current begin iterator for this partial view.
  void move_begin_to(iterator new_begin) { current_begin_ = new_begin; }

  // Move the current begin iterator to the end, i.e., make it empty.
  // This is needed in case 
  void move_begin_to_end() {
    if constexpr (std::is_same_v<iterator, sentinel>) {
      current_begin_ = end();
    } else if constexpr (partial_range_detail::has_size_v<BaseRange>) {
      current_begin_ = base_begin();
      std::advance(current_begin_, get_base().size());
    } else {
      // Last resort.
      while (current_begin_ != end()) {
        ++current_begin_;
      }
    }
  }

  // Get the actual begin iterator of the underlying range.
  iterator base_begin() const { return get_base().begin(); }

 private:
  using Storage = std::conditional_t<is_base_borrowable, BaseRange, std::unique_ptr<BaseRange>>;
  Storage base_range_;
  iterator current_begin_;

  const BaseRange& get_base() const {
    if constexpr (is_base_borrowable) {
      return base_range_;
    } else {
      return *base_range_;
    }
  }
  template <typename... Args>
  static Storage make_base(Args&&... args) {
    if constexpr (is_base_borrowable) {
      return BaseRange{std::forward<Args>(args)...};
    } else {
      return std::make_unique<BaseRange>(std::forward<Args>(args)...);
    }
  }
};

}  // namespace bagl

#endif  // BAGL_BAGL_PARTIAL_RANGE_H_
