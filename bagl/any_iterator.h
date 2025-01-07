// Copyright 2024 Mikael Persson

#ifndef BAGL_BAGL_ANY_ITERATOR_H_
#define BAGL_BAGL_ANY_ITERATOR_H_

#include <any>
#include <cassert>
#include <iterator>
#include <memory>
#include <ranges>
#include <type_traits>

namespace bagl {

// An iterator type that type-erase the underlying iterator.

// This is a dynamically dispatched classic iterator facade.
template <typename Reference = std::any>
class any_iterator_interface {
 public:
  any_iterator_interface() = default;
  any_iterator_interface(const any_iterator_interface&) = delete;
  any_iterator_interface(any_iterator_interface&&) = delete;
  any_iterator_interface& operator=(const any_iterator_interface&) = delete;
  any_iterator_interface& operator=(any_iterator_interface&&) = delete;
  virtual ~any_iterator_interface() = default;

  [[nodiscard]] virtual std::unique_ptr<any_iterator_interface> clone() const = 0;

  [[nodiscard]] virtual Reference dereference() const = 0;
  virtual void increment() = 0;
  virtual void decrement() = 0;
  virtual void advance(std::ptrdiff_t) = 0;
  [[nodiscard]] virtual bool is_equal(const any_iterator_interface&) const = 0;
  [[nodiscard]] virtual std::ptrdiff_t distance_from(const any_iterator_interface&) const = 0;
};

template <typename Reference = std::any>
class any_iterator {
 public:
  // Random access may be a lie (not O(1)).
  using iterator_category = std::random_access_iterator_tag;
  using value_type = std::remove_reference_t<Reference>;
  using difference_type = std::ptrdiff_t;
  using pointer = std::add_pointer_t<value_type>;
  using reference = Reference;

  any_iterator() = default;
  explicit any_iterator(std::unique_ptr<any_iterator_interface<Reference>> it) : it_(std::move(it)) {}
  any_iterator(const any_iterator& rhs)
      : it_(rhs.it_ == nullptr ? std::unique_ptr<any_iterator_interface<Reference>>() : rhs.it_->clone()) {}
  any_iterator& operator=(const any_iterator& rhs) {
    if (this == &rhs) {
      return *this;
    }
    it_ = rhs.it_ == nullptr ? std::unique_ptr<any_iterator_interface<Reference>>() : rhs.it_->clone();
    return *this;
  }
  any_iterator(any_iterator&&) noexcept = default;
  any_iterator& operator=(any_iterator&&) noexcept = default;
  ~any_iterator() = default;

  reference operator*() const { return it_->dereference(); }
  pointer operator->() const { return &it_->dereference(); }

  any_iterator& operator++() {
    it_->increment();
    return *this;
  }
  any_iterator operator++(int) {
    any_iterator result = *this;
    it_->increment();
    return result;
  }
  any_iterator& operator--() {
    it_->decrement();
    return *this;
  }
  any_iterator operator--(int) {
    any_iterator result = *this;
    it_->decrement();
    return result;
  }

  any_iterator operator+(difference_type n) const {
    any_iterator result = *this;
    result += n;
    return result;
  }
  any_iterator& operator+=(difference_type n) {
    it_->advance(n);
    return *this;
  }
  any_iterator operator-(difference_type n) const {
    any_iterator result = *this;
    result -= n;
    return result;
  }
  any_iterator& operator-=(difference_type n) {
    it_->advance(n);
    return *this;
  }
  value_type operator[](difference_type n) const {
    any_iterator result = *this;
    result += n;
    return *result;
  }
  any_iterator operator-(const any_iterator& rhs) const { return it_->distance_from(*rhs.it_); }

  bool operator==(const any_iterator& rhs) const {
    if (it_ == nullptr && rhs.it_ == nullptr) {
      return true;
    }
    if (it_ == nullptr || rhs.it_ == nullptr) {
      return false;
    }
    return it_->is_equal(*rhs.it_);
  }
  bool operator!=(const any_iterator& rhs) const { return !(*this == rhs); }
  bool operator<(const any_iterator& rhs) const { return (*this - rhs) < 0; }
  bool operator>(const any_iterator& rhs) const { return (*this - rhs) > 0; }
  bool operator<=(const any_iterator& rhs) const { return (*this - rhs) <= 0; }
  bool operator>=(const any_iterator& rhs) const { return (*this - rhs) >= 0; }

 private:
  std::unique_ptr<any_iterator_interface<Reference>> it_;
};

// Just to check.
static_assert(std::input_iterator<any_iterator<std::any>>);

template <std::input_iterator Iter, typename Reference = std::iter_reference_t<Iter>>
requires std::convertible_to<std::iter_reference_t<Iter>, Reference>
class any_iterator_wrapper : public any_iterator_interface<Reference> {
 public:
  using Base = any_iterator_interface<Reference>;
  using Self = any_iterator_wrapper<Iter, Reference>;
  any_iterator_wrapper() = default;
  explicit any_iterator_wrapper(Iter it) : it_(std::move(it)) {}

  [[nodiscard]] std::unique_ptr<Base> clone() const override { return std::make_unique<any_iterator_wrapper>(it_); }

  [[nodiscard]] Reference dereference() const override { return *it_; }
  void increment() override { ++it_; }
  void decrement() override {
    if constexpr (std::bidirectional_iterator<Iter>) {
      --it_;
    } else {
      assert(false && "The underlying iterator is not bidirectional!");
    }
  }
  void advance(std::ptrdiff_t d) override { std::advance(it_, d); }
  [[nodiscard]] bool is_equal(const Base& rhs) const override { return static_cast<const Self&>(rhs).it_ == it_; }
  [[nodiscard]] std::ptrdiff_t distance_from(const Base& rhs) const override {
    return std::distance(static_cast<const Self&>(rhs).it_, it_);
  }

 private:
  Iter it_;
};

template <typename Reference, std::input_iterator Iter>
auto make_any_iterator_to(Iter it) {
  return any_iterator<Reference>{std::make_unique<any_iterator_wrapper<Iter, Reference>>(std::move(it))};
}

template <std::input_iterator Iter>
auto make_any_iterator_to_any(Iter it) {
  return make_any_iterator_to<std::any>(std::move(it));
}

template <std::input_iterator Iter>
auto make_any_iterator(Iter it) {
  return make_any_iterator_to<std::iter_reference_t<Iter>>(std::move(it));
}

template <typename Reference = std::any>
class any_range_interface {
 public:
  any_range_interface() = default;
  any_range_interface(const any_range_interface&) = delete;
  any_range_interface(any_range_interface&&) = delete;
  any_range_interface& operator=(const any_range_interface&) = delete;
  any_range_interface& operator=(any_range_interface&&) = delete;
  virtual ~any_range_interface() = default;

  [[nodiscard]] virtual std::unique_ptr<any_range_interface<Reference>> clone() const = 0;
  [[nodiscard]] virtual any_iterator<Reference> begin() = 0;
  [[nodiscard]] virtual any_iterator<Reference> end() = 0;
};

template <typename Reference = std::any>
class any_range : public std::ranges::view_interface<any_range<Reference>> {
 public:
  any_range() = default;
  explicit any_range(std::unique_ptr<any_range_interface<Reference>> rg) : range_(std::move(rg)) {}
  any_range(const any_range& rhs)
      : range_(rhs.range_ == nullptr ? std::unique_ptr<any_range_interface<Reference>>() : rhs.range_->clone()) {}
  any_range& operator=(const any_range& rhs) {
    if (this == &rhs) {
      return *this;
    }
    range_ = rhs.range_ == nullptr ? std::unique_ptr<any_range_interface<Reference>>() : rhs.range_->clone();
    return *this;
  }
  any_range(any_range&&) noexcept = default;
  any_range& operator=(any_range&&) noexcept = default;
  ~any_range() = default;

  [[nodiscard]] any_iterator<Reference> begin() const { return range_->begin(); }
  [[nodiscard]] any_iterator<Reference> end() const { return range_->end(); }

 private:
  std::unique_ptr<any_range_interface<Reference>> range_;
};

// Check.
static_assert(std::ranges::input_range<any_range<std::any>>);

template <std::ranges::common_range BaseRange, typename Reference = std::ranges::range_reference_t<const BaseRange>>
class any_range_wrapper : public any_range_interface<Reference> {
 public:
  using Self = any_range_wrapper<BaseRange, Reference>;
  using Iter = std::ranges::iterator_t<BaseRange>;
  any_range_wrapper() = default;
  explicit any_range_wrapper(BaseRange rg) : range_(std::move(rg)) {}

  [[nodiscard]] std::unique_ptr<any_range_interface<Reference>> clone() const override {
    return std::make_unique<Self>(range_);
  }
  [[nodiscard]] any_iterator<Reference> begin() override {
    return any_iterator<Reference>{std::make_unique<any_iterator_wrapper<Iter, Reference>>(range_.begin())};
  }
  [[nodiscard]] any_iterator<Reference> end() override {
    return any_iterator<Reference>{std::make_unique<any_iterator_wrapper<Iter, Reference>>(range_.end())};
  }

 private:
  BaseRange range_;
};

template <typename Reference, std::ranges::common_range BaseRange>
auto make_any_range_to(BaseRange&& rg) {
  return any_range<Reference>{std::make_unique<any_range_wrapper<decltype(std::views::all(std::forward<BaseRange>(rg))), Reference>>(
      std::views::all(std::forward<BaseRange>(rg)))};
}

template <std::ranges::common_range BaseRange>
auto make_any_range_to_any(BaseRange&& rg) {
  return make_any_range_to<std::any>(std::forward<BaseRange>(rg));
}

template <std::ranges::common_range BaseRange>
auto make_any_range(BaseRange&& rg) {
  return make_any_range_to<std::ranges::range_reference_t<BaseRange>>(std::forward<BaseRange>(rg));
}

}  // namespace bagl

#endif  // BAGL_BAGL_ANY_ITERATOR_H_
