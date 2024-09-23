// Copyright 2024 Mikael Persson

#ifndef BAGL_BAGL_TRANSFORM_ITERATOR_H_
#define BAGL_BAGL_TRANSFORM_ITERATOR_H_

#include <iterator>

namespace bagl {

// A very simple, bare-bone transform iterator implementation
// since this is mostly what bagl needs, and it avoid bringing in
// a dependency for an iterators library.
// Derived class implementer just has to provide a public `dereference` function.
template <typename Derived, typename Value, typename BaseIter>
class transform_iterator {
 public:
  using self = transform_iterator<Derived, Value, BaseIter>;
  using iterator_type = BaseIter;
  using iterator_category = typename std::iterator_traits<BaseIter>::iterator_category;
  using value_type = Value;
  using difference_type = std::iter_difference_t<BaseIter>;
  using pointer = const Value*;
  using reference = const Value&;

  transform_iterator() = default;
  explicit transform_iterator(BaseIter curr) : curr_(curr) {}

  BaseIter base() const { return curr_; }
  value_type operator*() const { return derived().dereference(); }

  Derived& operator++() {
    ++curr_;
    return derived();
  }
  Derived operator++(int) {
    Derived result = derived();
    ++curr_;
    return result;
  }
  Derived& operator--() {
    --curr_;
    return derived();
  }
  Derived operator--(int) {
    Derived result = derived();
    --curr_;
    return result;
  }

  Derived operator+(difference_type n) const {
    Derived result = derived();
    result += n;
    return result;
  }
  Derived& operator+=(difference_type n) {
    curr_ += n;
    return derived();
  }
  Derived operator-(difference_type n) const {
    Derived result = derived();
    result -= n;
    return result;
  }
  Derived& operator-=(difference_type n) {
    curr_ -= n;
    return derived();
  }
  value_type operator[](difference_type n) const {
    Derived result = derived();
    result += n;
    return result.dereference();
  }

 protected:
  BaseIter curr_;

  Derived& derived() { return static_cast<Derived&>(*this); }
  const Derived& derived() const { return static_cast<const Derived&>(*this); }
};

}  // namespace bagl

#endif  // BAGL_BAGL_TRANSFORM_ITERATOR_H_
