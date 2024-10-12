// Copyright 2007-2009 Andrew Sutton
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_CONSTANT_PROPERTY_MAP_H_
#define BAGL_BAGL_CONSTANT_PROPERTY_MAP_H_

#include <utility>

#include "bagl/property_map.h"

namespace bagl
{

/**
 * A constant property is one, that regardless of the edge or vertex given,
 * will always return a constant value.
 */
template < typename Key, typename Value >
class constant_property_map : public put_get_helper<constant_property_map< Key, Value > >
{
public:
    using key_type = Key;
    using value_type = Value;
    using reference = const Value&;

    constant_property_map() = default;
    explicit constant_property_map(value_type value) : value_(std::move(value)) {}

    template <typename T>
    reference operator[](T&& /*unused*/) const { return value_; }
private:
    value_type value_;
};

template < typename Key, typename U >
auto make_constant_property(U&& value)
{
    return constant_property_map< Key, std::decay_t<U> >(std::forward<U>(value));
}

/**
 * Same as above, but pretends to be writable as well.
 */
template < typename Key, typename Value > class constant_writable_property_map
{
public:
    using key_type = Key;
    using value_type = Value;
    using reference = Value&;

    constant_writable_property_map() = default;
    explicit constant_writable_property_map(value_type value) : value_(std::move(value)) {}

    friend const Value& get(const constant_writable_property_map& me, const Key& /*unused*/)
    {
        return me.value_;
    }
    friend void put(const constant_writable_property_map& /*unused*/, const Key& /*unused*/, const Value& /*unused*/) {}

private:
    value_type value_;
};

template < typename Key, typename U >
auto make_constant_writable_property(U&& value)
{
    return constant_writable_property_map< Key, std::decay_t<U> >(std::forward<U>(value));
}

} // namespace bagl

#endif // BAGL_BAGL_CONSTANT_PROPERTY_MAP_H_
