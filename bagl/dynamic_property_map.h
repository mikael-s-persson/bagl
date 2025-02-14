// Copyright 2004-5 The Trustees of Indiana University.
//  Authors: Doug Gregor
//           Ronald Garcia
// Copyright 2024 Mikael Persson - Modernized to C++20

#ifndef BAGL_BAGL_DYNAMIC_PROPERTY_MAP_H_
#define BAGL_BAGL_DYNAMIC_PROPERTY_MAP_H_

//  dynamic_property_map.hpp -
//    Support for runtime-polymorphic property maps.  This header is factored
//  out of Doug Gregor's routines for reading GraphML files for use in reading
//  GraphViz graph files.

#include <any>
#include <charconv>
#include <exception>
#include <functional>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

#include "bagl/graph_traits.h"
#include "bagl/property_map.h"

namespace bagl {

namespace dynamic_pmap_detail {

// read_value -
//   A wrapper around lexical_cast, which does not behave as
//   desired for std::string types.
template <typename Value>
Value read_value(const std::string& value) {
  Value result{};
  if (std::from_chars(value.data(), value.data() + value.size(), result).ec == std::errc{}) {
    return result;
  }
  return Value{};
}

template <>
inline float read_value<float>(const std::string& value) {
  return std::stof(value);
}
template <>
inline double read_value<double>(const std::string& value) {
  return std::stod(value);
}

template <>
inline bool read_value<bool>(const std::string& value) {
  // Any other ways to say true?
  return value == "true" || value == "True" || value == "TRUE" || value == "1";
}

template <>
inline std::string read_value<std::string>(const std::string& value) {
  return value;
}
template <>
inline std::string_view read_value<std::string_view>(const std::string& value) {
  return value;
}

}  // namespace dynamic_pmap_detail

// dynamic_property_map -
//  This interface supports polymorphic manipulation of property maps.
class dynamic_property_map {  // NOLINT
 public:
  virtual ~dynamic_property_map() = default;

  [[nodiscard]] virtual std::any get_value(const std::any& key) = 0;
  [[nodiscard]] virtual void* get_pointer(const std::any& key) = 0;
  [[nodiscard]] virtual const void* get_const_pointer(const std::any& key) = 0;
  [[nodiscard]] virtual std::string get_string(const std::any& key) = 0;
  virtual void put_value(const std::any& key, const std::any& value) = 0;
  virtual void put_value(const std::any& key, const std::string& value) = 0;
  [[nodiscard]] virtual const std::type_info& key() const = 0;
  [[nodiscard]] virtual bool is_key_of_type(const std::type_info& info) const = 0;
  [[nodiscard]] virtual const std::type_info& value() const = 0;
  [[nodiscard]] virtual bool is_value_of_type(const std::type_info& info) const = 0;

  template <typename Value>
  Value get_as(const std::any& key) {
    if constexpr (std::is_reference_v<Value>) {
      if constexpr (std::is_const_v<std::remove_reference_t<Value>>) {
        return *(static_cast<std::remove_reference_t<Value>*>(get_const_pointer(key)));
      } else {
        return *(static_cast<std::remove_reference_t<Value>*>(get_pointer(key)));
      }
    } else {
      return std::any_cast<Value>(get_value(key));
    }
  }
};

template <typename ValueType>
struct dynamic_property_map_by_ptr : public bagl::put_get_helper<dynamic_property_map_by_ptr<ValueType>> {
  using self = dynamic_property_map_by_ptr<ValueType>;
  using value_type = ValueType;

  std::shared_ptr<dynamic_property_map> p_parent;

  explicit dynamic_property_map_by_ptr(std::shared_ptr<dynamic_property_map> p_ptr) : p_parent(std::move(p_ptr)) {}
  value_type operator[](const std::any& k) const {
    if constexpr (std::is_const_v<value_type>) {
      return *(static_cast<std::remove_reference_t<value_type>*>(p_parent->get_const_pointer(k)));
    } else {
      return *(static_cast<std::remove_reference_t<value_type>*>(p_parent->get_pointer(k)));
    }
  }
};

template <typename ValueType>
struct dynamic_property_map_by_any : public bagl::put_get_helper<dynamic_property_map_by_any<ValueType>> {
  using self = dynamic_property_map_by_any<ValueType>;
  using value_type = ValueType;

  std::shared_ptr<dynamic_property_map> p_parent;

  explicit dynamic_property_map_by_any(std::shared_ptr<dynamic_property_map> p_ptr) : p_parent(std::move(p_ptr)) {}
  value_type operator[](const std::any& k) const { return std::any_cast<value_type>(p_parent->get_value(k)); }
};

// Property map exceptions

struct dynamic_property_exception : public std::exception {};

struct property_not_found : public dynamic_property_exception {
  std::string property;
  mutable std::string statement;
  explicit property_not_found(std::string a_property) : property(std::move(a_property)) {}

  [[nodiscard]] const char* what() const noexcept override {
    if (statement.empty()) {
      statement = std::string("Property not found: ") + property + ".";
    }

    return statement.c_str();
  }
};

struct dynamic_get_failure : public dynamic_property_exception {
  std::string property;
  mutable std::string statement;
  explicit dynamic_get_failure(std::string a_property) : property(std::move(a_property)) {}

  [[nodiscard]] const char* what() const noexcept override {
    if (statement.empty()) {
      statement = std::string("dynamic property get cannot retrieve value for property: ") + property + ".";
    }

    return statement.c_str();
  }
};

struct dynamic_const_put_error : public dynamic_property_exception {
  [[nodiscard]] const char* what() const noexcept override {
    return "Attempt to put a value into a const property map: ";
  }
};

namespace dynamic_pmap_detail {

// dynamic_property_map_adaptor -
//   property-map adaptor to support runtime polymorphism.
template <typename Key, typename PropertyMap>
class dynamic_property_map_adaptor : public dynamic_property_map {
  using key_type = Key;
  using value_type = property_traits_value_t<PropertyMap>;

  template <bool ForGetter = false>
  static decltype(auto) cast_key(const std::any& in_key) {
    if constexpr (std::is_pointer_v<key_type>) {
      // Handle both const T* and T*.
      using key_value = std::remove_pointer_t<key_type>;
      using non_const_key = std::remove_cv_t<key_value>;
      using ptr_to_const_key = std::add_pointer_t<std::add_const_t<non_const_key>>;
      using ptr_to_non_const_key = std::add_pointer_t<non_const_key>;
      if constexpr (std::is_const_v<key_value> || ForGetter) {
        if (in_key.type() == typeid(ptr_to_const_key)) {
          return *std::any_cast<ptr_to_const_key>(in_key);
        }
        if (in_key.type() == typeid(ptr_to_non_const_key)) {
          return std::as_const(*std::any_cast<ptr_to_non_const_key>(in_key));
        }
        // Probably will throw a bad_any_cast exception.
        return *std::any_cast<ptr_to_const_key>(in_key);
      } else {
        // Might throw a bad_any_cast exception.
        return *std::any_cast<ptr_to_non_const_key>(in_key);
      }
    } else {
      return std::any_cast<key_type>(in_key);
    }
  }

  // do_put - overloaded dispatches from the put() member function.
  //   Attempts to "put" to a property map that does not model
  //   WritablePropertyMap result in a runtime exception.

  //   in_value must either hold an object of value_type or a string that
  //   can be converted to value_type via iostreams.
  void do_put(const std::any& in_key, const std::any& in_value) {
    if (in_value.type() == std::type_index{typeid(value_type)}) {
      put(property_map_, cast_key(in_key), std::any_cast<value_type>(in_value));
    } else {
      //  if in_value is an empty string, put a default constructed value_type.
      auto v = std::any_cast<std::string>(in_value);
      if (v.empty()) {
        put(property_map_, cast_key(in_key), value_type());
      } else {
        put(property_map_, cast_key(in_key), dynamic_pmap_detail::read_value<value_type>(v));
      }
    }
  }

 public:
  explicit dynamic_property_map_adaptor(PropertyMap property_map) : property_map_(std::move(property_map)) {}

  [[nodiscard]] decltype(auto) get_value_impl(const std::any& in_key) {
    if constexpr (std::is_pointer_v<key_type>) {
      return property_map_[cast_key<true>(in_key)];
    } else {
      return get(property_map_, cast_key<true>(in_key));
    }
  }

  [[nodiscard]] std::any get_value(const std::any& in_key) override { return get_value_impl(in_key); }

  [[nodiscard]] void* get_pointer(const std::any& in_key) override {
    using ref_type = decltype(get_value_impl(in_key));
    if constexpr (std::is_reference_v<ref_type>) {
      if constexpr (!std::is_const_v<std::remove_reference_t<ref_type>>) {
        return static_cast<void*>(&get_value_impl(in_key));
      } else {
        throw dynamic_const_put_error();
      }
    } else {
      throw dynamic_const_put_error();
    }
  }

  [[nodiscard]] const void* get_const_pointer(const std::any& in_key) override {
    using ref_type = decltype(get_value_impl(in_key));
    if constexpr (std::is_reference_v<ref_type>) {
      return static_cast<const void*>(&get_value_impl(in_key));
    } else {
      throw dynamic_const_put_error();
    }
  }

  [[nodiscard]] std::string get_string(const std::any& in_key) override {
    std::ostringstream out;
    out << get_value_impl(in_key);
    return out.str();
  }

  void put_value(const std::any& in_key, const std::any& in_value) override {
    if constexpr (!concepts::WritablePropertyMap<PropertyMap, decltype(cast_key(in_key))>) {
      throw dynamic_const_put_error();
    } else {
      do_put(in_key, in_value);
    }
  }
  void put_value(const std::any& in_key, const std::string& in_value) override {
    if constexpr (!concepts::WritablePropertyMap<PropertyMap, decltype(cast_key(in_key))>) {
      throw dynamic_const_put_error();
    } else {
      put(property_map_, cast_key(in_key), dynamic_pmap_detail::read_value<value_type>(in_value));
    }
  }

  [[nodiscard]] const std::type_info& key() const override { return typeid(key_type); }
  [[nodiscard]] const std::type_info& value() const override { return typeid(value_type); }
  [[nodiscard]] bool is_key_of_type(const std::type_info& info) const override {
    if constexpr (std::is_pointer_v<key_type>) {
      // Handle both const T* and T*.
      using non_const_key = std::remove_cv_t<std::remove_pointer_t<key_type>>;
      using ptr_to_const_key = std::add_pointer_t<std::add_const_t<non_const_key>>;
      using ptr_to_non_const_key = std::add_pointer_t<non_const_key>;
      return (info == typeid(ptr_to_const_key)) || (info == typeid(ptr_to_non_const_key));
    } else {
      return (info == typeid(key_type));
    }
  }
  [[nodiscard]] bool is_value_of_type(const std::type_info& info) const override {
    if constexpr (std::is_pointer_v<value_type>) {
      // Handle both const T* and T*.
      using non_const_key = std::remove_cv_t<std::remove_pointer_t<value_type>>;
      using ptr_to_const_key = std::add_pointer_t<std::add_const_t<non_const_key>>;
      using ptr_to_non_const_key = std::add_pointer_t<non_const_key>;
      return (info == typeid(ptr_to_const_key)) || (info == typeid(ptr_to_non_const_key));
    } else {
      return (info == typeid(value_type));
    }
  }

  [[nodiscard]] PropertyMap& base() { return property_map_; }
  [[nodiscard]] const PropertyMap& base() const { return property_map_; }

 private:
  PropertyMap property_map_;
};

}  // namespace dynamic_pmap_detail

//
// dynamic_properties -
//   container for dynamic property maps
//
struct dynamic_properties {
  using property_maps_type = std::unordered_multimap<std::string, std::shared_ptr<dynamic_property_map>>;
  using generate_fn_type =
      std::function<std::shared_ptr<dynamic_property_map>(const std::string&, const std::any&, const std::any&)>;

 public:
  using iterator = property_maps_type::iterator;
  using const_iterator = property_maps_type::const_iterator;

  dynamic_properties() = default;
  explicit dynamic_properties(generate_fn_type g) : generate_fn_(std::move(g)) {}

  template <typename Key, typename PropertyMap>
  dynamic_properties& property(std::string_view name, PropertyMap property_map_) {
    property_maps_.emplace(
        std::string{name},
        std::make_shared<dynamic_pmap_detail::dynamic_property_map_adaptor<Key, PropertyMap>>(property_map_));

    return *this;
  }

  // Add property maps from tags:
  // property(name, vertex_index, g) == property<vertex_descriptor>(name, get(vertex_index, g))
  template <typename PropertyTag, typename Graph>
  dynamic_properties& property(std::string_view name, PropertyTag property_tag, Graph&& g) {
    if constexpr (is_vertex_property_kind_v<PropertyTag>) {
      return property<graph_vertex_descriptor_t<Graph>>(std::move(name), get(property_tag, std::forward<Graph>(g)));
    } else if constexpr (is_edge_property_kind_v<PropertyTag>) {
      return property<graph_edge_descriptor_t<Graph>>(std::move(name), get(property_tag, std::forward<Graph>(g)));
    } else {
      return property<std::remove_reference_t<Graph>*>(std::move(name), get(property_tag, std::forward<Graph>(g)));
    }
  }
  // Add property maps from tags:
  // property(vertex_index, g) == property<vertex_descriptor>("vertex_index", get(vertex_index, g))
  template <typename PropertyTag, typename Graph>
  std::enable_if_t<!std::is_convertible_v<PropertyTag, std::string_view>, dynamic_properties&> property(
      PropertyTag property_tag, Graph&& g) {
    return property(std::string{PropertyTag::name}, property_tag, std::forward<Graph>(g));
  }

  template <typename Key, typename PropertyMap>
  [[nodiscard]] dynamic_properties property(std::string_view name, PropertyMap property_map_) const {
    dynamic_properties result = *this;
    result.property<Key>(std::string{name}, property_map_);
    return result;
  }

  // Add property maps from tags:
  // property(name, vertex_index, g) == property<vertex_descriptor>(name, get(vertex_index, g))
  template <typename PropertyTag, typename Graph>
  dynamic_properties property(std::string_view name, PropertyTag property_tag, Graph&& g) const {
    dynamic_properties result = *this;
    result.property(property_tag, std::forward<Graph>(g));
    return result;
  }
  // Add property maps from tags:
  // property(vertex_index, g) == property<vertex_descriptor>("vertex_index", get(vertex_index, g))
  template <typename PropertyTag, typename Graph>
  std::enable_if_t<!std::is_convertible_v<PropertyTag, std::string_view>, dynamic_properties> property(
      PropertyTag property_tag, Graph&& g) const {
    dynamic_properties result = *this;
    result.property(std::string{PropertyTag::name}, property_tag, std::forward<Graph>(g));
    return result;
  }

  [[nodiscard]] iterator begin() { return property_maps_.begin(); }
  [[nodiscard]] const_iterator begin() const { return property_maps_.begin(); }
  [[nodiscard]] iterator end() { return property_maps_.end(); }
  [[nodiscard]] const_iterator end() const { return property_maps_.end(); }

  [[nodiscard]] auto equal_range(const std::string& name) { return property_maps_.equal_range(name); }
  [[nodiscard]] auto equal_range(const std::string& name) const { return property_maps_.equal_range(name); }

  void insert(std::string name, std::shared_ptr<dynamic_property_map> pm) {
    property_maps_.emplace(std::move(name), std::move(pm));
  }

  template <typename Key, typename Value>
  [[nodiscard]] std::shared_ptr<dynamic_property_map> generate(const std::string& name, const Key& key,
                                                               const Value& value) {
    if (!generate_fn_) {
      throw property_not_found(name);
    }
    return generate_fn_(name, key, value);
  }

 private:
  property_maps_type property_maps_;
  generate_fn_type generate_fn_;
};

template <typename Value>
auto get_dynamic_property_map(const std::string& name, const dynamic_properties& dp) {
  for (auto [i, i_end] = dp.equal_range(name); i != i_end; ++i) {
    if (i->second->is_value_of_type(typeid(std::decay_t<Value>))) {
      if constexpr (std::is_reference_v<Value>) {
        return dynamic_property_map_by_ptr<Value>(i->second);
      } else {
        return dynamic_property_map_by_any<Value>(i->second);
      }
    }
  }
  throw dynamic_get_failure(name);
}

template <typename Key, typename Value>
bool put(const std::string& name, dynamic_properties& dp, const Key& key, const Value& value) {
  for (auto [i, i_end] = dp.equal_range(name); i != i_end; ++i) {
    if constexpr (std::is_same_v<Key, std::any>) {
      if (i->second->is_key_of_type(key.type())) {
        i->second->put_value(key, value);
        return true;
      }
    } else {
      if (i->second->is_key_of_type(typeid(key))) {
        i->second->put_value(key, value);
        return true;
      }
    }
  }

  auto new_map = dp.generate(name, key, value);
  if (new_map != nullptr) {
    new_map->put_value(key, value);
    dp.insert(name, std::move(new_map));
    return true;
  }
  return false;
}

template <typename Value, typename Key>
Value get(const std::string& name, const dynamic_properties& dp, const Key& key) {
  for (auto [i, i_end] = dp.equal_range(name); i != i_end; ++i) {
    if constexpr (std::is_same_v<Key, std::any>) {
      if (i->second->is_key_of_type(key.type())) {
        return i->second->get_as<Value>(key);
      }
    } else {
      if (i->second->is_key_of_type(typeid(key))) {
        return i->second->get_as<Value>(key);
      }
    }
  }
  throw dynamic_get_failure(name);
}

template <typename Value, typename Key>
Value get(const std::string& name, const dynamic_properties& dp, const Key& key, Value* /*unused*/) {
  for (auto [i, i_end] = dp.equal_range(name); i != i_end; ++i) {
    if constexpr (std::is_same_v<Key, std::any>) {
      if (i->second->is_key_of_type(key.type())) {
        return i->second->get_as<Value>(key);
      }
    } else {
      if (i->second->is_key_of_type(typeid(key))) {
        return i->second->get_as<Value>(key);
      }
    }
  }
  throw dynamic_get_failure(name);
}

template <typename Key>
std::string get(const std::string& name, const dynamic_properties& dp, const Key& key) {
  for (auto [i, i_end] = dp.equal_range(name); i != i_end; ++i) {
    if constexpr (std::is_same_v<Key, std::any>) {
      if (i->second->is_key_of_type(key.type())) {
        return i->second->get_string(key);
      }
    } else {
      if (i->second->is_key_of_type(typeid(key))) {
        return i->second->get_string(key);
      }
    }
  }
  throw dynamic_get_failure(name);
}

// The easy way to ignore properties.
inline std::shared_ptr<dynamic_property_map> ignore_other_properties(const std::string& /*unused*/,
                                                                     const std::any& /*unused*/,
                                                                     const std::any& /*unused*/) {
  return {};
}

}  // namespace bagl

#endif  // BAGL_BAGL_DYNAMIC_PROPERTY_MAP_H_
