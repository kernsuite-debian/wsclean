#ifndef AOCOMMON_VECTOR_MAP_H_
#define AOCOMMON_VECTOR_MAP_H_

#include <cstring>
#include <map>
#include <set>
#include <span>
#include <vector>

#include "uvector.h"

namespace aocommon {

/**
 * A map-like class whose key is a size_t, and the allowed keys are assumed to
 * span a small range around zero. It uses a vector indexed by the key, which
 * means that inserting a value with key 9 into the map will cause the map to
 * have 10 elements; 9 of which are default constructed (i.e., keys 0-8).
 *
 * It mixes functionality of a vector and a map, and aims at being fast.
 *
 * If it is necessary to know if a value exists in the map, the type can be
 * wrapped in a std::optional.
 */
template <typename T>
class VectorMap {
 public:
  using key_type = size_t;
  using mapped_type = T;
  using value_type = T;
  // bools won't work with std::vector, but non-trivial types don't work with
  // UVector, so use std::vector except for bools.
  using container = std::conditional_t<std::is_same_v<T, bool>, UVector<bool>,
                                       std::vector<T>>;
  using iterator = typename container::iterator;
  using const_iterator = typename container::const_iterator;

  VectorMap() = default;

  VectorMap(size_t size) : data_(size, T()) {}

  VectorMap(std::initializer_list<T> list) : data_(std::move(list)) {}

  /**
   * Construct from a std::map<size_t, T>. Afterwards, map[key] returns
   * the same as VectorMap[key].
   */
  VectorMap(const std::map<key_type, mapped_type>& map) {
    for (const std::pair<const key_type, mapped_type>& value : map) {
      AlwaysEmplace(value.first, value.second);
    }
  }

  /**
   * Get value belonging to a key. Unlike indexing into a std::map, the key must
   * exist in the vectormap.
   */
  T& operator[](key_type index) { return data_[index]; }

  const T& operator[](key_type index) const { return data_[index]; }

  /**
   * Construct a value in-place. The behaviour of this function is different
   * from std::map::emplace() in that it unconditionally assigns the value, even
   * when it already exists in the container. The reason for that is that a
   * VectorMap does not know if a value is explicitly inserted or created
   * because it was used as filling for inserting a higher key, i.e.:
   *
   *     vector_map.AlwaysEmplace(10, a);
   *     vector_map.AlwaysEmplace(5, b);
   *
   * If AlwaysEmplace() wouldn't store if the key exists, key 5 would be left
   * default constructed.
   */
  template <typename... Args>
  T& AlwaysEmplace(key_type index, Args&&... args) {
    if (index >= data_.size()) {
      Resize(index);
      return data_.emplace_back(std::forward<Args>(args)...);
    } else {
      data_[index] = T(std::forward<Args>(args)...);
      return data_[index];
    }
  }

  template <typename... Args>
  T& EmplaceBack(Args&&... args) {
    return data_.emplace_back(std::forward<Args>(args)...);
  }

  iterator begin() { return data_.begin(); }
  iterator end() { return data_.end(); }
  const_iterator begin() const { return data_.begin(); }
  const_iterator end() const { return data_.end(); }

  void Clear() { data_.clear(); }
  bool Empty() const { return data_.empty(); }
  /**
   * Like std::map::find(), returns an iterator to the value with the
   * given key, or returns end() otherwise.
   */
  iterator Find(key_type index) {
    return index < data_.size() ? data_.begin() + index : data_.end();
  }
  const_iterator Find(key_type index) const {
    return const_cast<VectorMap&>(*this).Find(index);
  }
  size_t Size() const { return data_.size(); }
  void Resize(size_t new_size) { data_.resize(new_size, T()); }

  /**
   * Returns a std::set that contains the keys for which the corresponding value
   * satisfies the given @p condition.
   * @param condition should be of the form bool(const value_type&).
   */
  template <typename Condition>
  std::set<size_t> GetKeysIf(Condition condition) const {
    std::set<size_t> result;
    for (size_t i = 0; i != data_.size(); ++i) {
      if (condition(data_[i])) result.emplace(i);
    }
    return result;
  }

  /**
   * Returns the count of values that satisfy the given @p condition.
   * @param condition should be of the form bool(const value_type&).
   */
  template <typename Condition>
  size_t CountKeysIf(Condition condition) const {
    size_t result = 0;
    for (size_t i = 0; i != data_.size(); ++i) {
      if (condition(data_[i])) ++result;
    }
    return result;
  }

 private:
  container data_;
};

}  // namespace aocommon

#endif
