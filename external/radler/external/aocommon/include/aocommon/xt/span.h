#ifndef AOCOMMON_XT_SPAN_H_
#define AOCOMMON_XT_SPAN_H_

#include <array>
#include <type_traits>
#include <utility>
#include <vector>

#include <xtensor/xadapt.hpp>

namespace aocommon {
namespace xt {
/**
 * Span is a generic non-owning XTensor view type.
 * Like std::span, this view allows modification of the underlying data
 * if the data type is non-const.
 * @tparam T Data type for the elements in the view. Using a const type makes
 * the view read-only.
 * @tparam Dimensions Number of dimensions in the view.
 */
template <typename T, size_t Dimensions>
using Span = decltype(::xt::adapt(std::add_pointer_t<T>{},
                                  std::array<size_t, Dimensions>{}));

/**
 * Create a Span from a raw pointer.
 * @tparam T Data type for the elements in the view.
 * @tparam Dimensions Number of dimensions in the view.
 * @param pointer Pointer to a contiguous list of elements.
 * @param shape The shape for the view.
 * @return A view of the elements behind the pointer.
 */
template <typename T, size_t Dimensions>
[[nodiscard]] Span<T, Dimensions> CreateSpan(
    T* pointer, const std::array<size_t, Dimensions>& shape) {
  // std::forward ensures that xt::adapt receives a pointer and not a
  // reference to a pointer.
  return ::xt::adapt(std::forward<T*>(pointer), shape);
}

/**
 * Creates a Span from a std::vector.
 * @tparam T Data type for the elements in the view.
 * @tparam Dimensions Number of dimensions in the view.
 * @param vector The underlying vector with elements.
 * @param shape The shape for the view. It must match the vector's size.
 * @return A view of the vector with modifiable elements.
 */
template <typename T, size_t Dimensions>
[[nodiscard]] Span<T, Dimensions> CreateSpan(
    std::vector<T>& vector, const std::array<size_t, Dimensions>& shape) {
  return ::xt::adapt(vector.data(), vector.size(), ::xt::no_ownership(), shape);
}

/**
 * Creates a Span from a const std::vector.
 * @tparam T Data type for the elements in the view.
 * @tparam Dimensions Number of dimensions in the view.
 * @param vector The underlying vector with elements.
 * @param shape The shape for the view. It must match the vector's size.
 * @return A view of the vector with non-modifiable elements.
 */
template <typename T, size_t Dimensions>
[[nodiscard]] Span<const T, Dimensions> CreateSpan(
    const std::vector<T>& vector, const std::array<size_t, Dimensions>& shape) {
  return ::xt::adapt(vector.data(), vector.size(), ::xt::no_ownership(), shape);
}

/**
 * Deleted overload, for avoiding creating spans on temporary objects.
 */
template <typename T, size_t Dimensions>
void CreateSpan(std::vector<T>&& vector,
                const std::array<size_t, Dimensions>& shape) = delete;

/**
 * Creates a Span from an xtensor or UTensor object.
 * @tparam TensorType xtensor or UTensor type.
 * @param tensor The underlying xtensor object.
 * @return A view of the tensor with modifiable elements.
 */
template <typename TensorType>
[[nodiscard]] Span<typename TensorType::value_type, TensorType::rank>
CreateSpan(TensorType& tensor) {
  return ::xt::adapt(tensor.data(), tensor.shape());
}

/**
 * Creates a Span from a const xtensor or const UTensor object.
 * @tparam TensorType xtensor or UTensor type.
 * @param tensor The underlying xtensor object.
 * @return A view of the tensor with non-modifiable elements.
 */
template <typename TensorType>
[[nodiscard]] Span<const typename TensorType::value_type, TensorType::rank>
CreateSpan(const TensorType& tensor) {
  return ::xt::adapt(tensor.data(), tensor.shape());
}

/**
 * Deleted overload, for avoiding creating spans on temporary objects.
 */
template <typename TensorType>
void CreateSpan(TensorType&& tensor) = delete;

}  // namespace xt
}  // namespace aocommon

#endif
