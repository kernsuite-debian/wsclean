#ifndef AOCOMMON_PY_UNIQUEPTR_H_
#define AOCOMMON_PY_UNIQUEPTR_H_

#include <memory>
#include <stdexcept>
#include <utility>

namespace aocommon::py {

/// Wrapper class around a unique pointer, which allows binding C++ functions
/// that have a std::unique_ptr argument, e.g., Step::process() in DP3.
/// The binding should use PyUniquePointer<T> as argument instead of
/// std::unique_ptr<T>, and pass the unique pointer to C++ using
/// PyUniquePointer::take().
/// Python bindings for T should bind PyUniquePointer<T> instead of T.
/// If Python code uses the object after passing it to C++, PyUniquePointer
/// throws a runtime error.
/// @tparam T The base type that should be wrapped, e.g., DPBuffer in DP3.
template <class T>
class PyUniquePointer {
 public:
  /// Constructor, which allows construction from an existing unique pointer.
  explicit PyUniquePointer(std::unique_ptr<T> t) : pointer_(std::move(t)) {}

  /// Constructor. Forwards all arguments to T's constructor.
  template <typename... Args>
  PyUniquePointer(Args&&... args)
      : pointer_(std::make_unique<T>(std::forward<Args>(args)...)) {}

  PyUniquePointer(const PyUniquePointer&) = delete;
  PyUniquePointer& operator=(const PyUniquePointer&) = delete;
  PyUniquePointer(PyUniquePointer&&) = default;
  PyUniquePointer& operator=(PyUniquePointer&&) = default;

  /// @return A reference to the wrapped object.
  /// @throw std::runtime_error If the PyUniquePointer holds a null pointer.
  [[nodiscard]] T& operator*() {
    CheckPointer();
    return *pointer_;
  }

  /// @return A const reference to the wrapped object.
  /// @throw std::runtime_error If the PyUniquePointer holds a null pointer.
  [[nodiscard]] const T& operator*() const {
    CheckPointer();
    return *pointer_;
  }

  /// @return A pointer to the wrapped object.
  /// @throw std::runtime_error If the PyUniquePointer holds a null pointer.
  [[nodiscard]] T* operator->() {
    CheckPointer();
    return pointer_.get();
  }

  /// @return A const pointer to the wrapped object.
  /// @throw std::runtime_error If the PyUniquePointer holds a null pointer.
  [[nodiscard]] const T* operator->() const {
    CheckPointer();
    return pointer_.get();
  }

  /// Extract the wrapped object and invalidate this wrapper by turning
  /// the wrapped pointer into a null pointer.
  /// @throw std::runtime_error If the PyUniquePointer holds a null pointer.
  [[nodiscard]] std::unique_ptr<T> take() {
    CheckPointer();
    return std::unique_ptr<T>(pointer_.release());
  }

 private:
  void CheckPointer() const {
    if (!pointer_) {
      throw std::runtime_error(
          "Object is no longer valid. Ownership was transferred from Python to "
          "C++.");
    }
  }

  /// Pointer to the wrapped object. If it is null, the wrapper is invalid and
  /// using it will throw an exception.
  std::unique_ptr<T> pointer_;
};

}  // namespace aocommon::py

#endif