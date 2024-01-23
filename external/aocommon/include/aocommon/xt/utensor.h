#ifndef AOCOMMON_XT_UTENSOR_H_
#define AOCOMMON_XT_UTENSOR_H_

#include <cstddef>

#include <xtensor/xtensor.hpp>

#include "../uvector.h"

namespace aocommon::xt {

// Define an xtensor type that uses aocommon::UVector for its storage.
// Like UVector, a UTensor does not initialize its contents.
// A UTensor is useful for value types without trivial default constructor,
// such as std::complex, since a regular xtensor only skips initialization for
// types that have a trivial default constructor.
template <class T, std::size_t N, ::xt::layout_type L = XTENSOR_DEFAULT_LAYOUT,
          class A = XTENSOR_DEFAULT_ALLOCATOR(T)>
using UTensor = ::xt::xtensor_container<aocommon::UVector<T, A>, N, L>;

}  // namespace aocommon::xt

#endif