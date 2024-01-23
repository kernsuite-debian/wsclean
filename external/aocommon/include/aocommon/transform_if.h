#ifndef AOCOMMON_TRANSFORM_IF_H_
#define AOCOMMON_TRANSFORM_IF_H_

namespace aocommon {

/**
 * A combination of std::copy_if and std::transform.
 *
 * There's no standard algorithm that can do a filtering transformation. This
 * function implements that.
 *
 * @param first The first element to transform.
 * @param last  Points one beyond the last element to transform.
 * @param out   The beginning of the output range.
 * @param pred  The unary predicate which returns @c true for the required
 *              arguments.
 * @param op    The unary operation to be applied to the required elements.
 *
 * @returns     One beyond the last element written to.
 *
 * Complexity:
 * - Exactly std::distance(first, last) applications of \a pred
 * - At most std::distance(first, last) applications of \a op
 */
template <class InputIt, class OutputIt, class UnaryPredicate,
          class UnaryOperation>
#if __cplusplus > 201402L
constexpr
#endif
    OutputIt
    transform_if(InputIt first, InputIt last, OutputIt out, UnaryPredicate pred,
                 UnaryOperation op) {
  while (first != last) {
    if (pred(*first)) {
      *out = op(*first);
      ++out;
    }
    ++first;
  }

  return out;
}

}  // namespace aocommon

#endif  // AO_TRANSFORM_IF_H
