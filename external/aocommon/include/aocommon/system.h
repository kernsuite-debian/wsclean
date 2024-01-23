#ifndef AOCOMMON_SYSTEM_H_
#define AOCOMMON_SYSTEM_H_

#include <unistd.h>
#include <sched.h>

#include <cstring>
#include <stdexcept>

namespace aocommon {
namespace system {

inline long TotalMemory() {
  // Solution from
  // https://stackoverflow.com/questions/2513505/how-to-get-available-memory-c-g
  return sysconf(_SC_PHYS_PAGES) * sysconf(_SC_PAGE_SIZE);
}

inline size_t ProcessorCount() {
#ifdef __APPLE__
  return sysconf(_SC_NPROCESSORS_ONLN);
#else
  cpu_set_t cs;
  CPU_ZERO(&cs);
  sched_getaffinity(0, sizeof(cs), &cs);
  return CPU_COUNT(&cs);
#endif
}

namespace detail {
inline char* HandleStrReturn(int value) {
  if (value != 0) throw std::runtime_error("strerror_r() reported an error");
  return nullptr;
}
inline char* HandleStrReturn(char* value) { return value; }
}  // namespace detail

/**
 * @brief Convert a Posix \c errno error number to a string. This function is
 * thread safe.
 */
inline std::string GetErrorString(int errnum) {
  // Because strerror_r() has different return values on different platforms,
  // two overloads of handle_strerror are used to make this compile and work
  // in either case of int or char*.
  char buffer[1024];
  char* ret = detail::HandleStrReturn(strerror_r(errnum, buffer, 1024));
  if (ret == nullptr) {
    return std::string(buffer);
  } else {
    return std::string(ret);
  }
}

}  // namespace system
}  // namespace aocommon
#endif
