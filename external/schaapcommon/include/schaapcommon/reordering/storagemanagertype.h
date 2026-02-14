#ifndef WSCLEAN_MAIN_STORAGE_MANAGER_TYPE_H_
#define WSCLEAN_MAIN_STORAGE_MANAGER_TYPE_H_

#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>

namespace schaapcommon::reordering {

enum class StorageManagerType { Default, StokesI };

inline StorageManagerType GetStorageManagerType(
    const std::string& type_string) {
  const std::string lowercase_type = boost::to_lower_copy(type_string);
  if (lowercase_type == "default")
    return StorageManagerType::Default;
  else if (lowercase_type == "stokes_i" || lowercase_type == "stokes-i")
    return StorageManagerType::StokesI;
  else
    throw std::runtime_error("Unknown storage manager type: " + type_string);
}

}  // namespace schaapcommon::reordering

#endif
