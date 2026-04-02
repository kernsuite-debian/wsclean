// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef SCHAAPCOMMON_H5PARM_H5PARM_H_
#define SCHAAPCOMMON_H5PARM_H5PARM_H_

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include "soltab.h"

namespace schaapcommon {
namespace h5parm {
class H5Parm : private H5::H5File {
 public:
  struct source_t {
    char name[128];
    /// in radians
    double dir[2];
  };

  /// Open existing H5Parm or create a new one
  /// Default name is given by sol_set_name, if that does not exist continue
  /// searching for sol000, sol001, etc.
  /// Only one solset of an H5Parm can be opened at once; this object only
  /// provides info about one SolSet (even though the file can contain more).
  explicit H5Parm(const std::string& filename, bool force_new = false,
                  bool force_new_sol_set = false,
                  const std::string& sol_set_name = "");

  /// Open existing H5 file and put the tables named in the second
  /// argument in the internal cache for fast access.
  explicit H5Parm(const std::string& filename,
                  const std::vector<std::string>& tables_in_cache);

  H5Parm();

  // Since sol_tabs_ contains std::unique_ptr, H5Parm is not copyable.
  H5Parm(const H5Parm&) = delete;
  H5Parm(H5Parm&&) = default;
  H5Parm& operator=(const H5Parm&) = delete;
  H5Parm& operator=(H5Parm&&) = default;

  ~H5Parm() override;

  /// Add metadata (directions on the sky in J2000) about named sources
  void AddSources(const std::vector<std::string>& names,
                  const std::vector<std::pair<double, double>>& dirs);

  /// Add metadata (positions on earth in ITRF) about antennas
  void AddAntennas(const std::vector<std::string>& names,
                   const std::vector<std::array<double, 3>>& positions);

  /// Add metadata about polarizations
  void AddPolarizations(const std::vector<std::string>& polarizations);

  /// Add a version stamp in the attributes of the group
  static void AddVersionStamp(H5::Group& node) {
    SolTab::AddVersionStamp(node);
  };

  /// Create and return a new soltab. Type is the type as used in BBS
  SolTab& CreateSolTab(const std::string& name, const std::string& type,
                       const std::vector<AxisInfo>& axes);

  SolTab& GetSolTab(const std::string& name);

  /// Get the name of the one SolSet used in this H5Parm
  std::string GetSolSetName() const;

  /// Get the number of SolTabs in the active solset of this h5parm
  size_t NumSolTabs() const { return sol_tabs_.size(); }

  /// Is the given soltab resent in the active solset of this h5parm
  bool HasSolTab(const std::string& sol_tab_name) const;

  /// Read the number of sources from the file.
  size_t GetNumSources() const;

  /// Read sources from the file.
  std::vector<source_t> GetSources() const;

  /// Get name of the source with coordinates closest to provided (ra, dec)
  /// coordinates (in J2000 and in radians)
  std::string GetNearestSource(double ra, double dec) const;

  /// Is the H5Parm implementation thread safe? If not, H5Parm does not support
  /// concurrent calls from multiple threads.
  static bool IsThreadSafe();

 private:
  struct antenna_t {
    char name[16];
    float position[3];
  };

  struct polarization_t {
    char name[2];
  };

  std::vector<source_t> ReadSourceTable(const H5::Group& sol_set) const;

  std::map<std::string, std::unique_ptr<SolTab>> sol_tabs_;

  H5::Group sol_set_;
};
}  // namespace h5parm
}  // namespace schaapcommon

#endif
