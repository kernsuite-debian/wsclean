// Copyright (C) 2020 ASTRON (Netherlands Institute for Radio Astronomy)
// SPDX-License-Identifier: GPL-3.0-or-later

#include <boost/test/unit_test.hpp>

#include <array>
#include <sstream>
#include <vector>

#include "h5parm.h"
#include "h5cache.h"
#include "soltab.h"

using schaapcommon::h5parm::AxisInfo;
using schaapcommon::h5parm::H5Cache;
using schaapcommon::h5parm::H5Parm;
using schaapcommon::h5parm::SolTab;
using schaapcommon::h5parm::TimesAndFrequencies;

BOOST_AUTO_TEST_SUITE(h5parm)

namespace {

const size_t kNumAntennas = 3;
const size_t kNumFrequencies = 4;
const size_t kNumTimes = 7;
const size_t kNumDirections = 4;
const size_t kNumPolarizations = 2;
const size_t kNumValues = kNumAntennas * kNumFrequencies * kNumTimes *
                          kNumDirections * kNumPolarizations;
const size_t kNumWeights = kNumValues;

void CheckAxes(SolTab& soltab, size_t ntimes) {
  BOOST_CHECK_EQUAL(soltab.NumAxes(), size_t{3});
  BOOST_CHECK(soltab.HasAxis("ant"));
  BOOST_CHECK(soltab.HasAxis("time"));
  BOOST_CHECK(soltab.HasAxis("bla"));
  BOOST_CHECK_EQUAL(soltab.GetAxis(0).name, "ant");
  BOOST_CHECK_EQUAL(soltab.GetAxis(1).name, "time");
  BOOST_CHECK_EQUAL(soltab.GetAxis(2).name, "bla");
  BOOST_CHECK_EQUAL(soltab.GetAxis(0).size, size_t{3});
  BOOST_CHECK_EQUAL(soltab.GetAxis(1).size, ntimes);
  BOOST_CHECK_EQUAL(soltab.GetAxis(2).size, size_t{1});
}

void InitializeH5(H5Parm& h5parm) {
  // Add some metadata
  std::vector<std::string> antNames;
  std::vector<std::array<double, 3>> antPositions;
  for (size_t i = 0; i < 5; ++i) {
    std::stringstream antNameStr;
    antNameStr << "Antenna" << i;
    antNames.push_back(antNameStr.str());
    antPositions.emplace_back();
  }
  h5parm.AddAntennas(antNames, antPositions);
  h5parm.AddSources({"aaa", "bbb", "ccc", "ddd"},
                    {std::make_pair(0.0, 0.0), std::make_pair(0.0, 1.0),
                     std::make_pair(1.0, 1.0), std::make_pair(1.0, 0.0)});

  std::vector<AxisInfo> axes;
  axes.push_back(AxisInfo{"ant", 3});
  axes.push_back(AxisInfo{"time", kNumTimes});
  axes.push_back(AxisInfo{"bla", 1});
  h5parm.CreateSolTab("mysol", "mytype", axes);

  std::vector<AxisInfo> axes_freq;
  axes_freq.push_back(AxisInfo{"ant", kNumAntennas});
  axes_freq.push_back(AxisInfo{"freq", kNumFrequencies});
  h5parm.CreateSolTab("mysolwithfreq", "mytype", axes_freq);

  std::vector<AxisInfo> axes_time_first;
  axes_time_first.push_back(AxisInfo{"ant", kNumAntennas});
  axes_time_first.push_back(AxisInfo{"time", kNumTimes});
  axes_time_first.push_back(AxisInfo{"freq", kNumFrequencies});
  h5parm.CreateSolTab("timefreq", "mytype", axes_time_first);

  std::vector<AxisInfo> axes_freq_first;
  axes_freq_first.push_back(AxisInfo{"ant", kNumAntennas});
  axes_freq_first.push_back(AxisInfo{"freq", kNumFrequencies});
  axes_freq_first.push_back(AxisInfo{"time", kNumTimes});
  h5parm.CreateSolTab("freqtime", "mytype", axes_freq_first);

  std::vector<AxisInfo> axes_complete_set;
  axes_complete_set.push_back(AxisInfo{"ant", kNumAntennas});
  axes_complete_set.push_back(AxisInfo{"freq", kNumFrequencies});
  axes_complete_set.push_back(AxisInfo{"time", kNumTimes});
  axes_complete_set.push_back(AxisInfo{"dir", kNumDirections});
  axes_complete_set.push_back(AxisInfo{"pol", kNumPolarizations});
  h5parm.CreateSolTab("completeset", "mytype", axes_complete_set);

  std::vector<AxisInfo> axes_canonical_set;
  axes_canonical_set.push_back(AxisInfo{"time", kNumTimes});
  axes_canonical_set.push_back(AxisInfo{"freq", kNumFrequencies});
  axes_canonical_set.push_back(AxisInfo{"ant", kNumAntennas});
  axes_canonical_set.push_back(AxisInfo{"dir", kNumDirections});
  axes_canonical_set.push_back(AxisInfo{"pol", kNumPolarizations});
  h5parm.CreateSolTab("canonical", "mytype", axes_canonical_set);

  std::vector<AxisInfo> axes_non_canonical_set;
  axes_non_canonical_set.push_back(AxisInfo{"pol", kNumPolarizations});
  axes_non_canonical_set.push_back(AxisInfo{"time", kNumTimes});
  axes_non_canonical_set.push_back(AxisInfo{"dir", kNumDirections});
  axes_non_canonical_set.push_back(AxisInfo{"ant", kNumAntennas});
  axes_non_canonical_set.push_back(AxisInfo{"freq", kNumFrequencies});
  h5parm.CreateSolTab("noncanonical", "mytype", axes_non_canonical_set);
}

void SetSolTabMeta(SolTab& soltab, bool set_freq_meta, bool set_time_meta) {
  // Add metadata for stations
  const std::vector<std::string> someAntNames = {"Antenna1", "Antenna12",
                                                 "Antenna123"};
  soltab.SetAntennas(someAntNames);

  if (set_freq_meta) {
    // Add metadata for freqs;
    const std::vector<double> freqs{130e6, 131e6, 135e6, 137e6};
    soltab.SetFreqs(freqs);
  }

  if (set_time_meta) {
    // Add metadata for times
    std::vector<double> times;
    for (size_t time = 0; time < kNumTimes; ++time) {
      times.push_back(57878.5 + 2.0 * time);
    }
    soltab.SetTimes(times);
  }
}

void SetCompleteDataMeta(SolTab& soltab) {
  // Add metadata for antennas
  const std::vector<std::string> some_antenna_names = {"Tile11", "Tile18",
                                                       "Tile27"};
  soltab.SetAntennas(some_antenna_names);

  // Add metadata for frequencies
  std::vector<double> frequencies(kNumFrequencies);
  std::iota(frequencies.begin(), frequencies.end(), 130.0e6);
  soltab.SetFreqs(frequencies);

  // Add metadata for times
  std::vector<double> times(kNumTimes);
  std::iota(times.begin(), times.end(), 57800.0);
  soltab.SetTimes(times);

  // Add metadata for directions
  const std::vector<std::string> some_direction_names = {
      "direction_00", "direction_01", "direction_02", "direction_03"};
  soltab.SetSources(some_direction_names);

  // Add metadata for polarizations
  const std::vector<std::string> some_polarization_names = {"X", "Y"};
  soltab.SetPolarizations(some_polarization_names);
}

void FillData(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("mysol");

  // Add some data
  std::vector<double> vals(kNumAntennas * kNumTimes);
  std::vector<double> weights(kNumAntennas * kNumTimes);
  for (size_t ant = 0; ant < kNumAntennas; ++ant) {
    for (size_t time = 0; time < kNumTimes; ++time) {
      vals[ant * kNumTimes + time] = 10 * ant + time;
      weights[ant * kNumTimes + time] = 0.4;
    }
  }

  soltab.SetValues(vals, weights, "CREATE with SchaapCommon-test");

  SetSolTabMeta(soltab, false, true);
}

void FillDataWithFreqAxis(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("mysolwithfreq");

  // Add some data
  const std::vector<double> vals(kNumAntennas * kNumFrequencies, 1.0);
  const std::vector<double> weights(kNumAntennas * kNumFrequencies, 1.0);

  soltab.SetValues(vals, weights, "CREATE with SchaapCommon-test");
  SetSolTabMeta(soltab, true, false);
}

void FillDataTimeFirst(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("timefreq");

  // Add some data
  std::vector<double> vals(kNumAntennas * kNumTimes * kNumFrequencies);
  std::vector<double> weights(kNumAntennas * kNumTimes * kNumFrequencies, 0);
  for (size_t ant = 0; ant < kNumAntennas; ++ant) {
    for (size_t time = 0; time < kNumTimes; ++time) {
      for (size_t freq = 0; freq < kNumFrequencies; ++freq) {
        vals[ant * kNumTimes * kNumFrequencies + time * kNumFrequencies +
             freq] = ant * time * freq;
      }
    }
  }

  // Add some data
  soltab.SetValues(vals, weights, "CREATE with SchaapCommon-test");
  SetSolTabMeta(soltab, true, true);
}

void FillDataFreqFirst(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("freqtime");

  // Add some data
  std::vector<double> vals(kNumAntennas * kNumTimes * kNumFrequencies);
  std::vector<double> weights(kNumAntennas * kNumTimes * kNumFrequencies, 0);
  for (size_t ant = 0; ant < kNumAntennas; ++ant) {
    for (size_t freq = 0; freq < kNumFrequencies; ++freq) {
      for (size_t time = 0; time < kNumTimes; ++time) {
        vals[ant * kNumTimes * kNumFrequencies + freq * kNumTimes + time] =
            ant * time * freq;
      }
    }
  }

  // Add some data
  soltab.SetValues(vals, weights, "CREATE with SchaapCommon-test");
  SetSolTabMeta(soltab, true, true);
}

void FillDataCompleteSet(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("completeset");

  // Create some data
  std::vector<double> values(kNumValues, 1);
  std::vector<double> weights(kNumWeights, 1);

  // Add some data
  soltab.SetValues(values, weights, "CREATE with SchaapCommon-test");
  SetCompleteDataMeta(soltab);
}

void FillDataCanonicalSet(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("canonical");

  // Create 5-dimensional data. Each value encodes the index information.
  // For instance: The value 52143 was generated for indices:
  // time         -> 5
  // frequency    -> 2
  // antenna      -> 1
  // direction    -> 4
  // polarization -> 3
  // This encoding allows to verify that we recover the correct values
  // independent of the order in which they are stored.
  // In this case it has the canonical order: [time, freq, ant, dir, pol].
  xt::xtensor<double, 5> original_values;
  const std::array<size_t, 5> shape = {kNumTimes, kNumFrequencies, kNumAntennas,
                                       kNumDirections, kNumPolarizations};
  original_values.resize(shape);

  for (size_t time = 0; time < kNumTimes; ++time) {
    for (size_t frequency = 0; frequency < kNumFrequencies; ++frequency) {
      for (size_t antenna = 0; antenna < kNumAntennas; ++antenna) {
        for (size_t direction = 0; direction < kNumDirections; ++direction) {
          for (size_t polarization = 0; polarization < kNumPolarizations;
               ++polarization) {
            original_values(time, frequency, antenna, direction, polarization) =
                time * 10000 + frequency * 1000 + antenna * 100 +
                direction * 10 + polarization;
          }
        }
      }
    }
  }

  std::vector<double> values(kNumValues, 1);
  std::copy(original_values.begin(), original_values.end(), values.begin());

  std::vector<double> weights(kNumWeights, 1);

  // Add some data
  soltab.SetValues(values, weights, "CREATE with SchaapCommon-test");
  SetCompleteDataMeta(soltab);
}

void FillDataNonCanonicalSet(H5Parm& h5parm) {
  SolTab soltab = h5parm.GetSolTab("noncanonical");

  // Create 5-dimensional data. Each value encodes the index information.
  // For instance: The value 52143 was generated for indices:
  // time         -> 5
  // frequency    -> 2
  // antenna      -> 1
  // direction    -> 4
  // polarization -> 3
  // This encoding allows to verify that we recover the correct values
  // independent of the order in which they are stored.
  // In this case it was done on purpose to be different from the
  // canonical order.
  xt::xtensor<double, 5> original_values;
  const std::array<size_t, 5> shape = {kNumPolarizations, kNumTimes,
                                       kNumDirections, kNumAntennas,
                                       kNumFrequencies};
  original_values.resize(shape);

  for (size_t polarization = 0; polarization < kNumPolarizations;
       ++polarization) {
    for (size_t time = 0; time < kNumTimes; ++time) {
      for (size_t direction = 0; direction < kNumDirections; ++direction) {
        for (size_t antenna = 0; antenna < kNumAntennas; ++antenna) {
          for (size_t frequency = 0; frequency < kNumFrequencies; ++frequency) {
            original_values(polarization, time, direction, antenna, frequency) =
                time * 10000 + frequency * 1000 + antenna * 100 +
                direction * 10 + polarization;
          }
        }
      }
    }
  }

  // Dump the 5-dimensional array into a 1-dimensiona array.
  // Soltab keep the information about axis order.
  std::vector<double> values(kNumValues, 1);
  std::copy(original_values.begin(), original_values.end(), values.begin());
  std::vector<double> weights(kNumWeights, 1);

  // Add the data into the soltab.
  soltab.SetValues(values, weights, "CREATE with SchaapCommon-test");
  SetCompleteDataMeta(soltab);
}

struct H5Fixture {
  H5Fixture() {
    H5Parm h5parm("tH5Parm_tmp.h5", true);
    InitializeH5(h5parm);
    FillData(h5parm);
    FillDataWithFreqAxis(h5parm);
    FillDataTimeFirst(h5parm);
    FillDataFreqFirst(h5parm);
    FillDataCompleteSet(h5parm);
    FillDataCanonicalSet(h5parm);
    FillDataNonCanonicalSet(h5parm);
  }

  ~H5Fixture() { remove("tH5Parm_tmp.h5"); }
};
}  // namespace

BOOST_AUTO_TEST_CASE(create) {
  // Create a new H5Parm
  H5Parm h5parm("tH5Parm_tmp.h5", true);

  // Check that something is created
  BOOST_CHECK_EQUAL(((H5::H5File&)(h5parm)).getNumObjs(), 1u);

  // Check the name of the new solset "sol000"
  BOOST_CHECK_EQUAL(h5parm.GetSolSetName(), "sol000");

  InitializeH5(h5parm);

  // Check that the soltab exists
  BOOST_CHECK_EQUAL(h5parm.NumSolTabs(), 7);
  BOOST_CHECK(h5parm.HasSolTab("mysol"));

  // Check the axes
  SolTab soltab = h5parm.GetSolTab("mysol");
  BOOST_CHECK_EQUAL(soltab.GetType(), "mytype");
  CheckAxes(soltab, kNumTimes);
}

BOOST_FIXTURE_TEST_CASE(new_soltab, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", true, true, "harry");
  BOOST_CHECK_EQUAL(h5parm.GetSolSetName(), "harry");
  BOOST_CHECK_EQUAL(h5parm.NumSolTabs(), 0u);
}

BOOST_FIXTURE_TEST_CASE(existing_soltab, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");
  BOOST_CHECK_EQUAL(h5parm.GetSolSetName(), "sol000");
  BOOST_CHECK_EQUAL(h5parm.NumSolTabs(), 7u);
  BOOST_CHECK(h5parm.HasSolTab("mysol"));
  BOOST_CHECK(!h5parm.HasSolTab("nonexistingsol"));
}

BOOST_FIXTURE_TEST_CASE(axes, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");

  // Check the axes
  SolTab soltab = h5parm.GetSolTab("mysol");
  BOOST_CHECK_EQUAL(soltab.GetType(), "mytype");
  CheckAxes(soltab, kNumTimes);

  BOOST_CHECK_EQUAL(h5parm.GetNumSources(), 4u);

  // Return and check nearest source
  BOOST_CHECK_EQUAL(h5parm.GetNearestSource(0.0, 0.1), "aaa");
  BOOST_CHECK_EQUAL(h5parm.GetNearestSource(1.0, 0.51), "ccc");
  BOOST_CHECK_EQUAL(h5parm.GetNearestSource(1.0, 0.49), "ddd");

  double starttime = 57878.49999;
  hsize_t starttimeindex = soltab.GetTimeIndex(starttime);
  std::vector<double> val = soltab.GetValues("Antenna12", starttimeindex,
                                             kNumTimes, 1, 0, 4, 0, 4, 0);
  BOOST_CHECK_CLOSE(val[0], 10., 1e-8);
  BOOST_CHECK_CLOSE(val[1], 11., 1e-8);
  BOOST_CHECK_CLOSE(val[2], 12., 1e-8);
  BOOST_CHECK_CLOSE(val[3], 13., 1e-8);

  starttime = 57880.5;
  starttimeindex = soltab.GetTimeIndex(starttime);
  BOOST_CHECK_EQUAL(starttimeindex, hsize_t{1});
  std::vector<double> val2 =
      soltab.GetValues("Antenna123", starttimeindex, 2, 2, 0, 4, 0, 4, 0);

  BOOST_CHECK_CLOSE(val2[0], 21., 1e-8);
  BOOST_CHECK_CLOSE(val2[1], 23., 1e-8);
  BOOST_CHECK_CLOSE(soltab.GetTimeInterval(), 2., 1e-8);

  const std::vector<std::string>& antennas = soltab.GetStringAxis("ant");
  BOOST_CHECK_EQUAL(antennas.size(), size_t{3});
  BOOST_CHECK_EQUAL(antennas[0], "Antenna1");
  BOOST_CHECK_EQUAL(antennas[1], "Antenna12");
  BOOST_CHECK_EQUAL(antennas[2], "Antenna123");
}

BOOST_FIXTURE_TEST_CASE(grid_interpolation, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");
  SolTab& soltab = h5parm.GetSolTab("mysol");

  const std::vector<double> freqs{130e6, 131e6};

  const std::vector<double> times1{57878.5, 57880.5, 57882.5, 57884.5,
                                   57886.5, 57888.5, 57890.5};

  std::vector<double> newgridvals =
      soltab.GetValues("Antenna1", times1, freqs, 0, 0, true);
  BOOST_REQUIRE_EQUAL(newgridvals.size(), times1.size() * freqs.size());
  size_t idx = 0;
  for (size_t time = 0; time < times1.size(); ++time) {
    for (size_t freq = 0; freq < freqs.size(); ++freq) {
      BOOST_CHECK_CLOSE(newgridvals[idx], time, 1e-8);
      ++idx;
    }
  }

  std::vector<double> times2;
  for (size_t time = 0; time < 3 * times1.size() + 2; ++time) {
    times2.push_back(57878.5 + 2.0 * time / 3.);
  }

  newgridvals = soltab.GetValues("Antenna1", times2, freqs, 0, 0, true);
  BOOST_REQUIRE_EQUAL(newgridvals.size(), times2.size() * freqs.size());
  idx = 0;
  for (size_t time = 0; time < times2.size(); ++time) {
    for (size_t freq = 0; freq < freqs.size(); ++freq) {
      BOOST_CHECK_CLOSE(newgridvals[idx],
                        std::min((time + 1) / 3, times1.size() - 1), 1e-8);
      ++idx;
    }
  }
}

BOOST_FIXTURE_TEST_CASE(interpolate_single_time, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");
  SolTab& soltab = h5parm.GetSolTab("mysol");

  const std::vector<double> freqs{130e6, 131e6};
  const std::vector<double> times{57000, 57878.7, 57880.3, 57890.3, 58000};
  const std::vector<double> expected_nearest{10, 10, 11, 16, 16};
  const std::vector<double> expected_bilinear{10.0, 10.1, 10.9, 15.9, 16.0};

  for (size_t time = 0; time < times.size(); ++time) {
    const std::vector<double> nearest_vals =
        soltab.GetValues("Antenna12", {times[time]}, freqs, 0, 0, true);
    const std::vector<double> bilinear_vals =
        soltab.GetValues("Antenna12", {times[time]}, freqs, 0, 0, false);

    BOOST_REQUIRE_EQUAL(nearest_vals.size(), 1 * freqs.size());
    BOOST_REQUIRE_EQUAL(bilinear_vals.size(), 1 * freqs.size());

    for (size_t freq = 0; freq < freqs.size(); ++freq) {
      BOOST_CHECK_CLOSE(nearest_vals[freq], expected_nearest[time], 1e-8);
      BOOST_CHECK_CLOSE(bilinear_vals[freq], expected_bilinear[time], 1e-8);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(freq_interval_and_index, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");
  const SolTab& soltab = h5parm.GetSolTab("mysolwithfreq");
  BOOST_CHECK_CLOSE(soltab.GetFreqInterval(0), 1.0e6, 1.0e-8);
  BOOST_CHECK_CLOSE(soltab.GetFreqInterval(1), 4.0e6, 1.0e-8);
  BOOST_CHECK_CLOSE(soltab.GetFreqInterval(2), 2.0e6, 1.0e-8);

  BOOST_CHECK_THROW(soltab.GetFreqIndex(128.0e6),
                    std::runtime_error);               // Too far from lowest
  BOOST_CHECK_EQUAL(soltab.GetFreqIndex(129.1e6), 0);  // closest to 130e6
  BOOST_CHECK_EQUAL(soltab.GetFreqIndex(130.4e6), 0);  // closest to 130e6
  BOOST_CHECK_EQUAL(soltab.GetFreqIndex(130.6e6), 1);  // closest to 131e6
  BOOST_CHECK_EQUAL(soltab.GetFreqIndex(136.1e6), 3);  // closest to 137e6
  BOOST_CHECK_EQUAL(soltab.GetFreqIndex(137.0e6), 3);  // closest to 137e6
  BOOST_CHECK_EQUAL(soltab.GetFreqIndex(137.8e6), 3);  // closest to 137e6
  BOOST_CHECK_THROW(soltab.GetFreqIndex(150.0e6),
                    std::runtime_error);  // Too far from highest
}

BOOST_FIXTURE_TEST_CASE(axis_ordering, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");

  SolTab soltab_tf = h5parm.GetSolTab("timefreq");
  SolTab soltab_ft = h5parm.GetSolTab("freqtime");

  const std::vector<double> freqs{130e6, 135e6, 131e6};
  const std::vector<double> times{57000, 57878.7, 57880.3, 57890.3, 58000};

  const std::vector<double> nearest_tf = soltab_tf.GetValuesOrWeights(
      "val", "Antenna12", times, freqs, 0, 0, true);
  const std::vector<double> nearest_ft = soltab_ft.GetValuesOrWeights(
      "val", "Antenna12", times, freqs, 0, 0, true);

  // GetValuesOrWeights should make sure that frequency is the fastest changing
  // index, even when in the underlying h5 array time was changing fastest
  BOOST_CHECK_EQUAL_COLLECTIONS(nearest_tf.begin(), nearest_tf.end(),
                                nearest_ft.begin(), nearest_ft.end());
}

BOOST_FIXTURE_TEST_CASE(times_freqs_generation, H5Fixture) {
  H5Parm h5parm("tH5Parm_tmp.h5", false, false, "sol000");

  SolTab soltab = h5parm.GetSolTab("timefreq");

  const std::vector<double> freqs{130.0e6, 135.0e6, 131.0e6};
  const std::vector<double> times{57000.0, 57878.7, 57880.3, 57890.3, 58000.0};

  TimesAndFrequencies times_and_frequencies =
      soltab.GetTimesAndFrequencies(times, freqs, 0, 0, true);

  BOOST_CHECK_EQUAL(times_and_frequencies.time_axis.size(), 7);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[0], 57878.5, 1.0e-8);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[1], 57880.5, 1.0e-8);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[2], 57882.5, 1.0e-8);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[3], 57884.5, 1.0e-8);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[4], 57886.5, 1.0e-8);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[5], 57888.5, 1.0e-8);
  BOOST_CHECK_CLOSE(times_and_frequencies.time_axis[6], 57890.5, 1.0e-8);
  BOOST_CHECK_EQUAL(times_and_frequencies.start_time_index, 0);
  BOOST_CHECK_EQUAL(times_and_frequencies.num_times, 7);

  BOOST_CHECK_EQUAL(times_and_frequencies.freq_axis.size(), 2);
  BOOST_CHECK_CLOSE(times_and_frequencies.freq_axis[0], 1.3e+08, 1.0e-12);
  BOOST_CHECK_CLOSE(times_and_frequencies.freq_axis[1], 1.31e+08, 1.0e-12);
  BOOST_CHECK_EQUAL(times_and_frequencies.start_freq_index, 0);
  BOOST_CHECK_EQUAL(times_and_frequencies.num_freqs, 2);
}

/**
 * It is not possible to replace an arbitrary soltab by H5Cache.
 * A soltab to be cacheable must have at least these axes:
 * ant, time, freq, dir, and pol, and their sizes multiplied
 * must be the size of the val and weight datasets. Therefore,
 * caching "completeset" is possible, but other soltab in the
 * H5 file are not.
 */

BOOST_FIXTURE_TEST_CASE(parameters_with_cache, H5Fixture) {
  std::vector<std::string> tables_in_cache = {"completeset"};
  H5Parm h5parm("tH5Parm_tmp.h5", tables_in_cache);

  SolTab& soltab = h5parm.GetSolTab("completeset");
  BOOST_CHECK(dynamic_cast<H5Cache*>(&soltab));
}

BOOST_FIXTURE_TEST_CASE(parameters_without_cache, H5Fixture) {
  std::vector<std::string> tables_in_cache = {"completeset"};
  H5Parm h5parm("tH5Parm_tmp.h5", tables_in_cache);

  SolTab& soltab = h5parm.GetSolTab("freqtime");
  BOOST_CHECK(dynamic_cast<H5Cache*>(&soltab) == nullptr);
}

BOOST_FIXTURE_TEST_CASE(caching_not_possible, H5Fixture) {
  std::vector<std::string> tables_in_cache = {"freqtime"};
  BOOST_CHECK_THROW(H5Parm h5parm("tH5Parm_tmp.h5", tables_in_cache),
                    std::runtime_error);
}

BOOST_FIXTURE_TEST_CASE(table_does_not_exist, H5Fixture) {
  std::vector<std::string> tables_in_cache = {"doesnotexist"};
  BOOST_CHECK_THROW(H5Parm h5parm("tH5Parm_tmp.h5", tables_in_cache),
                    std::runtime_error);
}

/**
 * Tests to verify the data integrity.
 * H5Cache is capable to ingest all H5 parameters stored
 * in any order in a file. The information must be preserved
 * no matter the order in which it is stored.
 */

BOOST_FIXTURE_TEST_CASE(canonical_order_cache, H5Fixture) {
  std::vector<std::string> tables_in_cache = {"canonical"};
  H5Parm h5parm("tH5Parm_tmp.h5", tables_in_cache);

  SolTab& soltab = h5parm.GetSolTab("canonical");
  BOOST_CHECK(dynamic_cast<H5Cache*>(&soltab));

  std::string antenna_name = "Tile11";

  std::vector<double> times(kNumTimes);
  std::iota(times.begin(), times.end(), 57800.0);

  std::vector<double> frequencies(kNumFrequencies);
  std::iota(frequencies.begin(), frequencies.end(), 130.0e6);

  size_t antenna = 0;
  size_t direction = 3;
  size_t polarization = 1;

  std::vector<double> sub_array = soltab.GetValues(
      antenna_name, times, frequencies, polarization, direction, true);

  xt::xtensor<double, 2> values;
  const std::array<size_t, 2> shape = {kNumTimes, kNumFrequencies};
  values.resize(shape);
  std::copy(sub_array.begin(), sub_array.end(), values.begin());

  // Recovering the data originally stored in the H5 file.
  // We verify that indices matches the ones encoded in each value
  // of the array. On success, this means that data integrity is
  // preserved.
  for (size_t time = 0; time < kNumTimes; ++time) {
    for (size_t frequency = 0; frequency < kNumFrequencies; ++frequency) {
      const double expected = time * 10000 + frequency * 1000 + antenna * 100 +
                              direction * 10 + polarization;
      const double actual = values(time, frequency);
      BOOST_CHECK_CLOSE(expected, actual, 1.0e-8);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(non_canonical_order_cache, H5Fixture) {
  std::vector<std::string> tables_in_cache = {"noncanonical"};
  H5Parm h5parm("tH5Parm_tmp.h5", tables_in_cache);

  SolTab& soltab = h5parm.GetSolTab("noncanonical");
  BOOST_CHECK(dynamic_cast<H5Cache*>(&soltab));

  std::string antenna_name = "Tile11";

  std::vector<double> times(kNumTimes);
  std::iota(times.begin(), times.end(), 57800.0);

  std::vector<double> frequencies(kNumFrequencies);
  std::iota(frequencies.begin(), frequencies.end(), 130.0e6);

  size_t antenna = 0;
  size_t direction = 3;
  size_t polarization = 1;

  std::vector<double> sub_array = soltab.GetValues(
      antenna_name, times, frequencies, polarization, direction, true);

  xt::xtensor<double, 2> values;
  const std::array<size_t, 2> shape = {kNumTimes, kNumFrequencies};
  values.resize(shape);
  std::copy(sub_array.begin(), sub_array.end(), values.begin());

  // Recovering the data originally stored in the H5 file.
  // We verify that indices matches the ones encoded in each value
  // of the array. On success, this means that data integrity is
  // preserved not matter it was stored with a different order
  // from the canonical one.
  for (size_t time = 0; time < kNumTimes; ++time) {
    for (size_t frequency = 0; frequency < kNumFrequencies; ++frequency) {
      const double expected = time * 10000 + frequency * 1000 + antenna * 100 +
                              direction * 10 + polarization;
      const double actual = values(time, frequency);
      BOOST_CHECK_CLOSE(expected, actual, 1.0e-8);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END()
