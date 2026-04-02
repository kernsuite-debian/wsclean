#include <boost/test/unit_test.hpp>

// See explanation in tbanddata.cpp
#define DISABLE_CASACORE_IN_BANDDATA
#include <aocommon/multibanddata.h>

#include <vector>

using aocommon::BandData;
using aocommon::ChannelInfo;
using aocommon::MultiBandData;

namespace {
// Band 1 has (purposely):
// - A higher frequency than Band 2
// - Fewer channels than Band 2
// - Have a different channel width
// MultiBandData should be able to handle this.
const std::vector<ChannelInfo> kChannels1{ChannelInfo(180e6, 10e6),
                                          ChannelInfo(190e6, 10e6)};
const std::vector<ChannelInfo> kChannels2{
    ChannelInfo(140e6, 5e6), ChannelInfo(145e6, 5e6), ChannelInfo(150e6, 5e6)};
constexpr size_t kDataDescId1 = 0;
constexpr size_t kDataDescId2 = 3;
constexpr double kReferenceFrequency1 = 185e6;
constexpr double kReferenceFrequency2 = 145e6;

void CheckDataDescIds(const MultiBandData& multi_band) {
  std::vector<size_t> data_desc_ids;
  for (size_t d : multi_band.DataDescIds()) {
    data_desc_ids.emplace_back(d);
  }
  BOOST_CHECK_EQUAL(data_desc_ids.size(), 2);
  BOOST_CHECK_EQUAL(data_desc_ids[0], kDataDescId1);
  BOOST_CHECK_EQUAL(data_desc_ids[1], kDataDescId2);
}

}  // namespace

BOOST_AUTO_TEST_SUITE(multi_band_data)

BOOST_AUTO_TEST_CASE(empty) {
  MultiBandData multiBand;
  BOOST_CHECK_EQUAL(multiBand.BandCount(), 0u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandStart(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandEnd(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.Bandwidth(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.CentreFrequency(), 0.0, 1e-6);
  BOOST_CHECK_EQUAL(multiBand.HighestDataDescId(), 0u);
  BOOST_CHECK_EQUAL(multiBand.HighestBandId(), 0u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.HighestFrequency(), 0.0, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.LowestFrequency(), 0.0, 1e-6);
  BOOST_CHECK_EQUAL(multiBand.MaxBandChannels(), 0u);
  BOOST_CHECK(!multiBand.HasDataDescId(0));
  BOOST_CHECK(multiBand.begin() == multiBand.end());
  BOOST_CHECK(multiBand.DataDescIds().begin() == multiBand.DataDescIds().end());

  const MultiBandData moved(std::move(multiBand));
  BOOST_CHECK_EQUAL(moved.BandCount(), 0u);

  const MultiBandData copied(multiBand);
  BOOST_CHECK_EQUAL(copied.BandCount(), 0u);

  aocommon::SerialOStream stream;
  copied.Serialize(stream);
  aocommon::SerialIStream input_stream(std::move(stream));
  MultiBandData serialized;
  serialized.Unserialize(input_stream);
  BOOST_CHECK_EQUAL(serialized.BandCount(), 0u);
}

BOOST_AUTO_TEST_CASE(irregular_bands) {
  MultiBandData multiBand;
  multiBand.SetBand(kDataDescId1, BandData(kChannels1, kReferenceFrequency1));
  multiBand.SetBand(kDataDescId2, BandData(kChannels2, kReferenceFrequency2));

  aocommon::SerialOStream stream;
  multiBand.Serialize(stream);
  aocommon::SerialIStream input_stream(std::move(stream));
  multiBand = MultiBandData();
  BOOST_CHECK_EQUAL(multiBand.BandCount(), 0u);
  multiBand.Unserialize(input_stream);

  BOOST_CHECK_EQUAL(multiBand.BandCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandStart(), 137.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.BandEnd(), 195e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.Bandwidth(), 195e6 - 137.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.CentreFrequency(),
                             0.5 * (195e6 + 137.5e6), 1e-6);
  BOOST_CHECK_EQUAL(multiBand.HighestDataDescId(), kDataDescId2);
  BOOST_CHECK_EQUAL(multiBand.HighestBandId(), 1u);
  BOOST_CHECK_EQUAL(multiBand.GetBandIndex(kDataDescId1), 0u);
  BOOST_CHECK_EQUAL(multiBand.GetBandIndex(kDataDescId2), 1u);
  BOOST_CHECK_EQUAL(multiBand.MaxBandChannels(), 3u);
  BOOST_CHECK_EQUAL(multiBand.DataDescIdWithMaxChannels(), kDataDescId2);
  BOOST_CHECK(multiBand.HasDataDescId(0));
  BOOST_CHECK(!multiBand.HasDataDescId(1));
  BOOST_CHECK(!multiBand.HasDataDescId(2));
  BOOST_CHECK(multiBand.HasDataDescId(3));
  BOOST_CHECK(!multiBand.HasDataDescId(4));
  BOOST_CHECK_CLOSE_FRACTION(multiBand.HighestFrequency(), 190e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(multiBand.LowestFrequency(), 140e6, 1e-6);
  BOOST_CHECK_EQUAL(multiBand[kDataDescId1].ChannelCount(), 2u);
  BOOST_CHECK_EQUAL(multiBand[kDataDescId2].ChannelCount(), 3u);

  CheckDataDescIds(multiBand);
}

BOOST_AUTO_TEST_CASE(partial_band_a) {
  MultiBandData multiBand;
  multiBand.SetBand(kDataDescId1, BandData(kChannels1, kReferenceFrequency1));
  multiBand.SetBand(kDataDescId2, BandData(kChannels2, kReferenceFrequency2));
  MultiBandData partial_band(multiBand, 1, 2);
  BOOST_CHECK_EQUAL(partial_band.BandCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.BandStart(), 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.BandEnd(), 195e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.Bandwidth(), 195e6 - 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.CentreFrequency(),
                             0.5 * (195e6 + 142.5e6), 1e-6);
  BOOST_CHECK_EQUAL(partial_band.HighestDataDescId(), kDataDescId2);
  BOOST_CHECK_EQUAL(partial_band.HighestBandId(), 1u);
  BOOST_CHECK_EQUAL(partial_band.GetBandIndex(kDataDescId1), 0u);
  BOOST_CHECK_EQUAL(partial_band.GetBandIndex(kDataDescId2), 1u);
  BOOST_CHECK_EQUAL(partial_band.MaxBandChannels(), 1u);
  BOOST_CHECK_EQUAL(partial_band.DataDescIdWithMaxChannels(), kDataDescId1);
  BOOST_CHECK(partial_band.HasDataDescId(0));
  BOOST_CHECK(!partial_band.HasDataDescId(1));
  BOOST_CHECK(!partial_band.HasDataDescId(2));
  BOOST_CHECK(partial_band.HasDataDescId(3));
  BOOST_CHECK(!partial_band.HasDataDescId(4));
  BOOST_CHECK_CLOSE_FRACTION(partial_band.HighestFrequency(), 190e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.LowestFrequency(), 145e6, 1e-6);
  BOOST_CHECK_EQUAL(partial_band[kDataDescId1].ChannelCount(), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partial_band[kDataDescId1].ChannelFrequency(0),
                             190e6, 1e-6);
  BOOST_CHECK_EQUAL(partial_band[kDataDescId2].ChannelCount(), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partial_band[kDataDescId2].ChannelFrequency(0),
                             145e6, 1e-6);

  CheckDataDescIds(partial_band);
}

BOOST_AUTO_TEST_CASE(partial_band_b) {
  MultiBandData multiBand;
  multiBand.SetBand(kDataDescId1, BandData(kChannels1, kReferenceFrequency1));
  multiBand.SetBand(kDataDescId2, BandData(kChannels2, kReferenceFrequency2));
  MultiBandData partial_band(multiBand, 1, 3);
  BOOST_CHECK_EQUAL(partial_band.BandCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.BandStart(), 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.BandEnd(), 195e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.Bandwidth(), 195e6 - 142.5e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.CentreFrequency(),
                             0.5 * (195e6 + 142.5e6), 1e-6);
  BOOST_CHECK_EQUAL(partial_band.HighestDataDescId(), kDataDescId2);
  BOOST_CHECK_EQUAL(partial_band.HighestBandId(), 1u);
  BOOST_CHECK_EQUAL(partial_band.GetBandIndex(kDataDescId1), 0u);
  BOOST_CHECK_EQUAL(partial_band.GetBandIndex(kDataDescId2), 1u);
  BOOST_CHECK_EQUAL(partial_band.MaxBandChannels(), 2u);
  BOOST_CHECK_EQUAL(partial_band.DataDescIdWithMaxChannels(), kDataDescId2);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.HighestFrequency(), 190e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band.LowestFrequency(), 145e6, 1e-6);
  BOOST_CHECK_EQUAL(partial_band[kDataDescId1].ChannelCount(), 1u);
  BOOST_CHECK_CLOSE_FRACTION(partial_band[kDataDescId1].ChannelFrequency(0),
                             190e6, 1e-6);
  BOOST_CHECK_EQUAL(partial_band[kDataDescId2].ChannelCount(), 2u);
  BOOST_CHECK_CLOSE_FRACTION(partial_band[kDataDescId2].ChannelFrequency(0),
                             145e6, 1e-6);
  BOOST_CHECK_CLOSE_FRACTION(partial_band[kDataDescId2].ChannelFrequency(1),
                             150e6, 1e-6);

  CheckDataDescIds(partial_band);
}

BOOST_AUTO_TEST_CASE(copy_and_move) {
  MultiBandData multi_band;
  multi_band.SetBand(kDataDescId1, BandData(kChannels1, kReferenceFrequency1));
  multi_band.SetBand(kDataDescId2, BandData(kChannels2, kReferenceFrequency2));

  MultiBandData copy(multi_band);
  BOOST_CHECK_EQUAL(copy.HighestDataDescId(), kDataDescId2);
  BOOST_CHECK_EQUAL(copy.BandCount(), 2u);
  CheckDataDescIds(copy);

  MultiBandData moved(std::move(multi_band));
  BOOST_CHECK_EQUAL(moved.HighestDataDescId(), kDataDescId2);
  BOOST_CHECK_EQUAL(moved.BandCount(), 2u);
  CheckDataDescIds(moved);
}

BOOST_AUTO_TEST_CASE(nonzero_first_data_desc_id) {
  MultiBandData multi_band;
  multi_band.SetBand(kDataDescId2, BandData(kChannels2, kReferenceFrequency2));
  BOOST_CHECK_EQUAL(multi_band.HighestDataDescId(), kDataDescId2);
  BOOST_CHECK_EQUAL(multi_band.BandCount(), 1u);
  BOOST_CHECK_EQUAL(multi_band.begin()->ChannelCount(), kChannels2.size());
  BOOST_CHECK_EQUAL(*multi_band.DataDescIds().begin(), kDataDescId2);
}

BOOST_AUTO_TEST_CASE(iterate) {
  MultiBandData multi_band;
  BOOST_CHECK(multi_band.begin() == multi_band.end());
  multi_band.SetBand(kDataDescId1, BandData(kChannels1, kReferenceFrequency1));
  BOOST_CHECK(++multi_band.begin() == multi_band.end());
  multi_band.SetBand(kDataDescId2, BandData(kChannels2, kReferenceFrequency2));
  BOOST_CHECK(++(++multi_band.begin()) == multi_band.end());
  constexpr auto check_function = [](const BandData& band, size_t& index) {
    BOOST_CHECK_LT(index, 2);
    if (index == 0) {
      BOOST_CHECK_EQUAL(band.ReferenceFrequency(), kReferenceFrequency1);
    } else {
      BOOST_CHECK_EQUAL(band.ReferenceFrequency(), kReferenceFrequency2);
    }
    ++index;
  };
  size_t index = 0;
  for (const BandData& band : multi_band) {
    check_function(band, index);
  }
  index = 0;
  for (BandData& band : multi_band) {
    check_function(band, index);
  }
  BOOST_CHECK_EQUAL(index, 2);
}

BOOST_AUTO_TEST_SUITE_END()
