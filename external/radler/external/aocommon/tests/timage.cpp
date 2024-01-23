#include <boost/test/unit_test.hpp>

#include <aocommon/image.h>

#include <aocommon/uvector.h>

#include <random>

using aocommon::DImage;
using aocommon::Image;
using aocommon::ImageBase;

BOOST_AUTO_TEST_SUITE(image)

struct ImageFixture {
  ImageFixture()
      : image(3, 2, {6.0, 5.0, 4.0, 3.0, 2.0, 1.0}), image_ref(image) {}

  Image image;
  const Image image_ref;
};

/// Generates an image with (y+1)*100 + (x+1) at position (x,y).
static Image CreateTestImage(int width, int height) {
  Image image(width, height);
  for (int y = 0; y < height; ++y) {
    for (int x = 0; x < width; ++x) {
      image[y * width + x] = (y + 1) * 100 + (x + 1);
    }
  }
  return image;
}

BOOST_AUTO_TEST_CASE(equality) {
  BOOST_CHECK(Image() == Image());
  BOOST_CHECK(!(Image() == Image(1, 1, 0.0)));
  BOOST_CHECK(Image(1, 2, 2.0) == Image(1, 2, 2.0));
  BOOST_CHECK(!(Image(1, 2, 2.0) == Image(2, 1, 2.0)));
  BOOST_CHECK(!(Image(1, 2, 2.0) == Image(1, 2, 3.0)));
  Image image = CreateTestImage(8, 8);
  *(image.end() - 1) = -1.0;
  BOOST_CHECK(!(image == CreateTestImage(8, 8)));
}

BOOST_AUTO_TEST_CASE(inequality) {
  BOOST_CHECK(!(Image() != Image()));
  BOOST_CHECK(Image() != Image(1, 1, 0.0));
  BOOST_CHECK(!(Image(1, 2, 2.0) != Image(1, 2, 2.0)));
  BOOST_CHECK(Image(1, 2, 2.0) != Image(2, 1, 2.0));
  BOOST_CHECK(Image(1, 2, 2.0) != Image(1, 2, 3.0));
}

BOOST_AUTO_TEST_CASE(median_empty) {
  BOOST_CHECK_EQUAL(Image::Median(nullptr, 0), 0.0f);
}

BOOST_AUTO_TEST_CASE(median_single) {
  aocommon::UVector<float> arr(1, 1.0);
  BOOST_CHECK_EQUAL(Image::Median(arr.data(), arr.size()), 1.0f);

  arr[0] = std::numeric_limits<float>::quiet_NaN();
  Image::Median(arr.data(),
                arr.size());  // undefined -- just make sure it doesn't crash
}

BOOST_AUTO_TEST_CASE(median_two_elements) {
  {
    aocommon::UVector<float> arr(2, 1.0);
    BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0,
                               1e-5);
  }

  {
    aocommon::UVector<float> arr(2, 0.0);
    arr[1] = 2.0;
    BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0,
                               1e-5);
  }

  {
    aocommon::UVector<float> arr(2, 1.0);
    arr[1] = -1.0;
    BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 0.0,
                               1e-5);
  }

  {
    aocommon::UVector<float> arr(2, 13.0);
    arr[1] = std::numeric_limits<float>::quiet_NaN();
    BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 13.0,
                               1e-5);
  }
}

BOOST_AUTO_TEST_CASE(median_three_elements) {
  aocommon::UVector<float> arr(3, 1.0);
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 0.0;
  arr[1] = 1.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 3.0;
  arr[1] = -3.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 2.0, 1e-5);

  arr[1] = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 2.5, 1e-5);

  arr[0] = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 2.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(median_double) {
  aocommon::UVector<double> arr = {6.0, -3.0, 2.0};
  BOOST_CHECK_CLOSE_FRACTION(DImage::Median(arr.data(), arr.size()), 2.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(median_invalid) {
  aocommon::UVector<std::complex<float>> arr = {6.0, -3.0, 2.0};
  BOOST_CHECK_THROW(
      ImageBase<std::complex<float>>::Median(arr.data(), arr.size()),
      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(mad_empty) {
  BOOST_CHECK_EQUAL(Image::MAD(nullptr, 0), 0.0);
}

BOOST_AUTO_TEST_CASE(mad_single) {
  aocommon::UVector<float> arr(1, 1.0);
  BOOST_CHECK_EQUAL(Image::MAD(arr.data(), arr.size()), 0.0);
}

BOOST_AUTO_TEST_CASE(mad_two_elements) {
  aocommon::UVector<float> arr(2, 1.0);
  BOOST_CHECK_EQUAL(Image::MAD(arr.data(), arr.size()), 0.0);

  arr[0] = 0.0;
  arr[1] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 1.0;
  arr[1] = -1.0;
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 0.0, 1e-5);

  arr[0] = 13.0;
  arr[1] = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(Image::Median(arr.data(), arr.size()), 13.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(mad_three_elements) {
  aocommon::UVector<float> arr(3, 1.0);
  BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 0.0, 1e-5);

  arr[0] = 0.0;
  arr[1] = 1.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 1.0, 1e-5);

  arr[0] = 3.0;
  arr[1] = -3.0;
  arr[2] = 2.0;
  BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 1.0, 1e-5);

  arr[1] = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 0.5, 1e-5);

  arr[0] = std::numeric_limits<float>::quiet_NaN();
  BOOST_CHECK_CLOSE_FRACTION(Image::MAD(arr.data(), arr.size()), 0.0, 1e-5);
}

BOOST_AUTO_TEST_CASE(stddev_from_mad) {
  std::mt19937 rnd;
  std::normal_distribution<float> dist(1.0f, 5.0f);
  aocommon::UVector<float> data(10000);
  for (size_t i = 0; i != data.size(); ++i) data[i] = dist(rnd);
  BOOST_CHECK_CLOSE_FRACTION(Image::StdDevFromMAD(data.data(), data.size()),
                             5.0f, 0.05);
}

BOOST_AUTO_TEST_CASE(uninitialized_constructor) {
  const size_t kWidth = 10;
  const size_t kHeight = 20;
  Image image(kWidth, kHeight);
  BOOST_CHECK_EQUAL(image.Width(), kWidth);
  BOOST_CHECK_EQUAL(image.Height(), kHeight);
  BOOST_CHECK_EQUAL(image.Size(), kWidth * kHeight);
}

BOOST_AUTO_TEST_CASE(initializing_constructors) {
  const size_t kWidth = 10;
  const size_t kHeight = 20;
  Image image(kWidth, kHeight, 42.0);
  aocommon::UVector<float> ref(kWidth * kHeight, 42.0);
  BOOST_CHECK_EQUAL_COLLECTIONS(image.begin(), image.end(), ref.begin(),
                                ref.end());

  // Copy constructor
  Image image1(image);
  BOOST_CHECK_EQUAL_COLLECTIONS(image1.begin(), image1.end(), image.begin(),
                                image.end());

  // Move constructor
  Image image2(std::move(image1));
  BOOST_CHECK_EQUAL_COLLECTIONS(image2.begin(), image2.end(), image.begin(),
                                image.end());
  BOOST_CHECK(image1.Empty());
  BOOST_CHECK(image1.Data() == nullptr);

  // reset
  image2.Reset();
  BOOST_CHECK(image2.Empty());
  BOOST_CHECK(image2.Data() == nullptr);
}

BOOST_AUTO_TEST_CASE(non_owning_constructors) {
  const size_t kWidth = 10;
  const size_t kHeight = 20;
  std::vector<float> data(kWidth * kHeight, 42.0);
  Image image1(data.data(), kWidth, kHeight);
  BOOST_CHECK(!image1.Empty());
  BOOST_CHECK_EQUAL(image1.Width(), kWidth);
  BOOST_CHECK_EQUAL(image1.Height(), kHeight);
  // Image1 should borrow the data:
  BOOST_CHECK_EQUAL(image1.Data(), data.data());

  std::vector<float> ref(kWidth * kHeight, 42.0);
  BOOST_CHECK_EQUAL_COLLECTIONS(image1.begin(), image1.end(), ref.begin(),
                                ref.end());

  // Copy constructor from non-owning image
  Image image2(image1);
  BOOST_CHECK_EQUAL_COLLECTIONS(image2.begin(), image2.end(), image1.begin(),
                                image1.end());
  // Image2 should have its own copy:
  BOOST_CHECK_NE(image2.Data(), data.data());
  BOOST_CHECK_EQUAL_COLLECTIONS(image2.begin(), image2.end(), ref.begin(),
                                ref.end());

  // Move constructor from non-owning image
  Image image3(std::move(image1));
  BOOST_CHECK_EQUAL_COLLECTIONS(image3.begin(), image3.end(), ref.begin(),
                                ref.end());
  BOOST_CHECK_EQUAL(image3.Data(), data.data());
  BOOST_CHECK(image1.Empty());
  BOOST_CHECK(image1.Data() == nullptr);

  BOOST_CHECK_EQUAL_COLLECTIONS(image3.begin(), image3.end(), ref.begin(),
                                ref.end());

  // reset
  image3.Reset();
  BOOST_CHECK(image3.Empty());
  BOOST_CHECK(image3.Data() == nullptr);
}

BOOST_AUTO_TEST_CASE(initializer_list_constructor) {
  Image image(2, 4, {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f, 7.0f, 8.0f});
  BOOST_REQUIRE_EQUAL(image.Width(), 2);
  BOOST_REQUIRE_EQUAL(image.Height(), 4);
  for (size_t i = 0; i != 8; ++i) {
    BOOST_CHECK_CLOSE(image[i], i + 1.0f, 1e-6f);
  }
}

BOOST_FIXTURE_TEST_CASE(indexing_operator, ImageFixture) {
  const aocommon::UVector<float> arr = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_EQUAL(image[i], arr[i]);
  }
}

BOOST_FIXTURE_TEST_CASE(add_assign, ImageFixture) {
  image += image;
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], image_ref[i] + image_ref[i], 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(multiply_assign, ImageFixture) {
  image *= image;
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], image_ref[i] * image_ref[i], 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(subtract_assign, ImageFixture) {
  image -= image;
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], 0.0f, 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(divide_assign, ImageFixture) {
  const float factor = 2.0;
  image /= factor;
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], image_ref[i] / factor, 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(multiply_assign_factor, ImageFixture) {
  const float factor = 2.0;
  image *= factor;
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], factor * image_ref[i], 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(sqrt, ImageFixture) {
  image.Sqrt();
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], std::sqrt(image_ref[i]), 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(square, ImageFixture) {
  image.Square();
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], image_ref[i] * image_ref[i], 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(add_with_factor, ImageFixture) {
  image.AddWithFactor(image, 4.0f);
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], 5.0f * image_ref[i], 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(square_with_factor, ImageFixture) {
  const float factor = 2.0;
  image.SquareWithFactor(factor);
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], image_ref[i] * image_ref[i] * factor, 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(sqrt_with_factor, ImageFixture) {
  const float factor = 2.0;
  image.SqrtWithFactor(factor);
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], std::sqrt(image_ref[i]) * factor, 1e-4);
  }
}

BOOST_FIXTURE_TEST_CASE(negate, ImageFixture) {
  image.Negate();
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_EQUAL(image[i], -image_ref[i]);
  }
}

BOOST_AUTO_TEST_CASE(trim_and_untrim) {
  const size_t width = 5;
  const size_t height = 3;
  aocommon::UVector<float> row_values = {0., 1., 2., 3., 4.};
  aocommon::UVector<float> image_data;
  for (size_t i = 0; i != height; ++i) {
    image_data.insert(image_data.begin() + i * width, row_values.begin(),
                      row_values.end());
  }

  Image image(width, height);
  image.Assign(image_data.begin(), image_data.end());

  const size_t x1 = 1;
  const size_t y1 = 1;
  const size_t box_width = 3;
  const size_t box_height = 1;

  Image trimmed_image = image.TrimBox(x1, y1, box_width, box_height);
  // Offsets are only correct if box_height == 1
  BOOST_CHECK_EQUAL_COLLECTIONS(trimmed_image.begin(), trimmed_image.end(),
                                row_values.begin() + x1,
                                row_values.begin() + x1 + box_width);

  Image untrimmed_image = trimmed_image.Untrim(width, height);
  BOOST_CHECK_EQUAL(untrimmed_image.Size(), image.Size());
  BOOST_CHECK_EQUAL(
      std::accumulate(untrimmed_image.begin(), untrimmed_image.end(), 0.0f),
      std::accumulate(trimmed_image.begin(), trimmed_image.end(), 0.0f));
}

BOOST_AUTO_TEST_CASE(resize_equal_size) {
  const size_t kWidth = 3;
  const size_t kHeight = 5;
  const Image image = CreateTestImage(kWidth, kHeight);
  const Image resized = image.Resize(kWidth, kHeight);
  BOOST_REQUIRE_EQUAL(resized.Width(), kWidth);
  BOOST_REQUIRE_EQUAL(resized.Height(), kHeight);
  BOOST_CHECK_EQUAL_COLLECTIONS(image.begin(), image.end(), resized.begin(),
                                resized.end());
}

BOOST_AUTO_TEST_CASE(resize_odd_less_width_more_height) {
  const size_t kWidth = 4;
  const size_t kHeight = 4;
  const Image image = CreateTestImage(kWidth, kHeight);

  const Image resized = image.Resize(kWidth - 1, kHeight + 3);
  BOOST_REQUIRE_EQUAL(resized.Width(), kWidth - 1);
  BOOST_REQUIRE_EQUAL(resized.Height(), kHeight + 3);
  std::vector<float> kExpectedData{
      0,   0,   0,    // one added row at the top
      101, 102, 103,  // existing row 0, last item removed
      201, 202, 203,  // existing row 1, last item removed
      301, 302, 303,  // existing row 2, last item removed
      401, 402, 403,  // existing row 3, last item removed
      0,   0,   0,    // first added row at the bottom
      0,   0,   0     // second added row at the bottom
  };
  BOOST_CHECK_EQUAL_COLLECTIONS(kExpectedData.begin(), kExpectedData.end(),
                                resized.begin(), resized.end());
}

BOOST_AUTO_TEST_CASE(resize_odd_more_width_less_height) {
  const size_t kWidth = 4;
  const size_t kHeight = 4;
  const Image image = CreateTestImage(kWidth, kHeight);
  const Image resized = image.Resize(kWidth + 3, kHeight - 1);
  BOOST_REQUIRE_EQUAL(resized.Width(), kWidth + 3);
  BOOST_REQUIRE_EQUAL(resized.Height(), kHeight - 1);
  std::vector<float> kExpectedData{
      0, 101, 102, 103, 104, 0, 0,  // existing row 0 with padding
      0, 201, 202, 203, 204, 0, 0,  // existing row 1 with padding
      0, 301, 302, 303, 304, 0, 0   // existing row 2 with padding
  };
  BOOST_CHECK_EQUAL_COLLECTIONS(kExpectedData.begin(), kExpectedData.end(),
                                resized.begin(), resized.end());
}

BOOST_AUTO_TEST_CASE(resize_even_enlarge_with_fill) {
  const size_t kWidth = 2;
  const size_t kHeight = 2;
  const Image image = CreateTestImage(kWidth, kHeight);
  const Image resized = image.Resize(kWidth + 4, kHeight + 4, 42);
  BOOST_REQUIRE_EQUAL(resized.Width(), kWidth + 4);
  BOOST_REQUIRE_EQUAL(resized.Height(), kHeight + 4);
  std::vector<float> kExpectedData{
      42, 42, 42,  42,  42, 42,  // first added top row
      42, 42, 42,  42,  42, 42,  // second added top row
      42, 42, 101, 102, 42, 42,  // existing row 0 with padding
      42, 42, 201, 202, 42, 42,  // existing row 1 with padding
      42, 42, 42,  42,  42, 42,  // first added bottom row
      42, 42, 42,  42,  42, 42,  // second added bottom row
  };
  BOOST_CHECK_EQUAL_COLLECTIONS(kExpectedData.begin(), kExpectedData.end(),
                                resized.begin(), resized.end());
}

BOOST_AUTO_TEST_CASE(resize_even_reduce) {
  const size_t kWidth = 10;
  const size_t kHeight = 10;
  const Image image = CreateTestImage(kWidth, kHeight);
  const Image resized = image.Resize(kWidth - 4, kHeight - 6);
  BOOST_REQUIRE_EQUAL(resized.Width(), kWidth - 4);
  BOOST_REQUIRE_EQUAL(resized.Height(), kHeight - 6);
  std::vector<float> kExpectedData{
      403, 404, 405, 406, 407, 408,  // Remainder of fourth row.
      503, 504, 505, 506, 507, 508,  // Remainder of fifth row.
      603, 604, 605, 606, 607, 608,  // Remainder of sixth row.
      703, 704, 705, 706, 707, 708,  // Remainder of seventh row.
  };
  BOOST_CHECK_EQUAL_COLLECTIONS(kExpectedData.begin(), kExpectedData.end(),
                                resized.begin(), resized.end());
}

BOOST_AUTO_TEST_CASE(resize_to_nothing) {
  const size_t kWidth = 4;
  const size_t kHeight = 4;
  const Image image = CreateTestImage(kWidth, kHeight);
  const Image resized_0_0 = image.Resize(0, 0);
  BOOST_CHECK_EQUAL(resized_0_0.Width(), 0);
  BOOST_CHECK_EQUAL(resized_0_0.Height(), 0);

  const Image resized_0_1 = image.Resize(0, 1);
  BOOST_CHECK_EQUAL(resized_0_1.Width(), 0);
  BOOST_CHECK_EQUAL(resized_0_1.Height(), 1);

  const Image resized_1_0 = image.Resize(1, 0);
  BOOST_CHECK_EQUAL(resized_1_0.Width(), 1);
  BOOST_CHECK_EQUAL(resized_1_0.Height(), 0);
}

BOOST_AUTO_TEST_CASE(pad) {
  const size_t kWidth = 2;
  const size_t kHeight = 2;
  const Image image = CreateTestImage(kWidth, kHeight);
  BOOST_CHECK(image.Pad(0, 0, 0, 0) == image);

  const size_t kLeft = 1;
  const size_t kTop = 2;
  const size_t kRight = 3;
  const size_t kBottom = 4;
  const Image padded = image.Pad(kLeft, kTop, kRight, kBottom);
  const Image reference = Image(
      kWidth + kLeft + kRight, kHeight + kTop + kBottom,
      {
          0.0, 0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0, 0.0, 0.0,
          0.0, 101.0, 102.0, 0.0, 0.0, 0.0, 0.0, 201.0, 202.0, 0.0, 0.0, 0.0,
          0.0, 0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0, 0.0, 0.0,
          0.0, 0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0,   0.0,   0.0, 0.0, 0.0,
      });
  BOOST_REQUIRE_EQUAL(padded.Width(), reference.Width());
  BOOST_REQUIRE_EQUAL(padded.Height(), reference.Height());
  for (size_t i = 0; i != padded.Width() * padded.Height(); ++i) {
    BOOST_CHECK_CLOSE_FRACTION(padded[i], reference[i], 1e-6);
  }
  BOOST_CHECK(padded.TrimBox(kLeft, kTop, kWidth, kHeight) == image);
}

BOOST_AUTO_TEST_CASE(serialization) {
  const size_t width = 3;
  const size_t height = 2;
  aocommon::UVector<float> arr = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
  Image image_in(width, height);
  image_in.Assign(arr.begin(), arr.end());

  aocommon::SerialOStream ostr;
  image_in.Serialize(ostr);

  aocommon::SerialIStream istr(std::move(ostr));
  Image image_out;
  image_out.Unserialize(istr);

  BOOST_CHECK_EQUAL(image_in.Width(), image_out.Width());
  BOOST_CHECK_EQUAL(image_in.Height(), image_out.Height());
  BOOST_CHECK_EQUAL(image_in.Size(), image_out.Size());
  for (size_t i = 0; i != arr.size(); ++i) {
    BOOST_CHECK_EQUAL(image_out[i], arr[i]);
  }
}

BOOST_AUTO_TEST_CASE(add_sub_image) {
  Image image(4, 5, {10.0f, 11.0f, 12.0f, 13.0f, 14.0f, 15.0f, 16.0f,
                     17.0f, 18.0f, 19.0f, 20.0f, 21.0f, 22.0f, 23.0f,
                     24.0f, 25.0f, 26.0f, 27.0f, 28.0f, 29.0f});
  const Image sub_image(2, 2, {100.0f, 200.0f, 300.0f, 400.0f});
  image.AddSubImage(image.Data(), 1, 1, image.Width(), sub_image.Data(),
                    sub_image.Width(), sub_image.Height());
  const Image reference(4, 5,
                        {10.0f, 11.0f, 12.0f,  13.0f,  14.0f, 115.0f, 216.0f,
                         17.0f, 18.0f, 319.0f, 420.0f, 21.0f, 22.0f,  23.0f,
                         24.0f, 25.0f, 26.0f,  27.0f,  28.0f, 29.0f});
  BOOST_REQUIRE_EQUAL(image.Size(), reference.Size());
  for (size_t i = 0; i != image.Size(); ++i) {
    BOOST_CHECK_CLOSE(image[i], reference[i], 1e-6);
  }
}

BOOST_AUTO_TEST_CASE(remove_nans) {
  const Image reference_a(2, 3, {1.0f, 2.0f, 3.0f, 4.0f, 5.0f, 6.0f});

  Image test_image = reference_a;
  test_image.RemoveNans();
  BOOST_CHECK_EQUAL_COLLECTIONS(reference_a.begin(), reference_a.end(),
                                test_image.begin(), test_image.end());

  test_image[0] = std::numeric_limits<float>::infinity();
  test_image[3] = std::numeric_limits<float>::quiet_NaN();
  test_image[5] = -std::numeric_limits<float>::infinity();
  test_image.RemoveNans();
  const Image reference_b(2, 3, {0.0f, 2.0f, 3.0f, 0.0f, 5.0f, 0.0f});
  BOOST_CHECK_EQUAL_COLLECTIONS(reference_b.begin(), reference_b.end(),
                                test_image.begin(), test_image.end());
}

BOOST_AUTO_TEST_SUITE_END()
