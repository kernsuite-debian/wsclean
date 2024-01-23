#include <iostream>

#include <aocommon/lane.h>

#include <fstream>
#include <thread>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <boost/make_unique.hpp>

using aocommon::Lane;

// These are "first order" tests: none of these tests explicitly test
// multi-threading / thread-safety features of Lane

BOOST_AUTO_TEST_SUITE(lane)

BOOST_AUTO_TEST_CASE(construct_default) {
  const Lane<int> l;
  BOOST_CHECK_EQUAL(l.capacity(), 0);
  BOOST_CHECK_EQUAL(l.size(), 0);
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(construct_with_capacity) {
  const Lane<int> l(7);
  BOOST_CHECK_EQUAL(l.capacity(), 7);
  BOOST_CHECK_EQUAL(l.size(), 0);
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(read_write) {
  Lane<int> l(4);
  constexpr size_t n = 2;
  const int values[n] = {3, 4};
  l.write(values, n);
  l.write_end();
  BOOST_CHECK(!l.empty());
  BOOST_CHECK_EQUAL(l.size(), n);

  for (size_t i = 0; i != n; ++i) {
    int result;
    BOOST_REQUIRE(l.read(result));
    BOOST_CHECK_EQUAL(result, values[i]);
  }
  int result;
  BOOST_CHECK(l.empty());
  BOOST_CHECK(!l.read(result));
}

BOOST_AUTO_TEST_CASE(move_construct) {
  Lane<int> a(4);
  constexpr size_t n = 2;
  const int values[n] = {3, 4};
  a.write(values, n);

  Lane<int> b(std::move(a));
  BOOST_CHECK_EQUAL(a.capacity(), 0);
  BOOST_CHECK(a.empty());
  BOOST_CHECK_EQUAL(a.size(), 0);
  BOOST_CHECK_EQUAL(b.capacity(), 4);
  BOOST_CHECK(!b.empty());
  BOOST_CHECK_EQUAL(b.size(), n);
}

BOOST_AUTO_TEST_CASE(move_assign) {
  Lane<int> a(4);
  constexpr size_t n = 2;
  const int values_a[n] = {3, 4};
  a.write(values_a, n);

  Lane<int> b(3);
  const int values_b[3] = {5, 6, 7};
  b.write(values_b, 3);
  b = std::move(a);

  BOOST_CHECK_EQUAL(b.capacity(), 4);
  BOOST_CHECK(!b.empty());
  BOOST_CHECK_EQUAL(b.size(), n);
  for (size_t i = 0; i != n; ++i) {
    int result;
    BOOST_REQUIRE(b.read(result));
    BOOST_CHECK_EQUAL(result, values_a[i]);
  }
  BOOST_CHECK(b.empty());
}

BOOST_AUTO_TEST_CASE(emplace) {
  Lane<std::pair<int, int>> l(4);
  l.emplace(3, 30);
  l.emplace(4, 40);
  l.write_end();
  BOOST_CHECK(!l.empty());
  BOOST_CHECK_EQUAL(l.size(), 2);

  std::pair<int, int> result;
  BOOST_REQUIRE(l.read(result));
  BOOST_CHECK_EQUAL(result.first, 3);
  BOOST_CHECK_EQUAL(result.second, 30);
  BOOST_REQUIRE(l.read(result));
  BOOST_CHECK_EQUAL(result.first, 4);
  BOOST_CHECK_EQUAL(result.second, 40);

  BOOST_CHECK(!l.read(result));
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(emplace_noncopyable) {
  // This test is to make sure that emplace() with multiple parameters and
  // move-only operands compiles
  using ValueType = std::pair<std::unique_ptr<int>, std::unique_ptr<int>>;
  aocommon::Lane<ValueType> l(4);
  ValueType v(boost::make_unique<int>(5), boost::make_unique<int>(11));
  l.emplace(std::move(v));
  BOOST_CHECK(!l.empty());

  ValueType result;
  BOOST_REQUIRE(l.read(result));
  BOOST_REQUIRE(result.first);
  BOOST_REQUIRE(result.second);
  BOOST_CHECK_EQUAL(*result.first, 5);
  BOOST_CHECK_EQUAL(*result.second, 11);
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(noncopyable_read_write) {
  Lane<std::unique_ptr<int>> l(4);
  constexpr size_t n = 2;
  const int values[n] = {3, 4};
  for (size_t i = 0; i != n; ++i) {
    l.write(boost::make_unique<int>(values[i]));
  }
  l.write_end();
  BOOST_CHECK(!l.empty());
  BOOST_CHECK_EQUAL(l.size(), n);

  for (size_t i = 0; i != n; ++i) {
    std::unique_ptr<int> result;
    BOOST_REQUIRE(l.read(result));
    BOOST_REQUIRE(result);
    BOOST_CHECK_EQUAL(*result, values[i]);
  }
  std::unique_ptr<int> result;
  BOOST_CHECK(!l.read(result));
  BOOST_CHECK(l.empty());
  BOOST_CHECK(!result);
}

BOOST_AUTO_TEST_CASE(move_write) {
  Lane<std::unique_ptr<int>> l(4);
  std::vector<std::unique_ptr<int>> values;
  values.emplace_back();
  values.emplace_back(boost::make_unique<int>(7));
  l.move_write(values.data(), values.size());
  l.write_end();
  BOOST_CHECK(!l.empty());
  BOOST_CHECK_EQUAL(l.size(), 2);

  std::unique_ptr<int> result;
  BOOST_REQUIRE(l.read(result));
  BOOST_CHECK(!result);

  BOOST_REQUIRE(l.read(result));
  BOOST_REQUIRE(result);
  BOOST_CHECK_EQUAL(*result, 7);

  BOOST_CHECK(!l.read(result));
  BOOST_CHECK(l.empty());
  BOOST_REQUIRE(result);
  BOOST_CHECK_EQUAL(*result, 7);
}

BOOST_AUTO_TEST_CASE(discard) {
  Lane<int> l(4);
  constexpr size_t n = 3;
  const int values[n] = {3, 4, 5};
  l.write(values, n);
  l.write_end();
  l.discard(2);
  BOOST_CHECK(!l.empty());
  BOOST_CHECK_EQUAL(l.size(), 1);

  int result;
  BOOST_REQUIRE(l.read(result));
  BOOST_CHECK_EQUAL(result, values[2]);

  BOOST_CHECK(!l.read(result));
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(clear) {
  constexpr size_t n = 2;
  const int values_a[n] = {3, 4};
  Lane<int> l(3);
  l.write(values_a, n);
  l.write_end();
  l.clear();
  BOOST_CHECK_EQUAL(l.capacity(), 3);
  BOOST_CHECK(l.empty());
  BOOST_CHECK_EQUAL(l.size(), 0);
  const int values_b[n] = {5, 6};
  l.write(values_b, n);
  l.write_end();
  for (size_t i = 0; i != n; ++i) {
    int result;
    BOOST_REQUIRE(l.read(result));
    BOOST_CHECK_EQUAL(result, values_b[i]);
  }
  int result;
  BOOST_CHECK(!l.read(result));
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(resize) {
  Lane<int> l(1);
  l.write(17);
  l.write_end();
  l.resize(3);
  BOOST_CHECK_EQUAL(l.capacity(), 3);
  BOOST_CHECK(l.empty());
  BOOST_CHECK_EQUAL(l.size(), 0);
  constexpr size_t n = 2;
  const int values[n] = {3, 4};
  l.write(values, n);
  l.write_end();
  for (size_t i = 0; i != n; ++i) {
    int result;
    BOOST_REQUIRE(l.read(result));
    BOOST_CHECK_EQUAL(result, values[i]);
  }
  int result;
  BOOST_CHECK(!l.read(result));
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_CASE(wait_for_empty) {
  Lane<int> l(1);
  BOOST_CHECK_NO_THROW(l.wait_for_empty());
}

BOOST_AUTO_TEST_CASE(multiple_threads) {
  Lane<int> l(5);

  std::thread t([&l]() {
    for (size_t i = 100; i != 150; ++i) {
      int result;
      BOOST_REQUIRE(l.read(result));
      BOOST_CHECK_EQUAL(result, i);
    }
    l.discard(25);
    for (size_t i = 175; i != 250; ++i) {
      int result;
      BOOST_REQUIRE(l.read(result));
      BOOST_CHECK_EQUAL(result, i);
    }
  });

  for (size_t i = 100; i != 190; ++i) l.write(i);
  for (size_t i = 190; i != 240; ++i) l.emplace(i);
  std::vector<int> values;
  for (size_t i = 240; i != 250; ++i) values.emplace_back(i);
  l.write(values.data(), values.size());
  l.wait_for_empty();
  l.write_end();

  t.join();
  BOOST_CHECK(l.empty());
}

BOOST_AUTO_TEST_SUITE_END()
