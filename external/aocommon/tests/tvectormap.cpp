#include <aocommon/vectormap.h>

#include <thread>

#include <boost/test/unit_test.hpp>

using aocommon::VectorMap;

BOOST_AUTO_TEST_SUITE(vector_map)

BOOST_AUTO_TEST_CASE(empty) {
  const VectorMap<double> map;
  BOOST_CHECK(map.Empty());
  BOOST_CHECK_EQUAL(map.Size(), 0);
  BOOST_CHECK(map.begin() == map.end());
  BOOST_CHECK(map.Find(3) == map.end());
}

BOOST_AUTO_TEST_CASE(indexing) {
  VectorMap<char> map;
  BOOST_CHECK_EQUAL(map.AlwaysEmplace(1, 'c'), 'c');
  BOOST_CHECK(!map.Empty());
  BOOST_CHECK(map.begin() != map.end());
  BOOST_CHECK_EQUAL(*map.Find(1), 'c');
  BOOST_CHECK_EQUAL(map[1], 'c');
  BOOST_CHECK(map.Find(2) == map.end());

  map[0] = 'a';
  BOOST_CHECK_EQUAL(map[0], 'a');
  BOOST_CHECK_EQUAL(map[1], 'c');

  map[1] = 'b';
  BOOST_CHECK_EQUAL(map[0], 'a');
  BOOST_CHECK_EQUAL(map[1], 'b');

  map.EmplaceBack('c');

  BOOST_CHECK_EQUAL(*map.begin(), 'a');
  BOOST_CHECK_EQUAL(*(map.begin() + 1), 'b');
  BOOST_CHECK_EQUAL(*(map.begin() + 2), 'c');
  BOOST_CHECK(map.begin() + 3 == map.end());

  map.Clear();
  BOOST_CHECK(map.Empty());
}

BOOST_AUTO_TEST_CASE(conditional_methods) {
  const VectorMap<bool> empty_map;
  BOOST_CHECK_EQUAL(empty_map.GetKeysIf([](bool) { return true; }).size(), 0);
  BOOST_CHECK_EQUAL(empty_map.CountKeysIf([](bool) { return true; }), 0);

  const VectorMap<int> map{-3, -1, 0, 4, 5};
  BOOST_CHECK_EQUAL(map.GetKeysIf([](int) { return false; }).size(), 0);
  BOOST_CHECK_EQUAL(map.CountKeysIf([](int) { return false; }), 0);

  const std::set<size_t> keys_a = map.GetKeysIf([](int a) { return a >= 0; });
  const std::set<size_t> expected_a = {2, 3, 4};
  BOOST_CHECK(keys_a == expected_a);
  BOOST_CHECK_EQUAL(map.CountKeysIf([](int a) { return a >= 0; }), 3);

  const std::set<size_t> keys_b = map.GetKeysIf([](int a) { return true; });
  const std::set<size_t> expected_b = {0, 1, 2, 3, 4};
  BOOST_CHECK(keys_b == expected_b);
  BOOST_CHECK_EQUAL(map.CountKeysIf([](int a) { return true; }), 5);
}

BOOST_AUTO_TEST_SUITE_END()
