#include <aocommon/barrier.h>

#include <thread>

#include <boost/test/unit_test.hpp>

using aocommon::Barrier;

namespace {
void nop() {}
}  // namespace

BOOST_AUTO_TEST_SUITE(barrier)

BOOST_AUTO_TEST_CASE(basic_use) {
  BOOST_CHECK_NO_THROW(Barrier b(1, nop); b.wait());

  BOOST_CHECK_NO_THROW(Barrier c(2, nop); std::thread t([&]() { c.wait(); });
                       c.wait(); t.join(););
}

BOOST_AUTO_TEST_SUITE_END()
