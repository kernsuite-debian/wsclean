#include <aocommon/logger.h>

#include <sstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE(logger)

BOOST_AUTO_TEST_CASE(verbosity) {
  BOOST_CHECK_EQUAL(aocommon::Logger::IsVerbose(), false);
  aocommon::Logger::SetVerbosity(aocommon::Logger::kVerboseVerbosity);
  BOOST_CHECK_EQUAL(aocommon::Logger::IsVerbose(), true);
}

BOOST_AUTO_TEST_CASE(logwriter) {
  std::stringstream output;
  aocommon::Logger::LogWriter<aocommon::Logger::kInfoLevel> logwriter(output);

  std::string str = "is a";
  logwriter << "T"
            << "h"
            << "i"
            << "s"
            << " " << str;
  BOOST_CHECK_EQUAL(output.str(), "This is a");

  aocommon::Logger::SetVerbosity(aocommon::Logger::kQuietVerbosity);
  logwriter << " quiet ";
  BOOST_CHECK_EQUAL(output.str(), "This is a");

  aocommon::Logger::SetVerbosity(aocommon::Logger::kNormalVerbosity);

  logwriter << " test.";
  BOOST_CHECK_EQUAL(output.str(), "This is a test.");

  const size_t my_size_t = 10;
  const double my_double = 9.999;
  logwriter << " Numerical output? " << my_size_t << "/" << my_double << ".";
  BOOST_CHECK_EQUAL(output.str(),
                    "This is a test. Numerical output? 10/9.999.");
}

BOOST_AUTO_TEST_SUITE_END()
