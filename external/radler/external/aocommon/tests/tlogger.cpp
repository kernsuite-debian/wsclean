#include <aocommon/logger.h>

#include <sstream>

#include <boost/test/unit_test.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <iostream>

BOOST_AUTO_TEST_SUITE(logger)

BOOST_AUTO_TEST_CASE(verbosity) {
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kNormal);
  BOOST_CHECK_EQUAL(aocommon::Logger::IsVerbose(), false);
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kVerbose);
  BOOST_CHECK_EQUAL(aocommon::Logger::IsVerbose(), true);
}

BOOST_AUTO_TEST_CASE(log_memory) {
  aocommon::Logger::SetLogTime(false);
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kNormal);
  aocommon::Logger::SetLogMemory(true);
  std::stringstream output;
  aocommon::Logger::LogWriter<aocommon::Logger::kInfoLevel> logwriter(output);
  logwriter << "Test line.\n";
  BOOST_CHECK_NE(output.str().find(" GB] Test line."), std::string::npos);

  output = std::stringstream();
  aocommon::Logger::SetLogMemory(false);
  logwriter << "Test line.\n";
  BOOST_CHECK_EQUAL(output.str(), "Test line.\n");
}

BOOST_AUTO_TEST_CASE(logwriter) {
  aocommon::Logger::SetLogTime(false);
  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kNormal);
  aocommon::Logger::SetLogMemory(false);

  std::stringstream output;
  aocommon::Logger::LogWriter<aocommon::Logger::kInfoLevel> logwriter(output);

  std::string str = "is a";
  logwriter << "T" << 'h' << "is" << ' ' << str;
  BOOST_CHECK_EQUAL(output.str(), "This is a");

  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kQuiet);
  logwriter << " quiet ";
  BOOST_CHECK_EQUAL(output.str(), "This is a");

  aocommon::Logger::SetVerbosity(aocommon::LogVerbosityLevel::kNormal);

  logwriter << " test.";
  BOOST_CHECK_EQUAL(output.str(), "This is a test.");

  const size_t my_size_t = 10;
  const double my_double = 9.999;
  logwriter << " Numerical output? " << my_size_t << "/" << my_double << ".";
  BOOST_CHECK_EQUAL(output.str(),
                    "This is a test. Numerical output? 10/9.999.");
}

BOOST_AUTO_TEST_SUITE_END()
