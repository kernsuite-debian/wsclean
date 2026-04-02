#ifndef AOCOMMON_LOGGER_H_
#define AOCOMMON_LOGGER_H_

#include <sys/time.h>
#include <sys/resource.h>

#include <iostream>
#include <mutex>
#include <sstream>

#include <boost/algorithm/string.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>

namespace aocommon {

enum class LogVerbosityLevel { kQuiet, kNormal, kVerbose };

inline LogVerbosityLevel StringToLogVerbosityLevel(
    std::string verbosity_string) {
  boost::algorithm::to_upper(verbosity_string);
  if (verbosity_string == "NORMAL")
    return LogVerbosityLevel::kNormal;
  else if (verbosity_string == "QUIET")
    return LogVerbosityLevel::kQuiet;
  else if (verbosity_string == "VERBOSE")
    return LogVerbosityLevel::kVerbose;
  else
    throw std::runtime_error(
        "Invalid log verbosity level: expected 'normal', 'quiet' or 'verbose'");
}

class Logger {
 public:
  enum LoggerLevel {
    kNoLevel = 5,
    kFatalLevel = 4,
    kErrorLevel = 3,
    kWarningLevel = 2,
    kInfoLevel = 1,
    kDebugLevel = 0
  };

  template <LoggerLevel Level>
  class LogWriter {
   public:
    /**
     * @brief Construct a new Log Writer object. Requires explicitly specifying
     * the output stream object, to enable testing. Internal to the @ref Logger
     * class, only std::cout will be used.
     */
    explicit LogWriter(std::ostream& output)
        : output_(output), at_new_line_(true) {}

    LogWriter& operator<<(const std::string& str) {
      std::lock_guard<std::mutex> lock(mutex_);
      size_t start = 0;
      size_t end;
      while (std::string::npos != (end = str.find('\n', start))) {
        OutputLinePart(str.substr(start, end - start + 1), true);
        start = end + 1;
      }
      OutputLinePart(str.substr(start, str.size() - start), false);
      return *this;
    }
    LogWriter& operator<<(const char* str) {
      (*this) << std::string(str);
      return *this;
    }
    LogWriter& operator<<(const char c) {
      std::lock_guard<std::mutex> lock(mutex_);
      OutputLinePart(std::string(1, c), c == '\n');
      return *this;
    }
    template <typename S>
    LogWriter& operator<<(const S& str) {
      std::ostringstream stream;
      stream << str;
      (*this) << stream.str();
      return *this;
    }
    void Flush() {
      std::lock_guard<std::mutex> lock(mutex_);
      output_.flush();
    }

   private:
    std::mutex mutex_;
    std::ostream& output_;
    bool at_new_line_;

    void OutputLinePart(const std::string& str, bool ends_with_cr) {
      if ((int)cout_level_ <= (int)Level && !str.empty()) {
        if (at_new_line_) {
          if (log_time_) OutputTime(output_);
          if (log_memory_) OutputMemory(output_);
        }
        output_ << str;
        at_new_line_ = ends_with_cr;
      }
    }
  };

  static void SetVerbosity(LogVerbosityLevel verbosity_level) {
    switch (verbosity_level) {
      case LogVerbosityLevel::kQuiet:
        cout_level_ = kErrorLevel;
        break;
      case LogVerbosityLevel::kNormal:
        cout_level_ = kInfoLevel;
        break;
      case LogVerbosityLevel::kVerbose:
        cout_level_ = kDebugLevel;
        break;
    }
  }

  static bool IsVerbose() { return cout_level_ == kDebugLevel; }

  static void SetLogTime(bool log_time) { log_time_ = log_time; }

  /**
   * Adds a number in front of every line logged that specifies the
   * amount of GBs in use by the current process.
   */
  static void SetLogMemory(bool log_memory) { log_memory_ = log_memory; }

  static bool LogTime() { return log_time_; }

  static bool LogMemory() { return log_memory_; }

  static inline LogWriter<kDebugLevel> Debug{std::cout};
  static inline LogWriter<kInfoLevel> Info{std::cout};
  static inline LogWriter<kWarningLevel> Warn{std::cout};
  static inline LogWriter<kErrorLevel> Error{std::cout};
  static inline LogWriter<kFatalLevel> Fatal{std::cout};

 private:
  Logger() = default;

  static void OutputTime(std::ostream& output) {
    boost::posix_time::ptime t(boost::posix_time::microsec_clock::local_time());
    const std::string str = boost::posix_time::to_simple_string(t);
    output << str << ' ';
  }

  static void OutputMemory(std::ostream& output) {
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    output << '[' << usage.ru_maxrss / 1000000ul << " GB] ";
  }

  static inline LoggerLevel cout_level_ = LoggerLevel::kInfoLevel;

  static inline bool log_time_ = false;
  static inline bool log_memory_ = false;
};

class LogReceiver {
 public:
  LogReceiver()
      : Fatal(this), Error(this), Warn(this), Info(this), Debug(this) {}
  virtual ~LogReceiver() = default;

  template <enum Logger::LoggerLevel Level>
  class LevelReceiver {
   public:
    LevelReceiver(LogReceiver* parent) : parent_(parent) {}
    LevelReceiver& operator<<(const std::string& str) {
      size_t start = 0, end;
      while (std::string::npos != (end = str.find('\n', start))) {
        parent_->Output(Level, str.substr(start, end - start + 1));
        start = end + 1;
      }
      parent_->Output(Level, str.substr(start, str.size() - start));
      return *this;
    }
    LevelReceiver& operator<<(const char* str) {
      (*this) << std::string(str);
      return *this;
    }
    LevelReceiver& operator<<(const char c) {
      parent_->Output(Level, std::string(1, c));
      return *this;
    }
    template <typename S>
    LevelReceiver& operator<<(const S& str) {
      std::ostringstream stream;
      stream << str;
      (*this) << stream.str();
      return *this;
    }

   private:
    LogReceiver* parent_;
  };  // end of class LevelReceiver

  LevelReceiver<Logger::kFatalLevel> Fatal;
  LevelReceiver<Logger::kErrorLevel> Error;
  LevelReceiver<Logger::kWarningLevel> Warn;
  LevelReceiver<Logger::kInfoLevel> Info;
  LevelReceiver<Logger::kDebugLevel> Debug;

 protected:
  virtual void Output(enum Logger::LoggerLevel level,
                      const std::string& str) = 0;

  void Forward(enum Logger::LoggerLevel level, const std::string& str) {
    switch (level) {
      case Logger::kFatalLevel:
        Logger::Fatal << str;
        break;
      case Logger::kErrorLevel:
        Logger::Error << str;
        break;
      case Logger::kWarningLevel:
        Logger::Warn << str;
        break;
      case Logger::kInfoLevel:
        Logger::Info << str;
        break;
      case Logger::kDebugLevel:
        Logger::Debug << str;
        break;
      case Logger::kNoLevel:
        break;
    }
  }
};

class ForwardingLogReceiver final : public LogReceiver {
 protected:
  void Output(enum Logger::LoggerLevel level, const std::string& str) override {
    Forward(level, str);
  }
};

}  // namespace aocommon
#endif
