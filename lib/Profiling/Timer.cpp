#include "marco/Runtime/Profiling/Timer.h"

using namespace ::std::chrono;

namespace marco::runtime::profiling {
Timer::Timer() { reset(); }

void Timer::start() {
#ifdef THREADS_ENABLE
  std::lock_guard<std::mutex> lockGuard(mutex);
#endif

  if (running_++ == 0) {
    start_ = steady_clock::now();
  }
}

void Timer::stop() {
#ifdef THREADS_ENABLE
  std::lock_guard<std::mutex> lockGuard(mutex);
#endif

  if (--running_ == 0) {
    accumulatedTime_ += elapsed();
  }
}

void Timer::reset() {
#ifdef THREADS_ENABLE
  std::lock_guard<std::mutex> lockGuard(mutex);
#endif

  accumulatedTime_ = duration_values<nanoseconds>::zero();
}

nanoseconds Timer::elapsed() const {
  return std::chrono::duration_cast<std::chrono::nanoseconds>(
      steady_clock::now() - start_);
}

std::chrono::nanoseconds Timer::totalElapsed() const {
#ifdef THREADS_ENABLE
  std::lock_guard<std::mutex> lockGuard(mutex);
#endif

  if (running_) {
    return accumulatedTime_ + elapsed();
  }

  return accumulatedTime_;
}
} // namespace marco::runtime::profiling
