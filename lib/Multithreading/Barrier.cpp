#include "marco/Runtime/Multithreading/Barrier.h"

using namespace ::marco::runtime;

namespace marco::runtime {
Barrier::Barrier(size_t numOfThreads)
    : counter(numOfThreads), numOfThreads(numOfThreads) {}

Barrier::~Barrier() noexcept { assert(numOfThreads == 0u); }

void Barrier::wait() {
  std::unique_lock<std::mutex> lock(mutex);

  assert(counter != 0u);

  if (--counter == 0u) {
    cv.notify_all();
  } else {
    cv.wait(lock, [this]() { return counter == 0u; });
  }
}
} // namespace marco::runtime
