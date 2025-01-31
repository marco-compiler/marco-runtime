#ifndef MARCO_RUNTIME_MULTITHREADING_BARRIER_H
#define MARCO_RUNTIME_MULTITHREADING_BARRIER_H

#include <cassert>
#include <condition_variable>

namespace marco::runtime {
class Barrier {
public:
  Barrier(size_t numOfThreads);

  Barrier(const Barrier &barrier) = delete;

  Barrier(Barrier &&barrier) = delete;

  ~Barrier() noexcept;

  Barrier &operator=(const Barrier &barrier) = delete;

  Barrier &operator=(Barrier &&barrier) = delete;

  void wait();

  void reset();

private:
  std::mutex mutex;
  std::condition_variable cv;
  size_t counter;
  size_t numOfThreads;
};
} // namespace marco::runtime

#endif // MARCO_RUNTIME_MULTITHREADING_BARRIER_H
