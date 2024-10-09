#include "marco/Runtime/Support/MemoryManagement.h"
#include <cassert>
#include <cstdlib>
#include <memory>
#include <vector>

using namespace ::marco::runtime;

#ifdef MARCO_PROFILING

#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Profiling/Timer.h"
#include <iostream>
#include <map>
#include <mutex>

namespace marco::runtime::profiling
{
  class MemoryProfiler : public Profiler
  {
    public:
      MemoryProfiler() : Profiler("Memory management")
      {
        registerProfiler(*this);
      }

      void reset() override
      {
        std::lock_guard<std::mutex> lockGuard(mutex);

        mallocCalls = 0;
        reallocCalls = 0;
        freeCalls = 0;
        totalHeapMemory = 0;
        currentHeapMemory = 0;
        peakHeapMemory = 0;
        timer.reset();
      }

      void print() const override
      {
        std::lock_guard<std::mutex> lockGuard(mutex);

        std::cerr << "Number of 'malloc' invocations: " << mallocCalls << "\n";
        std::cerr << "Number of 'realloc' invocations: " << reallocCalls << "\n";
        std::cerr << "Number of 'free' invocations: " << freeCalls << "\n";

        if (mallocCalls > reallocCalls + freeCalls) {
          std::cerr << "[Warning] Possible memory leak detected\n";
        } else if (mallocCalls + reallocCalls < freeCalls) {
          std::cerr << "[Warning] Possible double 'free' detected\n";
        }

        std::cerr << "Total amount of heap allocated memory: " << totalHeapMemory << " bytes\n";
        std::cerr << "Peak of heap memory usage: " << peakHeapMemory << " bytes\n";
        std::cerr << "Time spent on heap memory management: " << time() << " ms\n";
      }

      void malloc(void* address, int64_t bytes)
      {
        std::lock_guard<std::mutex> lockGuard(mutex);

        ++mallocCalls;

        totalHeapMemory += bytes;
        currentHeapMemory += bytes;
        sizes[address] = bytes;

        if (currentHeapMemory > peakHeapMemory) {
          peakHeapMemory = currentHeapMemory;
        }
      }

      void realloc(void* previous, void* current, int64_t bytes)
      {
        std::lock_guard<std::mutex> lockGuard(mutex);

        ++reallocCalls;

        totalHeapMemory -= sizes[previous];
        currentHeapMemory -= sizes[previous];

        totalHeapMemory += bytes;
        currentHeapMemory += bytes;
        sizes[current] = bytes;

        if (currentHeapMemory > peakHeapMemory) {
          peakHeapMemory = currentHeapMemory;
        }
      }

      void free(void* address)
      {
        std::lock_guard<std::mutex> lockGuard(mutex);

        ++freeCalls;

        if (auto it = sizes.find(address); it != sizes.end()) {
          currentHeapMemory -= it->second;
          sizes.erase(it);
        }
      }

      void startTimer()
      {
        std::lock_guard<std::mutex> lockGuard(mutex);
        timer.start();
      }

      void stopTimer()
      {
        std::lock_guard<std::mutex> lockGuard(mutex);
        timer.stop();
      }

    private:
      double time() const
      {
        return timer.totalElapsedTime<std::milli>();
      }

    private:
      size_t mallocCalls;
      size_t reallocCalls;
      size_t freeCalls;
      int64_t totalHeapMemory;
      int64_t currentHeapMemory;
      int64_t peakHeapMemory;
      std::map<void*, int64_t> sizes;
      Timer timer;

      mutable std::mutex mutex;
  };
}

namespace
{
  marco::runtime::profiling::MemoryProfiler& profiler()
  {
    static marco::runtime::profiling::MemoryProfiler obj;
    return obj;
  }
}

#endif

void* marco_malloc(int64_t sizeInBytes)
{
#ifdef MARCO_PROFILING
  ::profiler().startTimer();
#endif

  void* result = sizeInBytes == 0 ? nullptr : std::malloc(sizeInBytes);

#ifdef MARCO_PROFILING
  ::profiler().stopTimer();
  ::profiler().malloc(result, sizeInBytes);
#endif

  return result;
}

void* marco_realloc(void* ptr, int64_t sizeInBytes)
{
#ifdef MARCO_PROFILING
  ::profiler().startTimer();
#endif

  void* result = sizeInBytes == 0 ? nullptr : std::realloc(ptr, sizeInBytes);

#ifdef MARCO_PROFILING
  ::profiler().stopTimer();
  ::profiler().realloc(ptr, result, sizeInBytes);
#endif

  return result;
}

void marco_free(void* ptr)
{
#ifdef MARCO_PROFILING
  ::profiler().free(ptr);
  ::profiler().startTimer();
#endif

  if (ptr != nullptr) {
    std::free(ptr);
  }

#ifdef MARCO_PROFILING
  ::profiler().stopTimer();
#endif
}
#include <iostream>
namespace marco::runtime
{

MemoryPool::~MemoryPool() {
  for (double* buffer : buffers) {
    if (buffer != nullptr) {
      std::free(buffer);
    }
  }
}

double* MemoryPool::get(uint64_t id) const {
  assert(id < buffers.size());
  return buffers[id];
}

uint64_t MemoryPool::create(size_t numOfElements) {
  uint64_t id = buffers.size();
  buffers.push_back(static_cast<double*>(std::malloc(sizeof(double) * numOfElements)));
  return id;
}

MemoryPoolManager& MemoryPoolManager::getInstance() {
  static MemoryPoolManager instance;
  return instance;
}

MemoryPool& MemoryPoolManager::get(uint64_t pool) const {
  assert(pool < pools.size());
  return *pools[pool];
}

uint64_t MemoryPoolManager::create() {
  uint64_t id = pools.size();
  pools.push_back(std::make_unique<MemoryPool>());
  return id;
}
}

namespace {
void* memoryPoolGet_pvoid(uint64_t pool, uint64_t buffer) {
  MemoryPoolManager& manager = MemoryPoolManager::getInstance();
  return static_cast<void*>(manager.get(pool).get(buffer));
}
}

RUNTIME_FUNC_DEF(memoryPoolGet, PTR(void), int64_t, int64_t)
