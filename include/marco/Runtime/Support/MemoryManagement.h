#ifndef MARCO_RUNTIME_MEMORYMANAGEMENT_H
#define MARCO_RUNTIME_MEMORYMANAGEMENT_H

#include "marco/Runtime/Support/Mangling.h"
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <vector>

extern "C" {
void *marco_malloc(int64_t size);
void *marco_realloc(void *ptr, int64_t size);
void marco_free(void *ptr);
};

namespace marco::runtime {
class MemoryPool {
public:
  MemoryPool() = default;

  MemoryPool(const MemoryPool &other) = delete;

  MemoryPool(MemoryPool &&other) = default;

  ~MemoryPool();

  MemoryPool &operator=(const MemoryPool &other) = delete;

  MemoryPool &operator=(MemoryPool &&other) = default;

  double *get(uint64_t id) const;

  uint64_t create(size_t numOfElements);

private:
  std::vector<double *> buffers;
};

class MemoryPoolManager {
private:
  MemoryPoolManager() = default;

public:
  static MemoryPoolManager &getInstance();

  MemoryPoolManager(const MemoryPoolManager &other) = delete;

  MemoryPoolManager(MemoryPoolManager &&other) = default;

  ~MemoryPoolManager() = default;

  MemoryPoolManager &operator=(const MemoryPoolManager &other) = delete;

  MemoryPoolManager &operator=(MemoryPoolManager &&other) = default;

  MemoryPool &get(uint64_t pool) const;

  uint64_t create();

private:
  std::vector<std::unique_ptr<MemoryPool>> pools{};
};
} // namespace marco::runtime

RUNTIME_FUNC_DECL(memoryPoolGet, PTR(void), int64_t, int64_t)

#endif // MARCO_RUNTIME_MEMORYMANAGEMENT_H
