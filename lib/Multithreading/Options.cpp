#ifdef THREADS_ENABLE

#include "marco/Runtime/Multithreading/Options.h"

namespace marco::runtime::multithreading {
MultithreadingOptions &multithreadingOptions() {
  static MultithreadingOptions obj;
  return obj;
}
} // namespace marco::runtime::multithreading

#endif // THREADS_ENABLE
