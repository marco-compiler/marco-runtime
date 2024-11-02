#ifndef MARCO_RUNTIME_PRINTERS_DOUBLEBUFFER_H
#define MARCO_RUNTIME_PRINTERS_DOUBLEBUFFER_H

#include "marco/Runtime/Printers/Printer.h"
#include <functional>

namespace marco::runtime::printing {
class DoubleBuffer {
public:
  DoubleBuffer(uint64_t lineElementsCount, size_t bufferSizeBytes,
               std::function<void(const double *, uint64_t)> callback);

  DoubleBuffer(const DoubleBuffer &other) = delete;

  DoubleBuffer(DoubleBuffer &&other) = default;

  DoubleBuffer &operator=(const DoubleBuffer &other) = delete;

  DoubleBuffer &operator=(DoubleBuffer &&other) = default;

  double *getActiveBuffer();

  void endLine();

  void flush();

private:
  uint64_t lineElementsCount{1};
  uint64_t lines{1};
  std::function<void(const double *, uint64_t)> callback;

  std::vector<double> inputBuffer;
  std::vector<double> outputBuffer;

  uint64_t currentLine{0};
};
} // namespace marco::runtime::printing

#endif // MARCO_RUNTIME_PRINTERS_DOUBLEBUFFER_H
