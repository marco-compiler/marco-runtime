#include "marco/Runtime/Printers/DoubleBuffer.h"

using namespace ::marco::runtime::printing;

namespace marco::runtime::printing {
DoubleBuffer::DoubleBuffer(
    uint64_t lineElementsCount, size_t bufferSizeBytes,
    std::function<void(const double *, uint64_t)> callback)
    : lineElementsCount(lineElementsCount), callback(callback) {
  uint64_t requestedLines =
      bufferSizeBytes / (sizeof(double) * lineElementsCount);

  lines = std::max(static_cast<uint64_t>(1), requestedLines);
  inputBuffer.resize(lines * lineElementsCount, 0xDEADBEEF);
  outputBuffer.resize(lines * lineElementsCount, 0xDEADBEEF);
}

double *DoubleBuffer::getActiveBuffer() {
  return inputBuffer.data() + currentLine * lineElementsCount;
}

void DoubleBuffer::endLine() {
  if (++currentLine >= lines) {
    flush();
  }
}

void DoubleBuffer::flush() {
  if (currentLine > 0) {
    std::swap(inputBuffer, outputBuffer);

    for (uint64_t line = 0; line < lines; ++line) {
      callback(outputBuffer.data() + line * lineElementsCount,
               lineElementsCount);
    }

    currentLine = 0;
  }
}
} // namespace marco::runtime::printing
