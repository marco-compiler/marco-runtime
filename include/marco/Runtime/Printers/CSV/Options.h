#ifndef MARCO_RUNTIME_PRINTING_CONFIG_H
#define MARCO_RUNTIME_PRINTING_CONFIG_H

#include <cstddef>

namespace marco::runtime::printing {
struct PrintOptions {
  bool scientificNotation = false;
  unsigned int precision = 9;

  // Buffer size in bytes.
  size_t bufferSize = 10 * 1024 * 1024;
};

PrintOptions &printOptions();
} // namespace marco::runtime::printing

#endif // MARCO_RUNTIME_PRINTING_CONFIG_H
