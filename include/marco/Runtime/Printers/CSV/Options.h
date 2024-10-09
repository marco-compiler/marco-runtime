#ifndef MARCO_RUNTIME_PRINTING_CONFIG_H
#define MARCO_RUNTIME_PRINTING_CONFIG_H

namespace marco::runtime::printing
{
  struct PrintOptions
  {
    bool scientificNotation = false;
    unsigned int precision = 9;

    // The frequency at which the simulation state should be printed.
    unsigned int snapshotFrequency = 1;
  };

  PrintOptions& printOptions();
}

#endif // MARCO_RUNTIME_PRINTING_CONFIG_H
