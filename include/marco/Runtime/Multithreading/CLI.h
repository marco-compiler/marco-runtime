#ifndef MARCO_RUNTIME_MULTITHREADING_CLI_H
#define MARCO_RUNTIME_MULTITHREADING_CLI_H

#ifdef CLI_ENABLE
#ifdef THREADS_ENABLE

#include "marco/Runtime/CLI/CLI.h"

namespace marco::runtime::multithreading
{
  class CommandLineOptions : public cli::Category
  {
    public:
      std::string getTitle() const override;

      void printCommandLineOptions(std::ostream& os) const override;

      void parseCommandLineOptions(const argh::parser& options) const override;
  };

  std::unique_ptr<cli::Category> getCLIOptions();
}

#endif // THREADS_ENABLE
#endif // CLI_ENABLE

#endif // MARCO_RUNTIME_MULTITHREADING_CLI_H
