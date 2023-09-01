#ifndef MARCO_RUNTIME_CLI_CATEGORY_H
#define MARCO_RUNTIME_CLI_CATEGORY_H

#ifdef CLI_ENABLE

#include "argh.h"
#include <string>

namespace marco::runtime::cli
{
  class Category
  {
    public:
      virtual ~Category();

      virtual std::string getTitle() const = 0;

      virtual void printCommandLineOptions(std::ostream& os) const = 0;

      virtual void parseCommandLineOptions(
          const argh::parser& options) const = 0;
  };
}

#endif // CLI_ENABLE

#endif // MARCO_RUNTIME_CLI_CATEGORY_H
