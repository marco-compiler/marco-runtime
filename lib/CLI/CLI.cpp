#ifdef CLI_ENABLE

#include "marco/Runtime/CLI/CLI.h"
#include <cassert>

namespace marco::runtime
{
  const cli::Category& CLI::operator[](size_t index) const
  {
    assert(index < categories.size());
    return *categories[index];
  }

  void CLI::operator+=(std::unique_ptr<cli::Category> category)
  {
    categories.push_back(std::move(category));
  }

  size_t CLI::size() const
  {
    return categories.size();
  }

  CLI& getCLI()
  {
    static CLI cli;
    return cli;
  }
}

#endif // CLI_ENABLE
