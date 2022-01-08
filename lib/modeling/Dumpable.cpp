#include <marco/modeling/Dumpable.h>

namespace marco::modeling::internal
{
  Dumpable::~Dumpable() = default;

  void Dumpable::dump() const
  {
    dump(std::clog);
  }
}