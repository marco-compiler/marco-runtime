#include "marco/Runtime/Solvers/IDA/Options.h"

#ifdef SUNDIALS_ENABLE

namespace marco::runtime::sundials::ida {
Options &getOptions() {
  static Options options;
  return options;
}
} // namespace marco::runtime::sundials::ida

#endif // SUNDIALS_ENABLE
