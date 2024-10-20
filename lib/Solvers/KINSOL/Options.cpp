#include "marco/Runtime/Solvers/KINSOL/Options.h"

#ifdef SUNDIALS_ENABLE

namespace marco::runtime::sundials::kinsol {
Options &getOptions() {
  static Options options;
  return options;
}
} // namespace marco::runtime::sundials::kinsol

#endif // SUNDIALS_ENABLE
