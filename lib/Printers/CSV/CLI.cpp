#ifdef CLI_ENABLE

#include "marco/Runtime/Printers/CSV/CLI.h"
#include "marco/Runtime/Printers/CSV/Options.h"

namespace marco::runtime::printing {
std::string CommandLineOptions::getTitle() const { return "Formatting"; }

void CommandLineOptions::printCommandLineOptions(std::ostream &os) const {
  // clang-format off
  os << "  --scientific-notation             Print the values using the scientific notation." << std::endl;
  os << "  --precision=<value>               Set the number of decimals to be printed. Defaults to " << printOptions().precision << "." << std::endl;
  os << "  --disable-printing                Disables output. Useful for testing performance without I/O overhead. Defaults to " << printOptions().disablePrinting << "." << std::endl;
  // clang-format on
}

void CommandLineOptions::parseCommandLineOptions(
    const argh::parser &options) const {
  // clang-format off
    printOptions().scientificNotation = options["scientific-notation"];
    options("precision") >> printOptions().precision;
    printOptions().disablePrinting = options["disable-printing"];
  // clang-format on
}
} // namespace marco::runtime::printing

#endif // CLI_ENABLE
