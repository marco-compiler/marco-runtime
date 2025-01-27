#include "marco/Runtime/Simulation/Runtime.h"
#include "marco/Runtime/CLI/CLI.h"
#include "marco/Runtime/Drivers/Driver.h"
#include "marco/Runtime/Drivers/KINSOL/CLI.h"
#include "marco/Runtime/Multithreading/CLI.h"
#include "marco/Runtime/Printers/Printer.h"
#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Simulation/CLI.h"
#include "marco/Runtime/Simulation/Options.h"
#include "marco/Runtime/Simulation/Profiler.h"
#include <cassert>
#include <iostream>

using namespace ::marco::runtime;

//===---------------------------------------------------------------------===//
// CLI
//===---------------------------------------------------------------------===//

static void printHelp() {
  std::cout << "Modelica simulation.\n";
  std::cout << "Model: " << getModelName() << "\n";
  std::cout << "Generated with MARCO compiler.\n\n";

#ifdef CLI_ENABLE
  std::cout << "OPTIONS:\n";
  std::cout << "  --help    Display the available options.\n\n";

  auto &cli = getCLI();

  for (size_t i = 0; i < cli.size(); ++i) {
    std::cout << cli[i].getTitle() << "\n";
    cli[i].printCommandLineOptions(std::cout);
    std::cout << "\n";
  }
#endif // CLI_ENABLE
}

//===---------------------------------------------------------------------===//
// Simulation
//===---------------------------------------------------------------------===//

namespace marco::runtime {
int64_t Simulation::getNumOfPrintableScalarVariables() const {
  int64_t result = 0;

  for (int64_t variable : variablesPrintOrder) {
    result += getNumOfPrintableScalarVariables(variable);
  }

  return result;
}

int64_t Simulation::getNumOfPrintableScalarVariables(int64_t variable) const {
  if (!printableVariables[variable]) {
    // The variable must not be printed.
    return 0;
  }

  int64_t rank = variablesRanks[variable];

  if (rank == 0) {
    return 1;
  }

  int64_t result = 0;

  for (const auto &range : variablesPrintableIndices[variable]) {
    size_t rangeSize = 1;

    for (int64_t dim = 0; dim < rank; ++dim) {
      rangeSize *= range[dim].end - range[dim].begin;
    }

    result += rangeSize;
  }

  return result;
}

Printer *Simulation::getPrinter() {
  assert(printer != nullptr);
  return printer;
}

const Printer *Simulation::getPrinter() const {
  assert(printer != nullptr);
  return printer;
}

void Simulation::setPrinter(Printer *newPrinter) {
  assert(newPrinter != nullptr);
  printer = newPrinter;
}
} // namespace marco::runtime

namespace {
Simulation runtimeInit() {
  Simulation result;

  // Number of array variables of the model (both state and algebraic ones).
  int64_t numOfVariables = getNumOfVariables();

  // Pre-fetch the names.
  result.variablesNames.resize(numOfVariables);

  for (int64_t var = 0; var < numOfVariables; ++var) {
    result.variablesNames[var] = getVariableName(var);
  }

  // Pre-fetch the ranks.
  result.variablesRanks.resize(numOfVariables);

  for (int64_t var = 0; var < numOfVariables; ++var) {
    result.variablesRanks[var] = getVariableRank(var);
  }

  // Pre-fetch whether each variable is printable.
  result.printableVariables.resize(numOfVariables);

  for (int64_t var = 0; var < numOfVariables; ++var) {
    result.printableVariables[var] = isPrintable(var);
  }

  // Pre-fetch the printable indices.
  result.variablesPrintableIndices.resize(numOfVariables);
  result.variablesPrintOrder.resize(numOfVariables);

  for (int64_t var = 0; var < numOfVariables; ++var) {
    int64_t numOfPrintableRanges = getVariableNumOfPrintableRanges(var);
    result.variablesPrintableIndices[var].resize(numOfPrintableRanges);

    for (int64_t range = 0; range < numOfPrintableRanges; ++range) {
      int64_t rank = result.variablesRanks[var];
      result.variablesPrintableIndices[var][range].resize(rank);

      for (int64_t dim = 0; dim < rank; ++dim) {
        result.variablesPrintableIndices[var][range][dim].begin =
            getVariablePrintableRangeBegin(var, range, dim);

        result.variablesPrintableIndices[var][range][dim].end =
            getVariablePrintableRangeEnd(var, range, dim);
      }
    }

    result.variablesPrintOrder[var] = var;
  }

  // Pre-fetch the derivatives map.
  result.derivativesMap.resize(numOfVariables);

  for (int64_t var = 0; var < numOfVariables; ++var) {
    result.derivativesMap[var] = -1;
  }

  for (int64_t var = 0; var < numOfVariables; ++var) {
    int64_t derivative = getDerivative(var);

    if (derivative != -1) {
      result.derivativesMap[derivative] = var;
    }
  }

  // Compute the derivative order of each variable.
  result.derOrders.resize(numOfVariables);

  for (int64_t var = 0; var < numOfVariables; ++var) {
    int64_t currentVar = result.derivativesMap[var];
    int64_t order = 0;

    while (currentVar != -1) {
      ++order;
      currentVar = result.derivativesMap[currentVar];
    }

    result.derOrders[var] = order;
  }

  // Determine the print ordering for the variables.
  std::sort(result.variablesPrintOrder.begin(),
            result.variablesPrintOrder.end(),
            [&](const int64_t &x, const int64_t &y) -> bool {
              int64_t xDerOrder = result.derOrders[x];
              int64_t yDerOrder = result.derOrders[y];

              if (xDerOrder < yDerOrder) {
                return true;
              }

              if (xDerOrder > yDerOrder) {
                return false;
              }

              int64_t first = x;
              int64_t second = y;

              if (xDerOrder != 0) {
                for (int64_t i = 0; i < xDerOrder; ++i) {
                  first = result.derivativesMap[first];
                }

                for (int64_t i = 0; i < yDerOrder; ++i) {
                  second = result.derivativesMap[second];
                }
              }

              return std::string_view(result.variablesNames[first]) <
                     std::string_view(result.variablesNames[second]);
            });

  return result;
}

void runtimeDeinit(Simulation &simulationInfo) {
#ifdef MARCO_PROFILING
  if (simulation::getOptions().profiling) {
    printProfilingStats();
  }
#endif
}
} // namespace

[[maybe_unused]] int runSimulation(int argc, char *argv[]) {
  // Initialize the runtime library.
  Simulation simulation = runtimeInit();

  // Instantiate the simulation driver.
  std::unique_ptr<Driver> driver = getDriver(&simulation);

  // Instantiate the data printer.
  std::unique_ptr<Printer> printer = getPrinter(&simulation);
  simulation.setPrinter(printer.get());

  // Parse the command-line arguments.
#ifdef CLI_ENABLE
  auto &cli = getCLI();
  cli += simulation::getCLIOptions();

#ifdef THREADS_ENABLE
  cli += multithreading::getCLIOptions();
#endif

  cli += driver->getCLIOptions();

#ifdef SHARED_DEPS
  cli += std::make_unique<sundials::kinsol::CommandLineOptions>();
#endif

  cli += printer->getCLIOptions();

  argh::parser cmdl(argc, argv);

  if (cmdl["help"]) {
    printHelp();
    return EXIT_SUCCESS;
  }

  for (size_t i = 0; i < cli.size(); ++i) {
    cli[i].parseCommandLineOptions(cmdl);
  }
#endif // CLI_ENABLE

  SIMULATION_PROFILER_INIT_START
  init();
  SIMULATION_PROFILER_INIT_STOP

  // Set the start time.
  setTime(simulation::getOptions().startTime);

  // Tell the printer that the simulation has begun.
  simulation.getPrinter()->simulationBegin();

  // Compute the initial conditions and print their values.
  SIMULATION_PROFILER_INITIAL_MODEL_START
  icModelBegin();
  solveICModel();
  icModelEnd();
  SIMULATION_PROFILER_INITIAL_MODEL_STOP

  simulation.getPrinter()->printValues();

  // Solve the dynamic model.
  SIMULATION_PROFILER_DYNAMIC_MODEL_START
  dynamicModelBegin();
  int result = driver->run();
  dynamicModelEnd();
  SIMULATION_PROFILER_DYNAMIC_MODEL_STOP

  // Tell the printer that the simulation has finished.
  simulation.getPrinter()->simulationEnd();

  deinit();

  // De-initialize the runtime library.
  runtimeDeinit(simulation);

  return result;
}
