#ifdef CLI_ENABLE

#include "marco/Runtime/Simulation/CLI.h"
#include "marco/Runtime/Simulation/Options.h"
#include <iostream>

namespace marco::runtime::simulation {
std::string CommandLineOptions::getTitle() const { return "General"; }

void CommandLineOptions::printCommandLineOptions(std::ostream &os) const {
  // clang-format off
  os << "  --debug                          Enable the debug messages." << std::endl;
  os << "  --profile                        Enable the profilers." << std::endl;

  os << "  --start-time=<value>             Set the start time (in seconds). Defaults to " << getOptions().startTime << "." << std::endl;
  os << "  --end-time=<value>               Set the end time (in seconds). Defaults to " << getOptions().endTime << "." << std::endl;
  os << "  --equations-partitioning-factor  Set the amount of equation partitions each thread would process in an ideal scenario where all the equations are independent from each other and have equal computational cost. Defaults to " << getOptions().equationsPartitioningFactor << "." << std::endl;
  os << "  --scheduler-calibration-runs     Set the amount of sequential and multithreaded executions used to decide the execution policy. Defaults to " << getOptions().schedulerCalibrationRuns << "." << std::endl;
  os << "  --scheduler-policy=<value>       Force the schedulers to adopt a certain execution policy (sequential / multithreaded)." << std::endl;
  // clang-format on
}

void CommandLineOptions::parseCommandLineOptions(
    const argh::parser &options) const {
  // clang-format off
  getOptions().debug = options["debug"];
  getOptions().profiling = options["profile"];

  options("start-time") >> getOptions().startTime;
  options("end-time") >> getOptions().endTime;
  options("equations-partitioning-factor") >> getOptions().equationsPartitioningFactor;
  options("scheduler-calibration-runs") >> getOptions().schedulerCalibrationRuns;

  std::string schedulerPolicy;
  options("scheduler-policy") >> schedulerPolicy;

  if (schedulerPolicy == "sequential") {
    getOptions().schedulerPolicy = SchedulerPolicy::Sequential;
  } else if (schedulerPolicy == "multithreaded") {
    getOptions().schedulerPolicy = SchedulerPolicy::Multithreaded;
  }

  // clang-format on
}

std::unique_ptr<cli::Category> getCLIOptions() {
  return std::make_unique<CommandLineOptions>();
}
} // namespace marco::runtime::simulation

#endif // CLI_ENABLE
