#ifndef MARCO_RUNTIME_SIMULATION_OPTIONS_H
#define MARCO_RUNTIME_SIMULATION_OPTIONS_H

#include <cstdint>

namespace marco::runtime::simulation
{
  struct Options
  {
    bool debug = false;

    double startTime = 0;
    double endTime = 10;

    // The factor multiplying the threads count when computing the total number
    // of equations chunks.
    // In other words, it is the amount of chunks each thread would process in
    // a perfectly balanced scenario.
    int64_t equationsChunksFactor = 10;

    // The number of steps to be executed by the scheduler, first in sequential
    // and then multithreaded mode, to decide on the execution strategy of the
    // scheduler itself.
    int64_t schedulerCalibrationRuns = 10;
  };

  Options& getOptions();
}

#endif // MARCO_RUNTIME_SIMULATION_OPTIONS_H
