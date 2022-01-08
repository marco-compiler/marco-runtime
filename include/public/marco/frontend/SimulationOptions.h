#ifndef MARCO_FRONTEND_SIMULATIONOPTIONS_H
#define MARCO_FRONTEND_SIMULATIONOPTIONS_H

namespace marco::frontend
{
  struct SimulationOptions
  {
    double startTime;
    double endTime;
    double timeStep;

    SimulationOptions();
  };
}

#endif // MARCO_FRONTEND_SIMULATIONOPTIONS_H
