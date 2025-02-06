#include "marco/Runtime/Simulation/Scheduler.h"
#include "marco/Runtime/Multithreading/Barrier.h"
#include "marco/Runtime/Multithreading/ThreadPool.h"
#include "marco/Runtime/Simulation/Options.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <iostream>
#include <optional>

using namespace ::marco::runtime;

//===---------------------------------------------------------------------===//
// Profiling
//===---------------------------------------------------------------------===//

#ifdef MARCO_PROFILING
#include "marco/Runtime/Simulation/Profiler.h"

SchedulerProfiler::SchedulerProfiler(int64_t schedulerId)
    : Profiler("Scheduler " + std::to_string(schedulerId)) {}

void SchedulerProfiler::createPartitionsGroupsCounters(size_t amount) {
  partitionsGroupsCounters.clear();
  partitionsGroupsCounters.resize(amount, 0);
}

void SchedulerProfiler::createPartitionsGroupsTimers(size_t amount) {
  partitionsGroups.clear();

  for (size_t i = 0; i < amount; ++i) {
    partitionsGroups.push_back(std::make_unique<profiling::Timer>());
  }
}

void SchedulerProfiler::reset() {
  std::lock_guard<std::mutex> lockGuard(mutex);

  addEquation.reset();
  initialization.reset();
  run.reset();
  sequentialRuns = 0;
  multithreadedRuns = 0;

  for (auto &partitionsGroup : partitionsGroups) {
    partitionsGroup->reset();
  }
}

void SchedulerProfiler::print() const {
  std::lock_guard<std::mutex> lockGuard(mutex);

  std::cerr << "Time spent on adding the equations: "
            << addEquation.totalElapsedTime<std::milli>() << " ms" << std::endl;

  std::cerr << "Time spent on initialization: "
            << initialization.totalElapsedTime<std::milli>() << " ms"
            << std::endl;

  std::cerr << "Time spent on 'run' method: "
            << run.totalElapsedTime<std::milli>() << " ms" << std::endl;

  std::cerr << "Number of sequential executions: " << sequentialRuns
            << std::endl;

  std::cerr << "Number of multithreaded executions: " << multithreadedRuns
            << std::endl;

  for (size_t i = 0, e = partitionsGroups.size(); i < e; ++i) {
    auto partitionsGroupsCounter = partitionsGroupsCounters[i];

    double averagePartitionsGroupTime =
        partitionsGroups[i]->totalElapsedTime<std::nano>() /
        static_cast<double>(partitionsGroupsCounter);

    std::cerr << "\n";

    std::cerr << "Time (total) spent by thread #" << i
              << " in processing equations: "
              << partitionsGroups[i]->totalElapsedTime<std::milli>() << " ms"
              << std::endl;

    std::cerr << "Time (average) spent by thread #" << i
              << " in processing equations: " << averagePartitionsGroupTime
              << " ns" << std::endl;

    std::cerr << "Number of partitions groups processed: "
              << partitionsGroupsCounter << std::endl;
  }
}

// clang-format off
#define SCHEDULER_PROFILER_ADD_EQUATION_START                                 \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->addEquation.start();                                            \
  }

#define SCHEDULER_PROFILER_ADD_EQUATION_STOP                                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->addEquation.stop();                                             \
  }

#define SCHEDULER_PROFILER_RUN_START                                          \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->run.start();                                                    \
  }

#define SCHEDULER_PROFILER_RUN_STOP                                           \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->run.stop();                                                     \
  }

#define SCHEDULER_PROFILER_INCREMENT_SEQUENTIAL_RUNS_COUNTER                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ++profiler->sequentialRuns;                                               \
  }

#define SCHEDULER_PROFILER_INCREMENT_MULTITHREADED_RUNS_COUNTER               \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    ++profiler->multithreadedRuns;                                            \
  }

#define SCHEDULER_PROFILER_INITIALIZATION_START                               \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->initialization.start();                                         \
  }

#define SCHEDULER_PROFILER_INITIALIZATION_STOP                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->initialization.stop();                                          \
  }

#define SCHEDULER_PROFILER_PARTITIONS_GROUP_START(thread)                     \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->partitionsGroups[thread]->start();                              \
    profiler->partitionsGroupsCounters[thread]++;                             \
  }

#define SCHEDULER_PROFILER_PARTITIONS_GROUP_STOP(thread)                      \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    profiler->partitionsGroups[thread]->stop();                               \
  }
// clang-format on

#else

#define SCHEDULER_PROFILER_DO_NOTHING static_assert(true);

#define SCHEDULER_PROFILER_ADD_EQUATION_START SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_ADD_EQUATION_STOP SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_RUN_START SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_RUN_STOP SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_INCREMENT_SEQUENTIAL_RUNS_COUNTER                   \
  SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_INCREMENT_MULTITHREADED_RUNS_COUNTER                \
  SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_INITIALIZATION_START SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_INITIALIZATION_STOP SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_PARTITIONS_GROUP_START(thread)                      \
  SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_PARTITIONS_GROUP_STOP(thread)                       \
  SCHEDULER_PROFILER_DO_NOTHING

#endif

//===---------------------------------------------------------------------===//
// EquationPartition
//===---------------------------------------------------------------------===//

namespace marco::runtime {
void Scheduler::EquationPartition::run() {
  assert(equation != nullptr);
  equation->function(ranges.data());
}
} // namespace marco::runtime

//===---------------------------------------------------------------------===//
// BackwardEquation
//===---------------------------------------------------------------------===//

namespace marco::runtime {

Scheduler::BackwardEquation::BackwardEquation(
    const Scheduler::Equation &equation)
    : equation(&equation) {}

const Scheduler::Equation &Scheduler::BackwardEquation::getEquation() const {
  assert(equation != nullptr);
  return *equation;
}

void Scheduler::BackwardEquation::run(
    unsigned int threadId,
    std::function<std::optional<int64_t>(const Equation &, unsigned int)>
        claimRowFn,
    std::vector<char> *currentRowState, std::vector<char> *previousRowState) {
  std::optional<int64_t> row = claimRowFn(*equation, threadId);

  while (row) {
    MultidimensionalRange rowRange = getRowBoundaries(*row);
    std::vector<int64_t> indices = getBeginIndices(*row);
    std::vector<int64_t> functionRanges;
    functionRanges.resize(indices.size() * 2, 0xDEADF00D);

    do {
      // Compute the iteration boundaries that the equation function will use.
      for (size_t dim = 0, rank = indices.size(); dim < rank; ++dim) {
        functionRanges[dim * 2] = indices[dim];
        functionRanges[dim * 2 + 1] = indices[dim] + 1;
      }

      int64_t column = indices[1];

      if (previousRowState) {
        // Enforce volatile semantics so to avoid any compiler optimization on
        // the readiness flag.
        while (static_cast<volatile char>((*previousRowState)[column]) != 1) {
          // Wait until the dependencies are satisfied.
          // Synchronization is purposely not performed.
          // Using chars guarantees consistency.
        }
      }

      // Call the equation function on the specific grid cell.
      equation->function(functionRanges.data());

      // Mark the cell as processed, so that other threads can continue.
      static_cast<volatile char&>((*currentRowState)[column]) = 1;
    } while (advanceEquationIndices(indices, rowRange));

    // Move to the next row.
    row = claimRowFn(*equation, threadId);
  }
}

MultidimensionalRange
Scheduler::BackwardEquation::getRowBoundaries(int64_t row) const {
  std::vector<Range> ranges;

  ranges.emplace_back(equation->indices[0].begin + row,
                      equation->indices[0].begin + row + 1);
  ranges.push_back(equation->indices[1]);

  return ranges;
}

std::vector<int64_t>
Scheduler::BackwardEquation::getBeginIndices(int64_t row) const {
  std::vector<int64_t> ranges;

  ranges.push_back(equation->indices[0].begin + row);
  ranges.push_back(equation->indices[1].begin);

  return ranges;
}
} // namespace marco::runtime

//===---------------------------------------------------------------------===//
// Scheduler
//===---------------------------------------------------------------------===//

static int64_t getUniqueSchedulerIdentifier() {
  static int64_t identifier = 0;
  return identifier++;
}

// The thread pool is shared by all the schedulers.
// Having multiple ones would waste resources in instantiating new thread
// groups which would anyway be used one at a time.
static ThreadPool &getSchedulersThreadPool() {
  static ThreadPool instance;
  return instance;
}

static uint64_t getEquationPartitionFlatSize(
    const Scheduler::ContiguousEquationPartition &partition) {
  int64_t result = 1;

  assert(partition.second.size() % 2 == 0);
  size_t rank = partition.second.size() / 2;

  for (size_t dim = 0; dim < rank; ++dim) {
    auto lowerBound = partition.second[dim * 2];
    auto upperBound = partition.second[dim * 2 + 1];
    auto size = upperBound - lowerBound;
    result *= size;
  }

  return result;
}

namespace marco::runtime {
Scheduler::Equation::Equation(EquationFunction function,
                              MultidimensionalRange indices,
                              DependencyKind dependencyKind)
    : function(function), indices(std::move(indices)),
      dependencyKind(dependencyKind) {}

Scheduler::Scheduler() {
  identifier = getUniqueSchedulerIdentifier();

#ifdef MARCO_PROFILING
  if (simulation::getOptions().profiling) {
    ThreadPool &threadPool = getSchedulersThreadPool();
    unsigned int numOfThreads = threadPool.getNumOfThreads();

    profiler = std::make_shared<SchedulerProfiler>(identifier);
    profiler->createPartitionsGroupsCounters(numOfThreads);
    profiler->createPartitionsGroupsTimers(numOfThreads);

    registerProfiler(profiler);
  }
#endif
}

void Scheduler::addEquation(EquationFunction function, uint64_t rank,
                            int64_t *ranges, DependencyKind dependencyKind) {
  SCHEDULER_PROFILER_ADD_EQUATION_START;
  std::vector<Range> indices;

  for (uint64_t dim = 0; dim < rank; ++dim) {
    indices.emplace_back(ranges[dim * 2], ranges[dim * 2 + 1]);
  }

  if (simulation::getOptions().debug) {
    std::cerr << "[Scheduler " << identifier << "] New equation added"
              << std::endl;

    std::cerr << "  - Rank: " << rank << std::endl;

    if (rank != 0) {
      std::cerr << "  - Ranges: ";

      for (uint64_t i = 0; i < rank; ++i) {
        std::cerr << "[" << indices[i].begin << ", " << indices[i].end << ")";
      }

      std::cerr << std::endl;
    }
  }

  equations.emplace_back(function, std::move(indices), dependencyKind);
  SCHEDULER_PROFILER_ADD_EQUATION_STOP;
}

void Scheduler::initialize() {
  SCHEDULER_PROFILER_INITIALIZATION_START;
  assert(!initialized && "Scheduler already initialized");

  if (auto forcedPolicy = simulation::getOptions().schedulerPolicy;
      !forcedPolicy || *forcedPolicy == SchedulerPolicy::Sequential) {
    // Compute the sequential schedule.
    for (const Equation &equation : equations) {
      std::vector<int64_t> functionArgs;

      for (const auto &range : equation.indices) {
        functionArgs.push_back(range.begin);
        functionArgs.push_back(range.end);
      }

      sequentialSchedule.emplace_back(equation, functionArgs);
    }

    assert(std::all_of(equations.begin(), equations.end(),
                       [&](const Equation &equation) {
                         return checkEquationScheduledExactlyOnce(
                             equation, {sequentialSchedule});
                       }) &&
           "Not all the equations are scheduled exactly once in the sequential "
           "schedule");

    assert(std::all_of(sequentialSchedule.begin(), sequentialSchedule.end(),
                       [&](const ContiguousEquationPartition &partition) {
                         return checkEquationIndicesExistence(partition);
                       }) &&
           "Some nonexistent equation indices have been scheduled in the "
           "sequential schedule");
  }

  if (auto forcedPolicy = simulation::getOptions().schedulerPolicy;
      !forcedPolicy || *forcedPolicy == SchedulerPolicy::Multithreaded) {
    // Compute the multithreaded schedule.
    ThreadPool &threadPool = getSchedulersThreadPool();
    unsigned int numOfThreads = threadPool.getNumOfThreads();

    int64_t partitionsFactor =
        simulation::getOptions().equationsPartitioningFactor;

    int64_t numOfPartitions = numOfThreads * partitionsFactor;

    uint64_t numOfScalarEquations = 0;

    for (const Equation &equation : equations) {
      numOfScalarEquations += getFlatSize(equation.indices);
    }

    size_t partitionsGroupMaxFlatSize =
        (numOfScalarEquations + numOfPartitions - 1) / numOfPartitions;

    if (simulation::getOptions().debug) {
      std::cerr << "[Scheduler " << identifier << "] Initializing" << std::endl
                << "  - Number of equations: " << numOfScalarEquations
                << std::endl
                << "  - Number of threads: " << numOfThreads << std::endl
                << "  - Partitioning factor: " << partitionsFactor << std::endl
                << "  - Number of partitions: " << numOfPartitions << std::endl
                << "  - Max flat size of each partitions group: "
                << partitionsGroupMaxFlatSize << std::endl;
    }

    ContiguousEquationsGroup partitionsGroup;
    size_t partitionsGroupFlatSize = 0;

    auto pushPartitionsGroupFn = [&]() {
      if (simulation::getOptions().debug) {
        std::cerr << "[Scheduler " << identifier
                  << "] Adding equation partitions group" << std::endl;

        std::cerr << "  - Number of partitions: " << partitionsGroup.size()
                  << std::endl;

        uint64_t totalSize = 0;
        std::cerr << "  - Equation partition flat sizes: [";

        for (size_t i = 0, e = partitionsGroup.size(); i < e; ++i) {
          if (i != 0) {
            std::cerr << ", ";
          }

          uint64_t partitionFlatSize =
              getEquationPartitionFlatSize(partitionsGroup[i]);

          std::cerr << partitionFlatSize;
          totalSize += partitionFlatSize;
        }

        std::cerr << "]" << std::endl;

        for (const auto &partition : partitionsGroup) {
          std::cerr << "    - Function: "
                    << reinterpret_cast<void *>(partition.first.function)
                    << std::endl;

          std::cerr << "    - Range: ";

          assert(partition.second.size() % 2 == 0);
          size_t rank = partition.second.size() / 2;

          for (size_t dim = 0; dim < rank; ++dim) {
            auto lowerBound = partition.second[dim * 2];
            auto upperBound = partition.second[dim * 2 + 1];
            std::cerr << "[" << lowerBound << ", " << upperBound << ")";
          }

          std::cerr << std::endl;
        }

        std::cerr << "  - Total size: " << totalSize << std::endl;
      }

      multithreadedSchedule.contiguousEquations.push_back(
          std::move(partitionsGroup));
      partitionsGroup.clear();
      partitionsGroupFlatSize = 0;
    };

    for (const Equation &equation : equations) {
      uint64_t flatSize = getFlatSize(equation.indices);

      size_t remainingSpace =
          partitionsGroupMaxFlatSize - partitionsGroupFlatSize;

      if (simulation::getOptions().debug) {
        std::cerr << "[Scheduler " << identifier << "] Partitioning equation"
                  << std::endl;

        std::cerr << "  - Function: "
                  << reinterpret_cast<void *>(equation.function) << std::endl;

        std::cerr << "  - Ranges: ";

        for (const auto &range : equation.indices) {
          std::cerr << "[" << range.begin << ", " << range.end << ")";
        }

        std::cerr << std::endl;
        std::cerr << "  - Flat size: " << flatSize << std::endl;
        std::cerr << "   - Dependency kind: ";

        if (equation.dependencyKind == DependencyKind::Sequential) {
          std::cerr << "sequential";
        } else if (equation.dependencyKind == DependencyKind::Backward) {
          std::cerr << "backward";
        } else if (equation.dependencyKind ==
                   DependencyKind::IndependentIndices) {
          std::cerr << "independent indices";
        } else {
          std::cerr << "unknown";
        }

        std::cerr << std::endl;
        std::cerr << "  - Remaining space: " << remainingSpace << std::endl;
      }

      if (equation.dependencyKind == DependencyKind::IndependentIndices) {
        uint64_t equationFlatIndex = 0;
        size_t equationRank = equation.indices.size();

        // Divide the ranges.
        while (equationFlatIndex < flatSize) {
          uint64_t beginFlatIndex = equationFlatIndex;

          uint64_t endFlatIndex = std::min(
              beginFlatIndex + static_cast<uint64_t>(remainingSpace), flatSize);

          assert(endFlatIndex > 0);
          --endFlatIndex;

          std::vector<int64_t> beginIndices;
          std::vector<int64_t> endIndices;

          getIndicesFromFlatIndex(beginFlatIndex, beginIndices,
                                  equation.indices);

          getIndicesFromFlatIndex(endFlatIndex, endIndices, equation.indices);

          assert(beginIndices.size() == equationRank);
          assert(endIndices.size() == equationRank);

          if (simulation::getOptions().debug) {
            std::cerr << "    - Begin indices: [";

            size_t rank = beginIndices.size();

            for (size_t dim = 0; dim < rank; ++dim) {
              if (dim != 0) {
                std::cerr << ", ";
              }

              std::cerr << beginIndices[dim];
            }

            std::cerr << "]" << std::endl;
            std::cerr << "      End indices: [";

            for (size_t dim = 0; dim < rank; ++dim) {
              if (dim != 0) {
                std::cerr << ", ";
              }

              std::cerr << endIndices[dim];
            }

            std::cerr << "]" << std::endl;
          }

          std::vector<std::vector<int64_t>> unwrappingBeginIndices;
          std::vector<std::vector<int64_t>> unwrappingEndIndices;

          // We need to detect if some of the dimensions do wrap around.
          // In this case, the indices must be split into multiple ranges.
          std::optional<size_t> increasingDimension = std::nullopt;

          for (size_t dim = 0; dim < equationRank; ++dim) {
            if (endIndices[dim] > beginIndices[dim] &&
                dim + 1 != equationRank) {
              increasingDimension = dim;
              break;
            }
          }

          if (increasingDimension) {
            if (simulation::getOptions().debug) {
              std::cerr << "    - Increasing dimension: "
                        << *increasingDimension << std::endl;
            }

            std::vector<int64_t> currentBeginIndices(beginIndices);
            std::vector<int64_t> currentEndIndices(beginIndices);
            currentEndIndices.back() = equation.indices.back().end - 1;

            unwrappingBeginIndices.push_back(currentBeginIndices);
            unwrappingEndIndices.push_back(currentEndIndices);

            for (size_t i = 0, e = equationRank - *increasingDimension - 2;
                 i < e; ++i) {
              currentBeginIndices[equationRank - i - 1] = 0;

              currentEndIndices[equationRank - i - 2] =
                  equation.indices[equationRank - i - 2].end - 1;

              if (currentBeginIndices[equationRank - i - 2] + 1 !=
                  equation.indices[equationRank - i - 2].end) {
                ++currentBeginIndices[equationRank - i - 2];

                unwrappingBeginIndices.push_back(currentBeginIndices);
                unwrappingEndIndices.push_back(currentEndIndices);
              }
            }

            currentBeginIndices[*increasingDimension + 1] = 0;

            if (endIndices[*increasingDimension] -
                    beginIndices[*increasingDimension] >
                1) {
              ++currentBeginIndices[*increasingDimension];

              currentEndIndices[*increasingDimension] =
                  endIndices[*increasingDimension] - 1;

              unwrappingBeginIndices.push_back(currentBeginIndices);
              unwrappingEndIndices.push_back(currentEndIndices);
            }

            for (size_t i = 0, e = equationRank - *increasingDimension - 1;
                 i < e; ++i) {
              currentBeginIndices[*increasingDimension + i] =
                  endIndices[*increasingDimension + i];

              currentEndIndices[*increasingDimension + i] =
                  endIndices[*increasingDimension + i];

              currentEndIndices[*increasingDimension + i + 1] =
                  endIndices[*increasingDimension + i + 1];

              if (currentEndIndices[*increasingDimension + i + 1] != 0) {
                --currentEndIndices[*increasingDimension + i + 1];
                unwrappingBeginIndices.push_back(currentBeginIndices);
                unwrappingEndIndices.push_back(currentEndIndices);
              }
            }

            currentBeginIndices.back() = endIndices.back();
            currentEndIndices.back() = endIndices.back();
            unwrappingBeginIndices.push_back(currentBeginIndices);
            unwrappingEndIndices.push_back(currentEndIndices);
          } else {
            if (simulation::getOptions().debug) {
              std::cerr << "    - Increasing dimension not found" << std::endl;
            }

            unwrappingBeginIndices.push_back(std::move(beginIndices));
            unwrappingEndIndices.push_back(std::move(endIndices));
          }

          assert(unwrappingBeginIndices.size() == unwrappingEndIndices.size());

          if (simulation::getOptions().debug) {
            for (size_t unwrappingIndex = 0;
                 unwrappingIndex < unwrappingBeginIndices.size();
                 ++unwrappingIndex) {
              std::cerr << "    - #" << unwrappingIndex
                        << " Unwrapping begin indices: [";

              size_t rank = unwrappingBeginIndices[unwrappingIndex].size();

              for (size_t dim = 0; dim < rank; ++dim) {
                if (dim != 0) {
                  std::cerr << ", ";
                }

                std::cerr << unwrappingBeginIndices[unwrappingIndex][dim];
              }

              std::cerr << "]" << std::endl;
              std::cerr << "      #" << unwrappingIndex
                        << " Unwrapping end indices:   [";

              for (size_t dim = 0; dim < rank; ++dim) {
                if (dim != 0) {
                  std::cerr << ", ";
                }

                std::cerr << unwrappingEndIndices[unwrappingIndex][dim];
              }

              std::cerr << "]" << std::endl;
            }
          }

          for (size_t i = 0, e = unwrappingBeginIndices.size(); i < e; ++i) {
            std::vector<int64_t> ranges;

            for (size_t j = 0; j < equationRank; ++j) {
              const auto &currentBeginIndices = unwrappingBeginIndices[i];
              const auto &currentEndIndices = unwrappingEndIndices[i];

              assert(currentBeginIndices[j] <= currentEndIndices[j]);
              ranges.push_back(currentBeginIndices[j]);
              ranges.push_back(currentEndIndices[j] + 1);
            }

            partitionsGroup.emplace_back(equation, ranges);
          }

          // Move to the next partition.
          endFlatIndex =
              getFlatIndex(unwrappingEndIndices.back(), equation.indices);

          equationFlatIndex = endFlatIndex + 1;

          // Create a new equations group if necessary.
          partitionsGroupFlatSize += equationFlatIndex - beginFlatIndex;

          if (partitionsGroupFlatSize >= partitionsGroupMaxFlatSize) {
            pushPartitionsGroupFn();
          }
        }
      } else if (equation.dependencyKind == DependencyKind::Backward) {
        multithreadedSchedule.backwardEquations.emplace_back(equation);
      } else {
        // All the indices must be visited by a single thread.
        std::vector<int64_t> ranges;

        for (const Range &range : equation.indices) {
          ranges.push_back(range.begin);
          ranges.push_back(range.end);
        }

        if (flatSize <= remainingSpace) {
          // There is still space in the current group.
          partitionsGroup.emplace_back(equation, ranges);
          partitionsGroupFlatSize += flatSize;

          if (partitionsGroupFlatSize >= partitionsGroupMaxFlatSize) {
            pushPartitionsGroupFn();
          }
        } else {
          if (flatSize >= partitionsGroupMaxFlatSize) {
            // Independent equations exceeding the maximum number of
            // equations inside a group.
            if (simulation::getOptions().debug) {
              std::cerr << "[Scheduler " << identifier
                        << "] Equation independently exceeds the maximum size "
                           "for a group"
                        << std::endl;
            }

            ContiguousEquationsGroup separateEquationsGroup;
            separateEquationsGroup.emplace_back(equation, ranges);

            multithreadedSchedule.contiguousEquations.push_back(
                std::move(separateEquationsGroup));
          } else {
            pushPartitionsGroupFn();
            partitionsGroup.emplace_back(equation, ranges);
            partitionsGroupFlatSize += flatSize;

            if (partitionsGroupFlatSize >= partitionsGroupMaxFlatSize) {
              pushPartitionsGroupFn();
            }
          }
        }
      }
    }

    if (partitionsGroupFlatSize != 0) {
      pushPartitionsGroupFn();
    }

    assert(std::all_of(equations.begin(), equations.end(),
                       [&](const Equation &equation) {
                         return checkEquationScheduledExactlyOnce(
                             equation,
                             multithreadedSchedule.contiguousEquations);
                       }) &&
           "Not all the equations are scheduled exactly once in the "
           "multithreaded schedule");

    assert(std::all_of(multithreadedSchedule.contiguousEquations.begin(),
                       multithreadedSchedule.contiguousEquations.end(),
                       [&](const ContiguousEquationsGroup &group) {
                         return std::all_of(
                             group.begin(), group.end(),
                             [&](const ContiguousEquationPartition &partition) {
                               return checkEquationIndicesExistence(partition);
                             });
                       }) &&
           "Some nonexistent equation indices have been scheduled in the "
           "multithreaded schedule");
  }

  // Set the execution policy, if forced by the user.
  if (auto forcedPolicy = simulation::getOptions().schedulerPolicy) {
    policy = *forcedPolicy;
  }

  // Print debug information.
  initialized = true;

  if (simulation::getOptions().debug) {
    std::cerr << "[Scheduler " << identifier << "] Initialized" << std::endl;

    if (auto forcedPolicy = simulation::getOptions().schedulerPolicy;
        !forcedPolicy || *forcedPolicy == SchedulerPolicy::Sequential) {
      std::cerr << "  - Sequential schedule size: " << sequentialSchedule.size()
                << std::endl;
    }

    if (auto forcedPolicy = simulation::getOptions().schedulerPolicy;
        !forcedPolicy || *forcedPolicy == SchedulerPolicy::Multithreaded) {
      std::cerr << "  - Multithreaded schedule size: " << std::endl;
      std::cerr << "    - Contiguous equations: "
                << multithreadedSchedule.contiguousEquations.size()
                << std::endl;
      std::cerr << "    - Backward equations: "
                << multithreadedSchedule.backwardEquations.size() << std::endl;
    }
  }

  SCHEDULER_PROFILER_INITIALIZATION_STOP;
}

bool Scheduler::checkEquationScheduledExactlyOnce(
    const Equation &equation,
    const std::vector<ContiguousEquationsGroup> &schedule) const {
  auto beginIndicesIt = MultidimensionalRangeIterator::begin(equation.indices);

  auto endIndicesIt = MultidimensionalRangeIterator::end(equation.indices);

  size_t rank = equation.indices.size();

  for (auto it = beginIndicesIt; it != endIndicesIt; ++it) {
    std::vector<int64_t> indices;

    for (size_t dim = 0; dim < rank; ++dim) {
      indices.push_back((*it)[dim]);
    }

    size_t count = 0;

    for (const ContiguousEquationsGroup &equationsGroup : schedule) {
      count += std::count_if(
          equationsGroup.begin(), equationsGroup.end(),
          [&](const ContiguousEquationPartition &partition) {
            if (partition.first.function != equation.function) {
              return false;
            }

            bool containsPoint = true;

            for (size_t dim = 0; dim < rank && containsPoint; ++dim) {
              if (!(indices[dim] >= partition.second[dim * 2] &&
                    indices[dim] < partition.second[dim * 2 + 1])) {
                containsPoint = false;
              }
            }

            return containsPoint;
          });
    }

    if (count != 1) {
      return false;
    }
  }

  return true;
}

bool Scheduler::checkEquationIndicesExistence(
    const ContiguousEquationPartition &partition) const {
  if (partition.first.dependencyKind == DependencyKind::Backward) {
    // Backward equations are not pre-partitioned.
    return true;
  }

  const auto &equationIndices = partition.first.indices;
  assert(partition.second.size() % 2 == 0);
  size_t rank = partition.second.size() / 2;

  for (size_t dim = 0; dim < rank; ++dim) {
    auto lowerBound = partition.second[dim * 2];
    auto upperBound = partition.second[dim * 2 + 1];

    if (lowerBound < equationIndices[dim].begin) {
      return false;
    }

    if (upperBound > equationIndices[dim].end) {
      return false;
    }
  }

  return true;
}

void Scheduler::run() {
  SCHEDULER_PROFILER_RUN_START;

  if (!initialized) {
    initialize();
  }

  if (policy) {
    // The execution policy has already been determined.
    if (*policy == SchedulerPolicy::Sequential) {
      runSequential();
    } else {
      runMultithreaded();
    }
  } else {
    // Perform calibration.
    int64_t calibrationRuns = simulation::getOptions().schedulerCalibrationRuns;

    bool isSequentialCalibrationRun = runsCounter < calibrationRuns;

    if (isSequentialCalibrationRun) {
      runSequentialWithCalibration();
    } else {
      runMultithreadedWithCalibration();
      bool isLastCalibrationRound = runsCounter == (calibrationRuns * 2 - 1);

      if (isLastCalibrationRound) {
        policy = sequentialRunsMinTime < multithreadedRunsMinTime
                     ? SchedulerPolicy::Sequential
                     : SchedulerPolicy::Multithreaded;

        if (simulation::getOptions().debug) {
          if (policy == SchedulerPolicy::Sequential) {
            std::cerr << "[Scheduler " << identifier
                      << "] Execution policy: sequential" << std::endl;
          } else if (policy == SchedulerPolicy::Multithreaded) {
            std::cerr << "[Scheduler " << identifier
                      << "] Execution policy: multithreaded" << std::endl;
          }
        }
      }
    }
  }

  ++runsCounter;
  SCHEDULER_PROFILER_RUN_STOP;
}

void Scheduler::runSequential() {
  SCHEDULER_PROFILER_INCREMENT_SEQUENTIAL_RUNS_COUNTER;
  SCHEDULER_PROFILER_PARTITIONS_GROUP_START(0);

  for (const ContiguousEquationPartition &partition : sequentialSchedule) {
    const Equation &equation = partition.first;
    const auto &ranges = partition.second;
    equation.function(ranges.data());
  }

  SCHEDULER_PROFILER_PARTITIONS_GROUP_STOP(0);
}

void Scheduler::runSequentialWithCalibration() {
  // Measure the time spent on a sequential computation.
  using namespace std::chrono;

  auto start = steady_clock::now();
  runSequential();
  auto end = steady_clock::now();
  auto elapsed = duration_cast<nanoseconds>(end - start).count();

  if (sequentialRunsMinTime == 0 || elapsed < sequentialRunsMinTime) {
    sequentialRunsMinTime = elapsed;
  }
}

void Scheduler::runMultithreaded() {
  SCHEDULER_PROFILER_INCREMENT_MULTITHREADED_RUNS_COUNTER;

  ThreadPool &threadPool = getSchedulersThreadPool();
  unsigned int numOfThreads = threadPool.getNumOfThreads();

  // Run contiguous equations.
  std::atomic_size_t equationsGroupIndex = 0;

  for (unsigned int thread = 0; thread < numOfThreads; ++thread) {
    threadPool.async([this, thread, &equationsGroupIndex]() {
      size_t assignedEquationsGroup;

      while ((assignedEquationsGroup = equationsGroupIndex++) <
             multithreadedSchedule.contiguousEquations.size()) {
        SCHEDULER_PROFILER_PARTITIONS_GROUP_START(thread);

        const auto &equationsGroup =
            multithreadedSchedule.contiguousEquations[assignedEquationsGroup];

        for (const ContiguousEquationPartition &partition : equationsGroup) {
          const Equation &equation = partition.first;
          const auto &ranges = partition.second;
          equation.function(ranges.data());
        }

        SCHEDULER_PROFILER_PARTITIONS_GROUP_STOP(thread);
      }
    });
  }

  threadPool.wait();

  // Run backward equations.
  std::mutex rowMutex;
  int64_t row;

  // Keeps track of the satisfied dependencies.
  std::vector<std::vector<char>> rowStates;
  rowStates.resize(2 * numOfThreads);

  // Keeps track of which row state a thread is using.
  std::vector<int64_t> threadToRowState;
  threadToRowState.resize(numOfThreads, -1);

  // Keeps track of whether a row is being used by a thread.
  std::vector<bool> rowInUse;
  rowInUse.resize(2 * numOfThreads, false);

  // The identifier of the last assigned row state.
  std::optional<size_t> lastAssignedRowState = std::nullopt;

  auto findAvailableRowState = [&]() -> size_t {
    for (size_t i = 0, e = rowInUse.size(); i < e; ++i) {
      if (!rowInUse[i]) {
        return i;
      }
    }

    assert(false && "No available row state found");
    return -1;
  };

  auto claimRowFn = [&](const Equation &equation,
                        unsigned int threadId) -> std::optional<int64_t> {
    std::lock_guard<std::mutex> lock(rowMutex);

    assert(equation.indices.size() == 2);

    if (row >= equation.indices[1].end) {
      // No more available rows.
      return std::nullopt;
    }

    size_t availableRow = findAvailableRowState();
    threadToRowState[threadId] = availableRow;
    rowInUse[availableRow] = true;

    return row++;
  };

  for (auto &backwardEquation : multithreadedSchedule.backwardEquations) {
    for (auto &rowState : rowStates) {
      size_t numOfColumns = backwardEquation.getEquation().indices[1].end -
                            backwardEquation.getEquation().indices[1].begin;
      rowState.resize(numOfColumns);
    }

    for (unsigned int thread = 0; thread < numOfThreads; ++thread) {
      threadPool.async([&]() {
        std::vector<char> *currentRowState =
            &rowStates[threadToRowState[thread]];

        std::fill(std::begin(*currentRowState), std::end(*currentRowState), 0);

        std::vector<char> *previousRowState = nullptr;

        if (lastAssignedRowState) {
          previousRowState = &rowStates[*lastAssignedRowState];
        }

        backwardEquation.run(thread, claimRowFn, currentRowState,
                             previousRowState);
      });
    }

    row = 0;
  }

  threadPool.wait();
}

void Scheduler::runMultithreadedWithCalibration() {
  // Measure the time spent on a multithreaded computation.
  using namespace std::chrono;

  auto start = steady_clock::now();
  runMultithreaded();
  auto end = steady_clock::now();
  auto elapsed = duration_cast<nanoseconds>(end - start).count();

  if (multithreadedRunsMinTime == 0 || elapsed < multithreadedRunsMinTime) {
    multithreadedRunsMinTime = elapsed;
  }
}
} // namespace marco::runtime

//===---------------------------------------------------------------------===//
// schedulerCreate

[[maybe_unused]] static void *schedulerCreate_pvoid() {
  auto *instance = new Scheduler();
  return static_cast<void *>(instance);
}

RUNTIME_FUNC_DEF(schedulerCreate, PTR(void))

//===---------------------------------------------------------------------===//
// schedulerDestroy

[[maybe_unused]] static void schedulerDestroy_void(void *scheduler) {
  assert(scheduler != nullptr);
  delete static_cast<Scheduler *>(scheduler);
}

RUNTIME_FUNC_DEF(schedulerDestroy, void, PTR(void))

//===---------------------------------------------------------------------===//
// schedulerAddEquation

[[maybe_unused]] static void
schedulerAddEquation_void(void *scheduler, void *equationFunction,
                          uint64_t rank, int64_t *ranges,
                          bool independentIndices) {
  if (independentIndices) {
    static_cast<Scheduler *>(scheduler)->addEquation(
        reinterpret_cast<Scheduler::EquationFunction>(equationFunction), rank,
        ranges, Scheduler::DependencyKind::IndependentIndices);
  } else {
    static_cast<Scheduler *>(scheduler)->addEquation(
        reinterpret_cast<Scheduler::EquationFunction>(equationFunction), rank,
        ranges, Scheduler::DependencyKind::Sequential);
  }
}

[[maybe_unused]] static void schedulerAddEquation_void(void *scheduler,
                                                       void *equationFunction,
                                                       uint64_t rank,
                                                       int64_t *ranges,
                                                       int64_t dependencyKind) {
  assert(scheduler != nullptr);

  Scheduler::DependencyKind dependency = Scheduler::DependencyKind::Sequential;

  if (dependencyKind == 1) {
    dependency = Scheduler::DependencyKind::Backward;
  } else if (dependencyKind == 2) {
    dependency = Scheduler::DependencyKind::IndependentIndices;
  }

  static_cast<Scheduler *>(scheduler)->addEquation(
      reinterpret_cast<Scheduler::EquationFunction>(equationFunction), rank,
      ranges, dependency);
}

RUNTIME_FUNC_DEF(schedulerAddEquation, void, PTR(void), PTR(void), uint64_t,
                 PTR(int64_t), bool)

RUNTIME_FUNC_DEF(schedulerAddEquation, void, PTR(void), PTR(void), uint64_t,
                 PTR(int64_t), int64_t)

//===---------------------------------------------------------------------===//
// schedulerRun

[[maybe_unused]] static void schedulerRun_void(void *scheduler) {
  assert(scheduler != nullptr);
  static_cast<Scheduler *>(scheduler)->run();
}

RUNTIME_FUNC_DEF(schedulerRun, void, PTR(void))
