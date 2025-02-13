#include "marco/Runtime/Simulation/Scheduler.h"
#include "marco/Runtime/Multithreading/Barrier.h"
#include "marco/Runtime/Multithreading/ThreadPool.h"
#include "marco/Runtime/Simulation/Options.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <cstring>
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

void SchedulerProfiler::reset() {
  addEquation.reset();
  initialization.reset();
  run.reset();
  sequentialSchedule.reset();
  multithreadedSchedule.reset();
}

void SchedulerProfiler::SequentialScheduleProfiling::reset() { executions = 0; }

void SchedulerProfiler::MultithreadedScheduleProfiling::reset() {
  executions = 0;

  for (auto &profiling : contiguousEquations) {
    profiling.reset();
  }

  for (auto &profiling : backwardEquations) {
    profiling.reset();
  }
}

void SchedulerProfiler::MultithreadedScheduleProfiling::
    ContiguousEquationsGroupsProfiling::reset() {
  groupsCounter = 0;
  timer->reset();
}

void SchedulerProfiler::MultithreadedScheduleProfiling::
    BackwardEquationProfiling::reset() {
  equationFunction->reset();
  dependencyWaiting->reset();
}

void SchedulerProfiler::MultithreadedScheduleProfiling::setNumOfThreads(
    uint64_t numOfThreads) {
  contiguousEquations.clear();
  backwardEquations.clear();

  for (size_t i = 0; i < numOfThreads; ++i) {
    contiguousEquations.emplace_back(ContiguousEquationsGroupsProfiling{
        0, std::make_unique<profiling::Timer>()});

    backwardEquations.emplace_back(
        BackwardEquationProfiling{0, std::make_unique<profiling::Timer>(),
                                  std::make_unique<profiling::Timer>()});
  }
}

void SchedulerProfiler::print() const {
  std::cerr << "Time spent on adding the equations: "
            << addEquation.totalElapsedTime<std::milli>() << " ms\n";

  std::cerr << "Time spent on initialization: "
            << initialization.totalElapsedTime<std::milli>() << " ms\n";

  std::cerr << "Time spent on 'run' method: "
            << run.totalElapsedTime<std::milli>() << " ms\n";

  if (sequentialSchedule.executions != 0) {
    std::cerr << "\n";
    sequentialSchedule.print();
  }

  if (multithreadedSchedule.executions != 0) {
    std::cerr << "\n";
    multithreadedSchedule.print();
  }
}

void SchedulerProfiler::SequentialScheduleProfiling::print() const {
  std::cerr << "Sequential schedule:\n";
  std::cerr << "  - Number of executions: " << executions << "\n";
}

void SchedulerProfiler::MultithreadedScheduleProfiling::print() const {
  std::cerr << "Multithreaded schedule:\n";
  std::cerr << "  - Number of executions: " << executions << "\n\n";

  std::cerr << "  - Contiguous equation groups:\n";

  for (size_t threadId = 0, e = contiguousEquations.size(); threadId < e;
       ++threadId) {
    auto contiguousGroupsCounter = contiguousEquations[threadId].groupsCounter;

    double averageContiguousGroupTime = 0;

    if (contiguousGroupsCounter != 0) {
      averageContiguousGroupTime =
          contiguousEquations[threadId].timer->totalElapsedTime<std::nano>() /
          static_cast<double>(contiguousGroupsCounter);
    }

    std::cerr << "    - Thread " << threadId << "\n";

    std::cerr
        << "      "
        << "Total time spent processing contiguous equations: "
        << contiguousEquations[threadId].timer->totalElapsedTime<std::milli>()
        << " ms\n";

    std::cerr << "      "
              << "Average time spent processing contiguous equation groups: "
              << averageContiguousGroupTime << " ns\n";

    std::cerr << "      "
              << "Number of contiguous equation groups processed: "
              << contiguousGroupsCounter << "\n";

    std::cerr << "\n";
  }

  std::cerr << "  - Backward equations:\n";

  for (size_t threadId = 0, e = contiguousEquations.size(); threadId < e;
       ++threadId) {
    std::cerr << "    - Thread " << threadId << "\n";

    auto rowCounter = backwardEquations[threadId].rowCounter;
    std::cerr << "      " << "Number of rows processed: " << rowCounter << "\n";

    std::cerr << "      "
              << "Total time spent executing equation functions: "
              << backwardEquations[threadId]
                     .equationFunction->totalElapsedTime<std::milli>()
              << " ms\n";

    double averageEquationFunctionRowTime = 0;

    if (rowCounter != 0) {
      averageEquationFunctionRowTime =
          backwardEquations[threadId]
              .equationFunction->totalElapsedTime<std::nano>() /
          static_cast<double>(rowCounter);
    }

    std::cerr << "      "
              << "Average time per row spent executing equation functions: "
              << averageEquationFunctionRowTime << " ns\n";

    std::cerr << "      "
              << "Total time spent waiting for dependencies: "
              << backwardEquations[threadId]
                     .dependencyWaiting->totalElapsedTime<std::milli>()
              << " ms\n";

    double averageWaitRowTime = 0;

    if (rowCounter != 0) {
      averageWaitRowTime =
          backwardEquations[threadId]
              .dependencyWaiting->totalElapsedTime<std::nano>() /
          static_cast<double>(rowCounter);
    }

    std::cerr << "      "
              << "Average time per row spent waiting for dependencies: "
              << averageWaitRowTime << " ns\n";

    std::cerr << "\n";
  }
}

// clang-format off
#define SCHEDULER_PROFILER_ADD_EQUATION_START                                 \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    profiler->addEquation.start();                                            \
  }

#define SCHEDULER_PROFILER_ADD_EQUATION_STOP                                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    profiler->addEquation.stop();                                             \
  }

#define SCHEDULER_PROFILER_RUN_START                                          \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    profiler->run.start();                                                    \
  }

#define SCHEDULER_PROFILER_RUN_STOP                                           \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    profiler->run.stop();                                                     \
  }

#define SCHEDULER_PROFILER_INITIALIZATION_START                               \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    profiler->initialization.start();                                         \
  }

#define SCHEDULER_PROFILER_INITIALIZATION_STOP                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    profiler->initialization.stop();                                          \
  }

#define SCHEDULER_PROFILER_INCREMENT_SEQUENTIAL_EXECUTIONS_COUNTER            \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    ++profiler->sequentialSchedule.executions;                                \
  }

#define SCHEDULER_PROFILER_INCREMENT_MULTITHREADED_EXECUTIONS_COUNTER         \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    ++profiler->multithreadedSchedule.executions;                             \
  }

#define SCHEDULER_PROFILER_MT_CONTIGUOUS_GROUP_START(thread)                     \
  if (::marco::runtime::simulation::getOptions().profiling) {                    \
    assert(profiler != nullptr);                                                 \
    ++profiler->multithreadedSchedule.contiguousEquations[thread].groupsCounter; \
    profiler->multithreadedSchedule.contiguousEquations[thread].timer->start();  \
  }

#define SCHEDULER_PROFILER_MT_CONTIGUOUS_GROUP_STOP(thread)                    \
  if (::marco::runtime::simulation::getOptions().profiling) {                  \
    assert(profiler != nullptr);                                               \
    profiler->multithreadedSchedule.contiguousEquations[thread].timer->stop(); \
  }

#define SCHEDULER_PROFILER_MT_BACKWARD_ROW_INCREMENT(thread)                  \
  if (::marco::runtime::simulation::getOptions().profiling) {                 \
    assert(profiler != nullptr);                                              \
    ++profiler->multithreadedSchedule.backwardEquations[thread].rowCounter;   \
  }

#define SCHEDULER_PROFILER_MT_BACKWARD_FUNC_START(thread)                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                            \
    assert(profiler != nullptr);                                                         \
    profiler->multithreadedSchedule.backwardEquations[thread].equationFunction->start(); \
  }

#define SCHEDULER_PROFILER_MT_BACKWARD_FUNC_STOP(thread)                                \
  if (::marco::runtime::simulation::getOptions().profiling) {                           \
    assert(profiler != nullptr);                                                        \
    profiler->multithreadedSchedule.backwardEquations[thread].equationFunction->stop(); \
  }

#define SCHEDULER_PROFILER_MT_BACKWARD_DEPENDENCY_START(thread)                           \
  if (::marco::runtime::simulation::getOptions().profiling) {                             \
    assert(profiler != nullptr);                                                          \
    profiler->multithreadedSchedule.backwardEquations[thread].dependencyWaiting->start(); \
  }

#define SCHEDULER_PROFILER_MT_BACKWARD_DEPENDENCY_STOP(thread)                           \
  if (::marco::runtime::simulation::getOptions().profiling) {                            \
    assert(profiler != nullptr);                                                         \
    profiler->multithreadedSchedule.backwardEquations[thread].dependencyWaiting->stop(); \
  }
// clang-format on

#else

// clang-format off
#define SCHEDULER_PROFILER_DO_NOTHING static_assert(true);

#define SCHEDULER_PROFILER_ADD_EQUATION_START SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_ADD_EQUATION_STOP SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_RUN_START SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_RUN_STOP SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_INITIALIZATION_START SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_INITIALIZATION_STOP SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_INCREMENT_SEQUENTIAL_EXECUTIONS_COUNTER SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_INCREMENT_MULTITHREADED_EXECUTIONS_COUNTER SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_MT_CONTIGUOUS_GROUP_START(thread) SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_MT_CONTIGUOUS_GROUP_STOP(thread) SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_MT_BACKWARD_ROW_INCREMENT(thread) SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_MT_BACKWARD_FUNC_START(thread) SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_MT_BACKWARD_FUNC_STOP(thread) SCHEDULER_PROFILER_DO_NOTHING

#define SCHEDULER_PROFILER_MT_BACKWARD_DEPENDENCY_START(thread) SCHEDULER_PROFILER_DO_NOTHING
#define SCHEDULER_PROFILER_MT_BACKWARD_DEPENDENCY_STOP(thread) SCHEDULER_PROFILER_DO_NOTHING
// clang-format on

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
    const Scheduler::Equation &equation, uint64_t numOfThreads,
    SchedulerProfiler *profiler)
    : equation(&equation), profiler(profiler) {
  readyStates.resize(2 * numOfThreads);
  size_t numOfColumns = equation.indices[1].end - equation.indices[1].begin;

  for (auto &readyState : readyStates) {
    readyState.resize(numOfColumns, ReadyState(false));
  }

  rowInUse.resize(2 * numOfThreads, ReadyState(false));
}

Scheduler::BackwardEquation::BackwardEquation(BackwardEquation &&other) {
  std::unique_lock<std::mutex> lock(other.rowMutex);
  equation = other.equation;
  profiler = other.profiler;
  readyStates = std::move(other.readyStates);
  rowInUse = std::move(other.rowInUse);
  nextRow = other.nextRow;
  lastAssignedRowState = other.lastAssignedRowState;
}

Scheduler::BackwardEquation::~BackwardEquation() = default;

Scheduler::BackwardEquation &
Scheduler::BackwardEquation::operator=(BackwardEquation &&other) {
  if (this != &other) {
    std::unique_lock<std::mutex> lhsLock(rowMutex);
    std::unique_lock<std::mutex> rhsLock(other.rowMutex);

    equation = other.equation;
    profiler = other.profiler;
    readyStates = std::move(other.readyStates);
    rowInUse = std::move(other.rowInUse);
    nextRow = other.nextRow;
    lastAssignedRowState = other.lastAssignedRowState;
  }

  return *this;
}

const Scheduler::Equation &Scheduler::BackwardEquation::getEquation() const {
  assert(equation != nullptr);
  return *equation;
}

void Scheduler::BackwardEquation::run(unsigned int threadId) {
  std::optional<Job> job = claimJob();

  while (job) {
    MultidimensionalRange rowRange = getRowBoundaries(job->row);
    std::vector<int64_t> indices = getBeginIndices(job->row);
    size_t rank = indices.size();

    std::vector<int64_t> functionRanges;
    functionRanges.resize(rank * 2, 0xDEADBEEF);

    do {
      SCHEDULER_PROFILER_MT_BACKWARD_ROW_INCREMENT(threadId)

      // Compute the iteration boundaries that the equation function will use.
      for (size_t dim = 0; dim < rank; ++dim) {
        functionRanges[dim * 2] = indices[dim];
        functionRanges[dim * 2 + 1] = indices[dim] + 1;
      }

      int64_t column = indices[1];

      if (job->previousRowState) {
        SCHEDULER_PROFILER_MT_BACKWARD_DEPENDENCY_START(threadId)

        while (!readyStates[*job->previousRowState][column].isReady()) {
          // Intentionally spin-wait until the dependencies are satisfied.
          // Locking a mutex is expensive, and the waiting time tends to zero as
          // the iteration space is progressively visited.
        }

        SCHEDULER_PROFILER_MT_BACKWARD_DEPENDENCY_STOP(threadId)
      }

      // Call the equation function on the specific grid cell.
      SCHEDULER_PROFILER_MT_BACKWARD_FUNC_START(threadId)
      equation->function(functionRanges.data());
      SCHEDULER_PROFILER_MT_BACKWARD_FUNC_STOP(threadId)

      // Mark the cell as processed, so that other threads can continue.
      readyStates[job->currentRowState][column].setReady();
    } while (advanceEquationIndices(indices, rowRange));

    // Move to the next row.
    if (job->previousRowState) {
      rowInUse[*job->previousRowState] = false;
    }

    job = claimJob();
  }
}

void Scheduler::BackwardEquation::reset() {
  for (auto &readyStateRow : readyStates) {
    std::memset(readyStateRow.data(), false, readyStateRow.size());
  }

  for (size_t i = 0, e = rowInUse.size(); i < e; ++i) {
    rowInUse[i] = false;
  }

  nextRow = 0;
  lastAssignedRowState = std::nullopt;
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

std::optional<Scheduler::BackwardEquation::Job>
Scheduler::BackwardEquation::claimJob() {
  std::lock_guard<std::mutex> lock(rowMutex);
  assert(equation->indices.size() == 2);

  if (nextRow + equation->indices[0].begin >= equation->indices[0].end) {
    // No more available rows.
    return std::nullopt;
  }

  size_t availableRow = getAvailableRowState();

  std::optional<size_t> previousReadyState = lastAssignedRowState;
  lastAssignedRowState = availableRow;

  std::memset(readyStates[availableRow].data(), false,
              readyStates[availableRow].size());

  return Job{nextRow++, availableRow, previousReadyState};
}

size_t Scheduler::BackwardEquation::getAvailableRowState() {
  while (true) {
    for (size_t i = 0, e = rowInUse.size(); i < e; ++i) {
      if (rowInUse[i].compareExchange(false, true)) {
        return i;
      }
    }
  }
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
    profiler->multithreadedSchedule.setNumOfThreads(numOfThreads);

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
        std::cerr << "  - Dependency kind: ";

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
        multithreadedSchedule.backwardEquations.emplace_back(
            equation, numOfThreads, profiler.get());
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
                             equation, multithreadedSchedule);
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
    const Equation &equation, const MultithreadedSchedule &schedule) const {
  switch (equation.dependencyKind) {
  case DependencyKind::Sequential:
  case DependencyKind::IndependentIndices:
    return checkEquationScheduledExactlyOnce(equation,
                                             schedule.contiguousEquations);

  case DependencyKind::Backward:
    return checkEquationScheduledExactlyOnce(equation,
                                             schedule.backwardEquations);
  }

  return false;
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

bool Scheduler::checkEquationScheduledExactlyOnce(
    const Equation &equation,
    const std::vector<BackwardEquation> &schedule) const {
  for (const BackwardEquation &backwardEquation : schedule) {
    if (backwardEquation.getEquation().function == equation.function) {
      return true;
    }
  }

  return false;
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
  SCHEDULER_PROFILER_INCREMENT_SEQUENTIAL_EXECUTIONS_COUNTER;

  for (const ContiguousEquationPartition &partition : sequentialSchedule) {
    const Equation &equation = partition.first;
    const auto &ranges = partition.second;
    equation.function(ranges.data());
  }
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
  SCHEDULER_PROFILER_INCREMENT_MULTITHREADED_EXECUTIONS_COUNTER;

  ThreadPool &threadPool = getSchedulersThreadPool();
  unsigned int numOfThreads = threadPool.getNumOfThreads();

  // Run the contiguous equations.
  std::atomic_size_t equationsGroupIndex = 0;

  for (unsigned int thread = 0; thread < numOfThreads; ++thread) {
    threadPool.async([this, thread, &equationsGroupIndex]() {
      size_t assignedEquationsGroup;

      while ((assignedEquationsGroup = equationsGroupIndex++) <
             multithreadedSchedule.contiguousEquations.size()) {
        SCHEDULER_PROFILER_MT_CONTIGUOUS_GROUP_START(thread);

        const auto &equationsGroup =
            multithreadedSchedule.contiguousEquations[assignedEquationsGroup];

        for (const ContiguousEquationPartition &partition : equationsGroup) {
          const Equation &equation = partition.first;
          const auto &ranges = partition.second;
          equation.function(ranges.data());
        }

        SCHEDULER_PROFILER_MT_CONTIGUOUS_GROUP_STOP(thread);
      }
    });
  }

  threadPool.wait();

  // Run the backward equations.
  for (auto &backwardEquation : multithreadedSchedule.backwardEquations) {
    for (unsigned int thread = 0; thread < numOfThreads; ++thread) {
      threadPool.async(
          [&backwardEquation, thread]() { backwardEquation.run(thread); });
    }

    threadPool.wait();
    backwardEquation.reset();
  }
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

[[maybe_unused]] static void
schedulerAddEquation_void(void *scheduler, void *equationFunction,
                          uint64_t rank, int64_t *ranges,
                          uint32_t dependencyKind) {
  assert(scheduler != nullptr);

  Scheduler::DependencyKind dependency = Scheduler::DependencyKind::Sequential;

  if (dependencyKind == 1) {
    dependency = Scheduler::DependencyKind::IndependentIndices;
  } else if (dependencyKind == 2) {
    dependency = Scheduler::DependencyKind::Backward;
  }

  static_cast<Scheduler *>(scheduler)->addEquation(
      reinterpret_cast<Scheduler::EquationFunction>(equationFunction), rank,
      ranges, dependency);
}

RUNTIME_FUNC_DEF(schedulerAddEquation, void, PTR(void), PTR(void), uint64_t,
                 PTR(int64_t), bool)

RUNTIME_FUNC_DEF(schedulerAddEquation, void, PTR(void), PTR(void), uint64_t,
                 PTR(int64_t), uint32_t)

//===---------------------------------------------------------------------===//
// schedulerRun

[[maybe_unused]] static void schedulerRun_void(void *scheduler) {
  assert(scheduler != nullptr);
  static_cast<Scheduler *>(scheduler)->run();
}

RUNTIME_FUNC_DEF(schedulerRun, void, PTR(void))
