#ifndef MARCO_RUNTIME_SIMULATION_SCHEDULER_H
#define MARCO_RUNTIME_SIMULATION_SCHEDULER_H

#include "marco/Runtime/Modeling/MultidimensionalRange.h"
#include "marco/Runtime/Profiling/Profiling.h"
#include "marco/Runtime/Profiling/Timer.h"
#include "marco/Runtime/Simulation/Options.h"
#include "marco/Runtime/Support/Mangling.h"
#include <atomic>
#include <cstdint>
#include <variant>
#include <vector>

#ifdef THREADS_ENABLE
#include <condition_variable>
#include <mutex>
#endif

namespace marco::runtime {
#ifdef MARCO_PROFILING
class SchedulerProfiler : public profiling::Profiler {
public:
  SchedulerProfiler(int64_t schedulerId);

  void createPartitionsGroupsCounters(size_t amount);

  void createPartitionsGroupsTimers(size_t amount);

  void reset() override;

  void print() const override;

public:
  profiling::Timer addEquation;
  profiling::Timer initialization;
  profiling::Timer run;
  int64_t sequentialRuns{0};
  int64_t multithreadedRuns{0};
  std::vector<int64_t> partitionsGroupsCounters;
  std::vector<std::unique_ptr<profiling::Timer>> partitionsGroups;

  mutable std::mutex mutex;
};
#endif

//===---------------------------------------------------------------------===//
// ReadyState
//===---------------------------------------------------------------------===//
struct ReadyState {
  explicit ReadyState() { flag.store(false); }

  explicit ReadyState(bool flag) { this->flag.store(flag); }

  ReadyState(const ReadyState &other) { *this = other; }
  ReadyState(ReadyState &&other) { *this = std::move(other); }

  ReadyState &operator=(const ReadyState &other) {
    this->flag.store(other.flag.load());
    return *this;
  }

  ReadyState &operator=(bool flag) {
    this->flag.store(flag);
    return *this;
  }

  ReadyState &operator=(ReadyState &&other) {
    const bool val = other.flag.load();
    this->flag.store(val);
    other.flag.store(false);
    return *this;
  }

  bool operator!=(bool other) { return flag.load() != other; }

  bool operator==(bool other) { return flag.load() == other; }

  bool isReady() { return *this == true; }

  void setReady() { flag.store(true); }

  void unsetReady() { flag.store(false); }

  bool compareExchange(bool expect, bool desired) {
    return flag.compare_exchange_weak(expect, desired);
  }

  std::atomic<bool> flag;
};

class Scheduler {
public:
  using EquationFunction = void (*)(const int64_t *);

  enum class DependencyKind { Sequential, Backward, IndependentIndices };

  struct Equation {
    EquationFunction function;
    MultidimensionalRange indices;
    DependencyKind dependencyKind;

    Equation(EquationFunction function, MultidimensionalRange indices,
             DependencyKind dependencyKind);
  };

  // An equation partition is composed of:
  //   - the equation descriptor.
  //   - the ranges information to be passed to the equation function.
  class EquationPartition {
    const Equation *equation;
    std::vector<int64_t> ranges;

    void run();
  };

  // A contiguous equation partition is composed of:
  //   - the equation descriptor.
  //   - the ranges information to be passed to the equation function.
  using ContiguousEquationPartition = std::pair<Equation, std::vector<int64_t>>;

  class BackwardEquation {
  private:
    const Equation *equation;

    // Keeps track of the satisfied dependencies.
    std::vector<std::vector<ReadyState>> readyStates;

    // Keeps track of whether a row is being used by a thread.
    std::vector<ReadyState> rowInUse;

    // The next row to be computed.
    std::mutex rowMutex;
    int64_t nextRow{0};

    // The identifier of the last assigned ready state row.
    std::optional<size_t> lastAssignedRowState{std::nullopt};

    struct Job {
      int64_t row;
      size_t currentRowState;
      std::optional<size_t> previousRowState;
    };

  public:
    explicit BackwardEquation(const Equation &equation, uint64_t numOfThreads);

    BackwardEquation(const BackwardEquation &other) = delete;

    BackwardEquation(BackwardEquation &&other);

    ~BackwardEquation();

    BackwardEquation &operator=(const BackwardEquation &other) = delete;

    BackwardEquation &operator=(BackwardEquation &&other);

    const Equation &getEquation() const;

    void run();

    void reset();

  private:
    MultidimensionalRange getRowBoundaries(int64_t row) const;

    std::vector<int64_t> getBeginIndices(int64_t row) const;

    std::optional<Job> claimJob();

    size_t getAvailableRowState();
  };

  using SequentialSchedule = std::vector<ContiguousEquationPartition>;

  // A group of contiguous equation partitions.
  using ContiguousEquationsGroup = std::vector<ContiguousEquationPartition>;

  struct MultithreadedSchedule {
    std::vector<ContiguousEquationsGroup> contiguousEquations;
    std::vector<BackwardEquation> backwardEquations;
  };

  Scheduler();

  void addEquation(EquationFunction function, uint64_t rank, int64_t *ranges,
                   DependencyKind dependencyKind);

  void run();

private:
  void initialize();

  [[maybe_unused, nodiscard]] bool checkEquationScheduledExactlyOnce(
      const Equation &equation, const MultithreadedSchedule &schedule) const;

  [[maybe_unused, nodiscard]] bool checkEquationScheduledExactlyOnce(
      const Equation &equation,
      const std::vector<ContiguousEquationsGroup> &schedule) const;

  [[maybe_unused, nodiscard]] bool checkEquationScheduledExactlyOnce(
      const Equation &equation,
      const std::vector<BackwardEquation> &schedule) const;

  [[maybe_unused, nodiscard]] bool checkEquationIndicesExistence(
      const ContiguousEquationPartition &equationPartition) const;

  void runSequential();

  void runSequentialWithCalibration();

  void runMultithreaded();

  void runMultithreadedWithCalibration();

private:
  int64_t identifier{0};
  bool initialized{false};
  std::vector<Equation> equations;

  // The list of equation partitions to be executed in case of sequential
  // execution policy.
  // The information is computed only once during the initialization.
  SequentialSchedule sequentialSchedule;

  // The list of equations groups the threads will process. Each thread
  // processes one group at a time.
  // The information is computed only once during the initialization.
  MultithreadedSchedule multithreadedSchedule;

  int64_t runsCounter{0};
  int64_t sequentialRunsMinTime{0};
  int64_t multithreadedRunsMinTime{0};
  std::optional<SchedulerPolicy> policy{std::nullopt};

#ifdef MARCO_PROFILING
  // Profiling.
  std::shared_ptr<SchedulerProfiler> profiler;
#endif
};
} // namespace marco::runtime

//===---------------------------------------------------------------------===//
// Exported functions
//===---------------------------------------------------------------------===//

RUNTIME_FUNC_DECL(schedulerCreate, PTR(void))

RUNTIME_FUNC_DECL(schedulerDestroy, void, PTR(void))

RUNTIME_FUNC_DECL(schedulerAddEquation, void, PTR(void), PTR(void), uint64_t,
                  PTR(int64_t), bool)

RUNTIME_FUNC_DECL(schedulerAddEquation, void, PTR(void), PTR(void), uint64_t,
                  PTR(int64_t), uint32_t)

RUNTIME_FUNC_DECL(schedulerRun, void, PTR(void))

#endif // MARCO_RUNTIME_SIMULATION_SCHEDULER_H
