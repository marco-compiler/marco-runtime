#ifndef MARCO_MODELING_SCHEDULING_H
#define MARCO_MODELING_SCHEDULING_H

#include "Dependency.h"

namespace marco::modeling
{
  namespace scheduling
  {
    /// The direction to be used by the equations loop to update its iteration variable.
    enum class Direction
    {
      None,     // [i]
      Forward,  // [i + n] with n > 0
      Backward, // [i + n] with n < 0
      Constant, // [n] with n constant value
      Mixed,    // mix of the previous cases
      Unknown
    };
  }

  namespace internal::scheduling
  {

  }

  template<typename VariableProperty, typename EquationProperty>
  class Scheduler
  {
    private:
      using MCIS = internal::MCIS;

      using VectorDependencyGraph = internal::VVarDependencyGraph<VariableProperty, EquationProperty>;
      using Equation = typename VectorDependencyGraph::Equation;
      using SCC = typename VectorDependencyGraph::SCC;
      using SCCDependencyGraph = internal::SCCDependencyGraph<SCC>;

    public:
      bool schedule(llvm::ArrayRef<EquationProperty> equations) const
      {
        VectorDependencyGraph vectorDependencyGraph(equations);
        auto SCCs = vectorDependencyGraph.getSCCs();
        SCCDependencyGraph sccDependencyGraph(SCCs);
        auto scheduledSCCs = sccDependencyGraph.postOrder();

        for (const auto& sccDescriptor : scheduledSCCs) {
          const SCC& scc = sccDependencyGraph[sccDescriptor];
          assert(scc.size() == 1 && "Loop among equations to be scheduled");

          auto accessesDirection = getAccessesDirection(vectorDependencyGraph, scc);

          if (accessesDirection == scheduling::Direction::Forward) {
            dumpSCC(scc);
            std::cout << "Forward\n";
            continue;
          }

          if (accessesDirection == scheduling::Direction::Backward) {
            dumpSCC(scc);
            std::cout << "Backward\n";
            continue;
          }

          std::cout << "Non trivial scheduling";

          const auto& originalGraph = scc.getGraph();

          for (const auto& equationDescriptor : scc) {
            std::cout << originalGraph[equationDescriptor].getId() << " ";
          }

          std::cout << "\n";
        }

        return true;
      }

    private:
      /// Given a SSC containing only one equation that may depend on itself, determine the access direction
      /// with respect to the variable that is written by the equation.
      ///
      /// @param vectorGraph  the graph upon which the SCC has been detected
      /// @param scc          SCC to be examined (consisting of only one equation with a loop on itself)
      /// @return access direction
      scheduling::Direction getAccessesDirection(VectorDependencyGraph& vectorGraph, const SCC& scc) const
      {
        if (!scc.hasCycle()) {
          // If there is no cycle, then the iteration variable of the equation is irrelevant.
          // We prefer the forward direction for simplicity.
          return scheduling::Direction::Forward;
        }

        // If all the dependencies have the same direction, then we can set
        // the induction variable to increase or decrease accordingly.

        const Equation& equation = vectorGraph[scc[0]];
        auto equationRange = equation.getIterationRanges();
        const auto& write = equation.getWrite();
        const auto& writtenVariable = write.getVariable();
        const AccessFunction& writeAccessFunction = write.getAccesFunction();
        MCIS writtenIndexes(writeAccessFunction.map(equationRange));

        auto direction = scheduling::Direction::Unknown;

        for (const auto& read : equation.getReads()) {
          // The access is considered only if it reads the same variable it is being defined by the
          // equation and the ranges overlap.
          const auto& readVariable = read.getVariable();

          if (writtenVariable != readVariable) {
            continue;
          }

          const AccessFunction& readAccessFunction = read.getAccessFunction();
          MCIS readIndexes(readAccessFunction.map(equationRange));

          if (!readIndexes.overlaps(writtenIndexes)) {
            continue;
          }

          // Determine the access direction of the access

          if (!writeAccessFunction.isInvertible()) {
            return scheduling::Direction::Unknown;
          }

          auto relativeAccess = writeAccessFunction.inverse().combine(readAccessFunction);
          auto accessDirection = getAccessFunctionDirection(relativeAccess);
          assert(accessDirection != scheduling::Direction::None && "Algebraic loop detected");

          if (direction == scheduling::Direction::Unknown) {
            direction = accessDirection;
          } else if (direction != accessDirection) {
            return scheduling::Direction::Mixed;
          }
        }

        return direction;
      }

      /// Get the access direction of an access function.
      /// For example, an access consisting in [i0 + 1][i1 + 2] has a forward direction, meaning that
      /// it requires variables that will be defined later in the loop execution. A [i0 - 1][i1 -2]
      /// access function has a backward direction and a [i0 + 1][i1 - 2] has a mixed one.
      /// The indexes of the above induction variables refer to the order in which the induction
      /// variables have been defined, meaning that i0 is the outer-most induction, i1 the second
      /// outer-most one, etc.
      ///
      /// @param accessFunction access function to be analyzed
      /// @return access direction
      scheduling::Direction getAccessFunctionDirection(const AccessFunction& accessFunction) const
      {
        auto direction = scheduling::Direction::Unknown;
        assert(accessFunction.size() != 0);

        for (const auto& dimensionAccess : llvm::enumerate(accessFunction)) {
          auto dimensionDirection = scheduling::Direction::None;

          if (dimensionAccess.value().isConstantAccess()) {
            dimensionDirection = scheduling::Direction::Constant;
          } else {
            if (dimensionAccess.value().getInductionVariableIndex() != dimensionAccess.index()) {
              // If the iteration indexes are out of order, then some accesses will refer to future
              // written variables and others to past written ones.
              return scheduling::Direction::Mixed;
            }

            // Examine the offset of the single dimension access
            auto offset = dimensionAccess.value().getOffset();

            if (offset == 0) {
              dimensionDirection = scheduling::Direction::None;
            } else if (offset > 0) {
              dimensionDirection = scheduling::Direction::Forward;
            } else if (offset < 0) {
              dimensionDirection = scheduling::Direction::Backward;
            }
          }

          if (direction == scheduling::Direction::Unknown) {
            direction = dimensionDirection;
          } else if (direction != dimensionDirection) {
            return scheduling::Direction::Mixed;
          }
        }

        assert(direction != scheduling::Direction::Unknown);
        return direction;
      }

      void dumpSCC(const SCC& scc) const
      {
        const auto& graph = scc.getGraph();

        for (const auto& equation : scc) {
          std::cout << graph[equation].getId() << " ";
        }

        std::cout << "\n";
      }
  };
}

#endif // MARCO_MODELING_SCHEDULING_H
