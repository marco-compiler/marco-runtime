#ifndef MARCO_CODEGEN_TRANSFORMS_EXTERNALSOLVER_H
#define MARCO_CODEGEN_TRANSFORMS_EXTERNALSOLVER_H

#include "marco/Codegen/Transforms/Model/Scheduling.h"
#include "mlir/Transforms/DialectConversion.h"
#include <memory>
#include <vector>

namespace marco::codegen
{
  class IDASolver;

  class ExternalSolver
  {
    public:
      ExternalSolver(mlir::TypeConverter* typeConverter);

      ~ExternalSolver();

      virtual bool isEnabled() const = 0;

      virtual void setEnabled(bool status) = 0;

      virtual bool containsEquation(ScheduledEquation* equation) const = 0;

      virtual mlir::Type getRuntimeDataType(mlir::MLIRContext* context) = 0;

      virtual mlir::LogicalResult processInitFunction(
          mlir::OpBuilder& builder,
          mlir::Value runtimeDataPtr,
          mlir::FuncOp initFunction,
          mlir::ValueRange variables,
          const Model<ScheduledEquationsBlock>& model,
          const mlir::BlockAndValueMapping& derivatives) = 0;

      virtual mlir::LogicalResult processDeinitFunction(
          mlir::OpBuilder& builder,
          mlir::Value runtimeDataPtr,
          mlir::FuncOp deinitFunction) = 0;

      virtual mlir::LogicalResult processUpdateStatesFunction(
          mlir::OpBuilder& builder,
          mlir::Value runtimeDataPtr,
          mlir::FuncOp updateStatesFunction,
          mlir::ValueRange variables,
          const mlir::BlockAndValueMapping& derivatives,
          double requestedTimeStep) = 0;

    protected:
      mlir::TypeConverter* getTypeConverter();

    private:
      mlir::TypeConverter* typeConverter;
  };

  class ExternalSolvers
  {
    private:
      using Container = std::vector<std::unique_ptr<ExternalSolver>>;

    public:
      using iterator = typename Container::iterator;
      using const_iterator = typename Container::const_iterator;

      void addSolver(std::unique_ptr<ExternalSolver> solver);

      bool containsEquation(ScheduledEquation* equation) const;

      iterator begin();
      const_iterator begin() const;

      iterator end();
      const_iterator end() const;

    public:
      Container solvers;
  };
}

#endif // MARCO_CODEGEN_TRANSFORMS_EXTERNALSOLVER_H
