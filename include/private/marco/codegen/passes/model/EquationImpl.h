#ifndef MARCO_CODEGEN_PASSES_MODEL_EQUATIONIMPL_H
#define MARCO_CODEGEN_PASSES_MODEL_EQUATIONIMPL_H

#include "marco/codegen/passes/model/Equation.h"
#include "marco/codegen/passes/model/Path.h"
#include "marco/modeling/AccessFunction.h"

namespace marco::codegen
{
  class BaseEquation : public Equation
  {
    public:
      BaseEquation(modelica::EquationOp equation, Variables variables);

      modelica::EquationOp getOperation() const override;

      const Variables& getVariables() const override;

      void setVariables(Variables variables) override;

      mlir::Value getValueAtPath(const EquationPath& path) const override;

      mlir::LogicalResult explicitate(
          mlir::OpBuilder& builder, const EquationPath& path) override;

      std::unique_ptr<Equation> cloneAndExplicitate(
          mlir::OpBuilder& builder, const EquationPath& path) const override;

      mlir::LogicalResult replaceInto(
          mlir::OpBuilder& builder,
          Equation& destination,
          const ::marco::modeling::AccessFunction& destinationAccessFunction,
          const EquationPath& destinationPath,
          const Access& sourceAccess) const override;

      virtual mlir::FuncOp createTemplateFunction(
          mlir::OpBuilder& builder,
          llvm::StringRef functionName,
          mlir::ValueRange vars,
          ::marco::modeling::scheduling::Direction iterationDirection) const override;

    protected:
      modelica::EquationSidesOp getTerminator() const;

      mlir::LogicalResult explicitate(
          mlir::OpBuilder& builder,
          size_t argumentIndex,
          EquationPath::EquationSide side);

      mlir::LogicalResult groupLeftHandSide(
          mlir::OpBuilder& builder,
          const Access& access);

      std::pair<mlir::Value, std::vector<mlir::Value>> collectSubscriptionIndexes(mlir::Value value) const;

      mlir::Value getMultiplyingFactor(
          mlir::OpBuilder& builder,
          mlir::Value value,
          mlir::Value variable,
          const ::marco::modeling::AccessFunction& accessFunction) const;

      virtual mlir::LogicalResult mapInductionVariables(
          mlir::OpBuilder& builder,
          mlir::BlockAndValueMapping& mapping,
          Equation& destination,
          const ::marco::modeling::AccessFunction& transformation) const = 0;

      virtual mlir::LogicalResult createTemplateFunctionBody(
          mlir::OpBuilder& builder,
          mlir::BlockAndValueMapping& mapping,
          mlir::ValueRange beginIndexes,
          mlir::ValueRange endIndexes,
          mlir::ValueRange steps,
          ::marco::modeling::scheduling::Direction iterationDirection) const = 0;

    private:
      mlir::Operation* equationOp;
      Variables variables;
  };
}

#endif // MARCO_CODEGEN_PASSES_MODEL_EQUATIONIMPL_H
