#ifndef MARCO_CODEGEN_TRANSFORMS_MODEL_MATCHING_H
#define MARCO_CODEGEN_TRANSFORMS_MODEL_MATCHING_H

#include "marco/Codegen/Transforms/Model/Model.h"
#include "marco/Modeling/Dependency.h"
#include <memory>
#include <vector>

namespace marco::codegen
{
  class MatchedEquation : public Equation
  {
    public:
      MatchedEquation(
          std::unique_ptr<Equation> equation,
          modeling::MultidimensionalRange matchedIndexes,
          EquationPath matchedPath);

      MatchedEquation(const MatchedEquation& other);

      ~MatchedEquation();

      MatchedEquation& operator=(const MatchedEquation& other);
      MatchedEquation& operator=(MatchedEquation&& other);

      friend void swap(MatchedEquation& first, MatchedEquation& second);

      std::unique_ptr<Equation> clone() const override;

      /// @name Forwarded methods
      /// {

      mlir::modelica::EquationOp cloneIR() const override;

      void eraseIR() override;

      void dumpIR() const override;

      void dumpIR(llvm::raw_ostream& os) const override;

      mlir::modelica::EquationOp getOperation() const override;

      Variables getVariables() const override;

      void setVariables(Variables variables) override;

      std::vector<Access> getAccesses() const override;

      modeling::DimensionAccess resolveDimensionAccess(
          std::pair<mlir::Value, long> access) const override;

      mlir::Value getValueAtPath(const EquationPath& path) const override;

      mlir::LogicalResult explicitate(
          mlir::OpBuilder& builder, const EquationPath& path) override;

      std::unique_ptr<Equation> cloneIRAndExplicitate(
          mlir::OpBuilder& builder, const EquationPath& path) const override;

      std::vector<mlir::Value> getInductionVariables() const override;

      mlir::LogicalResult replaceInto(
          mlir::OpBuilder& builder,
          Equation& destination,
          const modeling::AccessFunction& destinationAccessFunction,
          const EquationPath& destinationPath) const override;

      mlir::FuncOp createTemplateFunction(
          mlir::OpBuilder& builder,
          llvm::StringRef functionName,
          mlir::ValueRange vars,
          modeling::scheduling::Direction iterationDirection) const override;

      /// }
      /// @name Modified methods
      /// {

      size_t getNumOfIterationVars() const override;

      modeling::MultidimensionalRange getIterationRanges() const override;

      /// }

      std::vector<Access> getReads() const;

      Access getWrite() const;

      std::unique_ptr<Equation> cloneIRAndExplicitate(mlir::OpBuilder& builder) const;

    private:
      std::unique_ptr<Equation> equation;
      modeling::MultidimensionalRange matchedIndexes;
      EquationPath matchedPath;
  };
}

// Traits specializations for the modeling library
namespace marco::modeling::dependency
{
  template<>
  struct EquationTraits<::marco::codegen::MatchedEquation*>
  {
    using Equation = ::marco::codegen::MatchedEquation*;
    using Id = mlir::Operation*;

    static Id getId(const Equation* equation)
    {
      return (*equation)->getOperation().getOperation();
    }

    static size_t getNumOfIterationVars(const Equation* equation)
    {
      return (*equation)->getNumOfIterationVars();
    }

    static MultidimensionalRange getIterationRanges(const Equation* equation)
    {
      return (*equation)->getIterationRanges();
    }

    using VariableType = codegen::Variable*;

    using AccessProperty = codegen::EquationPath;

    static Access<VariableType, AccessProperty> getWrite(const Equation* equation)
    {
      auto write = (*equation)->getWrite();
      return Access(write.getVariable(), write.getAccessFunction(), write.getPath());
    }

    static std::vector<Access<VariableType, AccessProperty>> getReads(const Equation* equation)
    {
      std::vector<Access<VariableType, AccessProperty>> reads;

      for (const auto& read : (*equation)->getReads()) {
        reads.emplace_back(read.getVariable(), read.getAccessFunction(), read.getPath());
      }

      return reads;
    }
  };
}

namespace marco::codegen
{
  /// Match each scalar variable to a scalar equation.
  mlir::LogicalResult match(
      Model<MatchedEquation>& result, Model<Equation>& model, const mlir::BlockAndValueMapping& derivatives);
}

#endif // MARCO_CODEGEN_TRANSFORMS_MODEL_MATCHING_H