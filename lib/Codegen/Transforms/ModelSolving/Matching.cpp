#include "marco/Codegen/Transforms/ModelSolving/Matching.h"
#include "marco/Codegen/Transforms/ModelSolving/FilteredVariable.h"
#include "marco/Dialect/Modelica/ModelicaDialect.h"

using namespace ::marco::codegen;
using namespace ::marco::modeling;
using namespace ::mlir::modelica;

namespace marco::codegen
{
  MatchedEquation::MatchedEquation(
      std::unique_ptr<Equation> equation,
      modeling::MultidimensionalRange matchedIndexes,
      EquationPath matchedPath)
    : equation(std::move(equation)),
      matchedIndexes(std::move(matchedIndexes)),
      matchedPath(std::move(matchedPath))
  {
    assert(this->equation->getIterationRanges().contains(this->matchedIndexes));
  }

  MatchedEquation::MatchedEquation(const MatchedEquation& other)
    : equation(other.equation->clone()),
      matchedIndexes(other.matchedIndexes),
      matchedPath(other.matchedPath)
  {
  }

  MatchedEquation::~MatchedEquation() = default;

  MatchedEquation& MatchedEquation::operator=(const MatchedEquation& other)
  {
    MatchedEquation result(other);
    swap(*this, result);
    return *this;
  }

  MatchedEquation& MatchedEquation::operator=(MatchedEquation&& other) = default;

  void swap(MatchedEquation& first, MatchedEquation& second)
  {
    using std::swap;
    swap(first.equation, second.equation);
    swap(first.matchedIndexes, second.matchedIndexes);
    swap(first.matchedPath, second.matchedPath);
  }

  std::unique_ptr<Equation> MatchedEquation::clone() const
  {
    return std::make_unique<MatchedEquation>(*this);
  }

  EquationOp MatchedEquation::cloneIR() const
  {
    return equation->cloneIR();
  }

  void MatchedEquation::eraseIR()
  {
    equation->eraseIR();
  }

  void MatchedEquation::dumpIR() const
  {
    equation->dumpIR();
  }

  void MatchedEquation::dumpIR(llvm::raw_ostream& os) const
  {
    equation->dumpIR(os);
  }

  EquationOp MatchedEquation::getOperation() const
  {
    return equation->getOperation();
  }

  Variables MatchedEquation::getVariables() const
  {
    return equation->getVariables();
  }

  void MatchedEquation::setVariables(Variables variables)
  {
    equation->setVariables(std::move(variables));
  }

  std::vector<Access> MatchedEquation::getAccesses() const
  {
    return equation->getAccesses();
  }

  ::marco::modeling::DimensionAccess MatchedEquation::resolveDimensionAccess(
      std::pair<mlir::Value, long> access) const
  {
    return equation->resolveDimensionAccess(std::move(access));
  }

  mlir::FuncOp MatchedEquation::createTemplateFunction(
      mlir::OpBuilder& builder,
      llvm::StringRef functionName,
      mlir::ValueRange vars,
      ::marco::modeling::scheduling::Direction iterationDirection) const
  {
    return equation->createTemplateFunction(builder, functionName, vars, iterationDirection);
  }

  mlir::Value MatchedEquation::getValueAtPath(const EquationPath& path) const
  {
    return equation->getValueAtPath(path);
  }

  Access MatchedEquation::getAccessAtPath(const EquationPath& path) const
  {
    return equation->getAccessAtPath(path);
  }

  mlir::LogicalResult MatchedEquation::explicitate(
      mlir::OpBuilder& builder,
      const MultidimensionalRange& equationIndices,
      const EquationPath& path)
  {
    return equation->explicitate(builder, equationIndices, path);
  }

  std::unique_ptr<Equation> MatchedEquation::cloneIRAndExplicitate(
      mlir::OpBuilder& builder,
      const MultidimensionalRange& equationIndices,
      const EquationPath& path) const
  {
    return equation->cloneIRAndExplicitate(builder, equationIndices, path);
  }

  std::vector<mlir::Value> MatchedEquation::getInductionVariables() const
  {
    return equation->getInductionVariables();
  }

  mlir::LogicalResult MatchedEquation::replaceInto(
      mlir::OpBuilder& builder,
      const MultidimensionalRange& equationIndices,
      Equation& destination,
      const ::marco::modeling::AccessFunction& destinationAccessFunction,
      const EquationPath& destinationPath) const
  {
    return equation->replaceInto(builder, equationIndices, destination, destinationAccessFunction, destinationPath);
  }

  size_t MatchedEquation::getNumOfIterationVars() const
  {
    return matchedIndexes.rank();
  }

  modeling::MultidimensionalRange MatchedEquation::getIterationRanges() const
  {
    return matchedIndexes;
  }

  std::vector<Access> MatchedEquation::getReads() const
  {
    std::vector<Access> result;

    auto iterationRanges = getIterationRanges();

    auto writeAccess = getWrite();
    auto writtenVariable = writeAccess.getVariable();
    auto writtenIndices = writeAccess.getAccessFunction().map(iterationRanges);

    for (const auto& access : getAccesses()) {
      auto accessedVariable = access.getVariable();

      if (accessedVariable != writtenVariable) {
        result.push_back(access);
      } else {
        auto accessedIndices = access.getAccessFunction().map(iterationRanges);

        if (accessedIndices != writtenIndices) {
          result.push_back(access);
        }
      }
    }

    return result;
  }

  Access MatchedEquation::getWrite() const
  {
    return getAccessAtPath(matchedPath);
  }

  std::unique_ptr<Equation> MatchedEquation::cloneIRAndExplicitate(
      mlir::OpBuilder& builder,
      const MultidimensionalRange& equationIndices) const
  {
    return equation->cloneIRAndExplicitate(builder, equationIndices, getWrite().getPath());
  }

  std::unique_ptr<Equation> MatchedEquation::cloneIRAndExplicitate(mlir::OpBuilder& builder) const
  {
    return equation->cloneIRAndExplicitate(builder, getIterationRanges(), getWrite().getPath());
  }

  mlir::LogicalResult match(
      Model<MatchedEquation>& result,
      const Model<Equation>& model,
      std::function<IndexSet(const Variable&)> matchableIndicesFn)
  {
    Variables allVariables = model.getVariables();

    // Map the variables by their argument number for a faster lookup
    std::vector<mlir::BlockArgument> variablesByPosition;
    variablesByPosition.resize(allVariables.size());

    for (const auto& variable : allVariables) {
      auto argument = variable->getValue().cast<mlir::BlockArgument>();
      variablesByPosition[argument.getArgNumber()] = argument;
    }

    // Filter the variables. State and constant ones must not in fact
    // take part into the matching process as their values are already
    // determined (state variables depend on their derivatives, while
    // constants have a fixed value).

    Variables filteredVariables;

    for (const auto& variable : allVariables) {
      auto matchableIndices = matchableIndicesFn(*variable);

      for (const auto& range : matchableIndices) {
        auto filteredVariable = std::make_unique<FilteredVariable>(variable->clone(), IndexSet(range));
        filteredVariables.add(std::move(filteredVariable));
      }
    }

    Equations<Equation> filteredEquations;

    for (const auto& equation : model.getEquations()) {
      auto clone = equation->clone();
      clone->setVariables(filteredVariables);
      filteredEquations.add(std::move(clone));
    }

    Model<Equation> filteredModel(model.getOperation());
    filteredModel.setVariables(filteredVariables);
    filteredModel.setEquations(filteredEquations);

    // Create the matching graph. We use the pointers to the real nodes in order
    // to speed up the copies.
    MatchingGraph<Variable*, Equation*> matchingGraph;

    for (const auto& variable : filteredModel.getVariables()) {
      matchingGraph.addVariable(variable.get());
    }

    for (const auto& equation : filteredModel.getEquations()) {
      matchingGraph.addEquation(equation.get());
    }

    auto numberOfScalarEquations = matchingGraph.getNumberOfScalarEquations();
    auto numberOfScalarVariables = matchingGraph.getNumberOfScalarVariables();

    if (numberOfScalarEquations < numberOfScalarVariables) {
      model.getOperation().emitError(
          "Underdetermined model. Found " +
          std::to_string(numberOfScalarEquations) +
          " scalar equations and " +
          std::to_string(numberOfScalarVariables) +
          " scalar variables.");

      return mlir::failure();
    } else if (numberOfScalarEquations > numberOfScalarVariables) {
      model.getOperation().emitError(
          "Overdetermined model. Found " +
          std::to_string(numberOfScalarEquations) +
          " scalar equations and " +
          std::to_string(numberOfScalarVariables) +
          " scalar variables.");

      return mlir::failure();
    }

    // Apply the simplification algorithm to solve the obliged matches
    if (!matchingGraph.simplify()) {
      model.getOperation().emitError("Inconsistency found during the matching simplification process");
      return mlir::failure();
    }

    // Apply the full matching algorithm for the equations and variables that are still unmatched
    if (!matchingGraph.match()) {
      model.getOperation().emitError("Matching failed");
      return mlir::failure();
    }

    Equations<MatchedEquation> matchedEquations;

    for (auto& solution : matchingGraph.getMatch()) {
      auto clone = solution.getEquation()->clone();

      matchedEquations.add(std::make_unique<MatchedEquation>(
          std::move(clone), solution.getIndexes(), solution.getAccess()));
    }

    result.setVariables(model.getVariables());
    result.setEquations(matchedEquations);

    return mlir::success();
  }
}
