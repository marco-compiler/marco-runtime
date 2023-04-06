#include "marco/Codegen/Lowering/StandardFunctionLowerer.h"

using namespace ::marco;
using namespace ::marco::codegen;
using namespace ::mlir::modelica;

namespace marco::codegen::lowering
{
  StandardFunctionLowerer::StandardFunctionLowerer(BridgeInterface* bridge)
      : Lowerer(bridge)
  {
  }

  void StandardFunctionLowerer::declare(const ast::StandardFunction& function)
  {
    mlir::Location location = loc(function.getLocation());

    // Create the record operation.
    auto functionOp = builder().create<FunctionOp>(location, function.getName());

    mlir::OpBuilder::InsertionGuard guard(builder());
    builder().setInsertionPointToStart(functionOp.bodyBlock());

    // Declare the inner classes.
    for (const auto& innerClassNode : function.getInnerClasses()) {
      declare(*innerClassNode->cast<ast::Class>());
    }
  }

  void StandardFunctionLowerer::declareVariables(
      const ast::StandardFunction& function)
  {
    mlir::OpBuilder::InsertionGuard guard(builder());
    LookupScopeGuard lookupScopeGuard(&getContext());

    // Get the operation.
    auto functionOp = mlir::cast<FunctionOp>(getClass(function));
    pushLookupScope(functionOp);
    builder().setInsertionPointToEnd(functionOp.bodyBlock());

    // Declare the variables.
    for (const auto& variable : function.getVariables()) {
      declare(*variable->cast<ast::Member>());
    }

    // Declare the variables of inner classes.
    for (const auto& innerClassNode : function.getInnerClasses()) {
      declareVariables(*innerClassNode->cast<ast::Class>());
    }
  }

  void StandardFunctionLowerer::lower(const ast::StandardFunction& function)
  {
    mlir::OpBuilder::InsertionGuard guard(builder());

    Lowerer::VariablesScope varScope(getVariablesSymbolTable());
    LookupScopeGuard lookupScopeGuard(&getContext());

    // Get the operation.
    auto functionOp = mlir::cast<FunctionOp>(getClass(function));
    pushLookupScope(functionOp);
    builder().setInsertionPointToEnd(functionOp.bodyBlock());

    // Map the variables.
    insertVariable(
        "time",
        Reference::time(builder(), builder().getUnknownLoc()));

    for (VariableOp variableOp : functionOp.getVariables()) {
      insertVariable(
          variableOp.getSymName(),
          Reference::variable(
              builder(), variableOp->getLoc(),
              variableOp.getSymName(),
              variableOp.getVariableType().unwrap()));
    }

    // Lower the annotations.
    llvm::SmallVector<llvm::StringRef, 3> inputVarNames;

    for (VariableOp variable : functionOp.getVariables()) {
      if (variable.isInput()) {
        inputVarNames.emplace_back(variable.getName());
      }
    }

    llvm::SmallVector<llvm::StringRef, 1> outputVarNames;

    for (VariableOp variable : functionOp.getVariables()) {
      if (variable.isOutput()) {
        outputVarNames.emplace_back(variable.getName());
      }
    }

    if (function.hasAnnotation()) {
      const auto* annotation = function.getAnnotation();

      // Inline attribute.
      functionOp->setAttr(
          "inline",
          builder().getBoolAttr(
              function.getAnnotation()->getInlineProperty()));

      // Inverse functions attribute.
      auto inverseFunctionAnnotation =
          annotation->getInverseFunctionAnnotation();

      InverseFunctionsMap map;

      // Create a map of the function members indexes for faster retrieval.
      llvm::StringMap<unsigned int> indexes;

      for (const auto& name : llvm::enumerate(inputVarNames)) {
        indexes[name.value()] = name.index();
      }

      for (const auto& name : llvm::enumerate(outputVarNames)) {
        indexes[name.value()] = inputVarNames.size() + name.index();
      }

      mlir::StorageUniquer::StorageAllocator allocator;

      // Iterate over the input arguments and for each invertible one
      // add the function to the inverse map.
      for (const auto& arg : inputVarNames) {
        if (!inverseFunctionAnnotation.isInvertible(arg)) {
          continue;
        }

        auto inverseArgs = inverseFunctionAnnotation.getInverseArgs(arg);
        llvm::SmallVector<unsigned int, 3> permutation;

        for (const auto& inverseArg : inverseArgs) {
          assert(indexes.find(inverseArg) != indexes.end());
          permutation.push_back(indexes[inverseArg]);
        }

        map[indexes[arg]] = std::make_pair(
            inverseFunctionAnnotation.getInverseFunction(arg),
            allocator.copyInto(llvm::ArrayRef<unsigned int>(permutation)));
      }

      if (!map.empty()) {
        auto inverseFunctionAttribute =
            InverseFunctionsAttr::get(builder().getContext(), map);

        functionOp->setAttr("inverse", inverseFunctionAttribute);
      }

      if (annotation->hasDerivativeAnnotation()) {
        auto derivativeAnnotation = annotation->getDerivativeAnnotation();

        auto derivativeAttribute = DerivativeAttr::get(
            builder().getContext(),
            derivativeAnnotation.getName(),
            derivativeAnnotation.getOrder());

        functionOp->setAttr("derivative", derivativeAttribute);
      }
    }

    // Create the default values for variables.
    for (const auto& variable : function.getVariables()) {
      lowerVariableDefaultValue(*variable->cast<ast::Member>());
    }

    // Lower the body.
    lowerClassBody(function);

    // Lower the inner classes.
    for (const auto& innerClassNode : function.getInnerClasses()) {
      lower(*innerClassNode->cast<ast::Class>());
    }
  }

  void StandardFunctionLowerer::lowerVariableDefaultValue(
      const ast::Member& variable)
  {
    if (!variable.hasExpression()) {
      return;
    }

    const ast::Expression* expression = variable.getExpression();

    mlir::Location expressionLoc = loc(expression->getLocation());

    auto defaultOp = builder().create<DefaultOp>(
        expressionLoc, variable.getName());

    mlir::OpBuilder::InsertionGuard guard(builder());
    mlir::Block* bodyBlock = builder().createBlock(&defaultOp.getBodyRegion());
    builder().setInsertionPointToStart(bodyBlock);

    mlir::Value value = lower(*expression)[0].get(expressionLoc);
    builder().create<YieldOp>(expressionLoc, value);
  }
}
