#include "marco/Codegen/Transforms/Model/ModelConverter.h"
#include "marco/Codegen/Transforms/Model/IDA.h"
#include "marco/Dialect/IDA/IDADialect.h"
#include "mlir/Dialect/LLVMIR/FunctionCallUtils.h"
#include "mlir/Dialect/SCF/SCF.h"
#include "mlir/IR/AffineMap.h"

using namespace ::marco;
using namespace ::marco::codegen;
using namespace ::marco::codegen::modelica;
using namespace ::marco::modeling;

namespace
{
  /// Class to be used to uniquely identify an equation template function.
  /// Two templates are considered to be equal if they refer to the same EquationOp and have
  /// the same scheduling direction, which impacts on the function body itself due to the way
  /// the iteration indexes are updated.
  class EquationTemplateInfo
  {
    public:
      EquationTemplateInfo(EquationOp equation, modeling::scheduling::Direction schedulingDirection)
          : equation(equation.getOperation()), schedulingDirection(schedulingDirection)
      {
      }

      bool operator<(const EquationTemplateInfo& other) const {
        return equation < other.equation && schedulingDirection < other.schedulingDirection;
      }

    private:
      mlir::Operation* equation;
      modeling::scheduling::Direction schedulingDirection;
  };
}

/// Get the MLIR function with the given name, or declare it inside the module if not present.
static mlir::LLVM::LLVMFuncOp getOrInsertFunction(
    mlir::OpBuilder& builder,
    mlir::ModuleOp module,
    llvm::StringRef name,
    mlir::LLVM::LLVMFunctionType type)
{
  if (auto foo = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(name)) {
    return foo;
  }

  mlir::OpBuilder::InsertionGuard guard(builder);
  builder.setInsertionPointToStart(module.getBody());
  return builder.create<mlir::LLVM::LLVMFuncOp>(module.getLoc(), name, type);
}

/// Remove the unused arguments of a function and also update the function calls
/// to reflect the function signature change.
template<typename CallIt>
static std::vector<unsigned int> removeUnusedArguments(
    mlir::FuncOp function, CallIt callsBegin, CallIt callsEnd)
{
  std::vector<unsigned int> removedArgs;

  // Determine the unused arguments
  for (const auto& arg : function.getArguments()) {
    if (arg.getUsers().empty()) {
      removedArgs.push_back(arg.getArgNumber());
    }
  }

  if (!removedArgs.empty()) {
    // Erase the unused function arguments
    function.eraseArguments(removedArgs);

    // Update the function calls
    std::sort(removedArgs.begin(), removedArgs.end());

    for (auto callIt = callsBegin; callIt != callsEnd; ++callIt) {
      mlir::CallOp call = callIt->second;

      for (auto argIt = removedArgs.rbegin(); argIt != removedArgs.rend(); ++argIt) {
        call->eraseOperand(*argIt);
      }
    }
  }

  return removedArgs;
}

namespace marco::codegen
{
  ModelConverter::ModelConverter(SolveModelOptions options, mlir::LLVMTypeConverter& typeConverter)
      : options(std::move(options)),
        typeConverter(&typeConverter)
  {
  }

  mlir::LogicalResult ModelConverter::convert(
      mlir::OpBuilder& builder,
      const Model<ScheduledEquationsBlock>& model,
      const mlir::BlockAndValueMapping& derivatives) const
  {
    ModelOp modelOp = model.getOperation();

    // Convert the original derivatives map between values into a map between positions
    DerivativesPositionsMap derivativesPositions;

    for (size_t i = 0, e = modelOp.body().getNumArguments(); i < e; ++i) {
      mlir::Value var = modelOp.body().getArgument(i);

      if (derivatives.contains(var)) {
        mlir::Value derivative = derivatives.lookup(var);
        bool derivativeFound = false;
        unsigned int position = 0;

        for (size_t j = 0; j < e && !derivativeFound; ++j) {
          mlir::Value arg = modelOp.body().getArgument(j);

          if (arg == derivative) {
            derivativeFound = true;
            position = j;
          }
        }

        assert(derivativeFound && "Derivative not found among arguments");
        derivativesPositions[i] = position;
      }
    }

    // Create the external solvers
    ExternalSolvers solvers;

    ConversionInfo conversionInfo;

    // Determine which equations can be potentially processed by MARCO.
    // Those are the ones that can me bade explicit with respect to the matched variable and
    // the non-cyclic ones.
    for (auto& scheduledBlock : model.getScheduledBlocks()) {
      if (!scheduledBlock->hasCycle()) {
        for (auto& scheduledEquation : *scheduledBlock) {
          auto explicitClone = scheduledEquation->cloneIRAndExplicitate(builder);

          if (explicitClone == nullptr) {
            conversionInfo.implicitEquations.emplace(scheduledEquation.get());
          } else {
            auto& movedClone = *conversionInfo.explicitEquations.emplace(std::move(explicitClone)).first;
            conversionInfo.explicitEquationsMap[scheduledEquation.get()] = movedClone.get();
          }
        }
      } else {
        for (const auto& equation : *scheduledBlock) {
          conversionInfo.cyclicEquations.emplace(equation.get());
        }
      }
    }

    // Add the implicit equations to the set of equations managed by IDA, together with their
    // written variables.
    for (const auto& implicitEquation : conversionInfo.implicitEquations) {
      auto var = implicitEquation->getWrite().getVariable();
      solvers.ida->addVariable(var->getValue());
      solvers.ida->addEquation(implicitEquation);
    }

    // Add the cyclic equations to the set of equations managed by IDA, together with their
    // written variables.
    for (const auto& cyclicEquation : conversionInfo.cyclicEquations) {
      auto var = cyclicEquation->getWrite().getVariable();
      solvers.ida->addVariable(var->getValue());
      solvers.ida->addEquation(cyclicEquation);
    }

    for (const auto& scheduledBlock : model.getScheduledBlocks()) {
      for (auto& scheduledEquation : *scheduledBlock) {
        auto var = scheduledEquation->getWrite().getVariable();

        if (derivatives.contains(var->getValue())) {
          solvers.ida->addVariable(var->getValue());
          solvers.ida->addEquation(scheduledEquation.get());
        }
      }
    }

    // If any of the remaining equations manageable by MARCO does write on a variable managed
    // by IDA, then the equation must be passed to IDA even if not strictly necessary.
    // Avoiding this would require either memory duplication or a more severe restructuring
    // of the solving infrastructure, which would have to be able to split variables and equations
    // according to which runtime solver manages such variables.
    for (const auto& scheduledBlock : model.getScheduledBlocks()) {
      for (auto& scheduledEquation : *scheduledBlock) {
        auto var = scheduledEquation->getWrite().getVariable();

        if (solvers.ida->hasVariable(var->getValue())) {
          solvers.ida->addEquation(scheduledEquation.get());
        }
      }
    }

    // Create the various functions composing the simulation
    if (auto res = createInitFunction(builder, model, solvers, derivatives); failed(res)) {
      model.getOperation().emitError("Could not create the '" + initFunctionName + "' function");
      return res;
    }

    if (auto res = createDeinitFunction(builder, modelOp); failed(res)) {
      model.getOperation().emitError("Could not create the '" + deinitFunctionName + "' function");
      return res;
    }

    if (auto res = createUpdateNonStateVariablesFunction(builder, model, conversionInfo, solvers); mlir::failed(res)) {
      model.getOperation().emitError("Could not create the '" + updateNonStateVariablesFunctionName + "' function");
      return res;
    }

    if (auto res = createUpdateStateVariablesFunction(builder, modelOp, derivativesPositions); mlir::failed(res)) {
      model.getOperation().emitError("Could not create the '" + updateStateVariablesFunctionName + "' function");
      return res;
    }

    if (auto res = createIncrementTimeFunction(builder, model); mlir::failed(res)) {
      model.getOperation().emitError("Could not create the '" + incrementTimeFunctionName + "' function");
      return res;
    }

    if (auto res = createPrintHeaderFunction(builder, modelOp, derivativesPositions); failed(res)) {
      model.getOperation().emitError("Could not create the '" + printHeaderFunctionName + "' function");
      return res;
    }

    if (auto res = createPrintFunction(builder, modelOp, derivativesPositions); failed(res)) {
      model.getOperation().emitError("Could not create the '" + printFunctionName + "' function");
      return res;
    }

    if (options.emitMain) {
      if (auto res = createMainFunction(builder, model); mlir::failed(res)) {
        model.getOperation().emitError("Could not create the '" + mainFunctionName + "' function");
        return res;
      }
    }

    // Erase the model operation, which has been converted to algorithmic code
    model.getOperation().erase();

    return mlir::success();
  }

  mlir::Type ModelConverter::getVoidPtrType() const
  {
    return mlir::LLVM::LLVMPointerType::get(mlir::IntegerType::get(&typeConverter->getContext(), 8));
  }

  mlir::LLVM::LLVMFuncOp ModelConverter::lookupOrCreateHeapAllocFn(
      mlir::OpBuilder& builder, mlir::ModuleOp module) const
  {
    std::string name = "_MheapAlloc_pvoid_i64";

    if (auto foo = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(name)) {
      return foo;
    }

    mlir::PatternRewriter::InsertionGuard insertGuard(builder);
    builder.setInsertionPointToStart(module.getBody());
    auto llvmFnType = mlir::LLVM::LLVMFunctionType::get(getVoidPtrType(), builder.getI64Type());
    return builder.create<mlir::LLVM::LLVMFuncOp>(module->getLoc(), name, llvmFnType);
  }

  mlir::LLVM::LLVMFuncOp ModelConverter::lookupOrCreateHeapFreeFn(
      mlir::OpBuilder& builder, mlir::ModuleOp module) const
  {
    std::string name = "_MheapFree_void_pvoid";

    if (auto foo = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(name)) {
      return foo;
    }

    mlir::PatternRewriter::InsertionGuard guard(builder);
    builder.setInsertionPointToStart(module.getBody());
    mlir::Type voidType = mlir::LLVM::LLVMVoidType::get(module.getContext());
    auto llvmFnType = mlir::LLVM::LLVMFunctionType::get(voidType, getVoidPtrType());
    return builder.create<mlir::LLVM::LLVMFuncOp>(module->getLoc(), name, llvmFnType);
  }

  mlir::LogicalResult ModelConverter::createMainFunction(
      mlir::OpBuilder& builder, const Model<ScheduledEquationsBlock>& model) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    ModelOp modelOp = model.getOperation();
    mlir::Location loc = modelOp.getLoc();

    // Create the function inside the parent module
    auto module = modelOp->getParentOfType<mlir::ModuleOp>();
    builder.setInsertionPointToEnd(module.getBody());

    llvm::SmallVector<mlir::Type, 3> argsTypes;
    llvm::SmallVector<mlir::Type, 3> resultsTypes;

    argsTypes.push_back(builder.getI32Type());
    argsTypes.push_back(mlir::LLVM::LLVMPointerType::get(mlir::LLVM::LLVMPointerType::get(builder.getIntegerType(8))));
    resultsTypes.push_back(builder.getI32Type());

    auto function = builder.create<mlir::FuncOp>(
        loc, mainFunctionName, builder.getFunctionType(argsTypes, resultsTypes));

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    // Call the function to start the simulation.
    // Its definition lives within the runtime library.
    mlir::Type voidType = mlir::LLVM::LLVMVoidType::get(modelOp.getContext());

    auto runFunction = getOrInsertFunction(
        builder, module, runFunctionName, mlir::LLVM::LLVMFunctionType::get(voidType, llvm::None));

    builder.create<mlir::LLVM::CallOp>(loc, runFunction, llvm::None);

    // Create the return statement
    mlir::Value returnValue = builder.create<mlir::ConstantOp>(loc, builder.getI32IntegerAttr(0));
    builder.create<mlir::ReturnOp>(loc, returnValue);

    return mlir::success();
  }

  mlir::Value ModelConverter::loadDataFromOpaquePtr(
      mlir::OpBuilder& builder,
      mlir::Value ptr,
      mlir::TypeRange varTypes) const
  {
    mlir::Location loc = ptr.getLoc();
    llvm::SmallVector<mlir::Type, 3> structTypes;

    for (const auto& type : varTypes) {
      structTypes.push_back(typeConverter->convertType(type));
    }

    mlir::Type structType = mlir::LLVM::LLVMStructType::getLiteral(ptr.getContext(), structTypes);
    mlir::Type structPtrType = mlir::LLVM::LLVMPointerType::get(structType);
    mlir::Value structPtr = builder.create<mlir::LLVM::BitcastOp>(loc, structPtrType, ptr);
    mlir::Value structValue = builder.create<mlir::LLVM::LoadOp>(loc, structPtr);

    return structValue;
  }

  mlir::Value ModelConverter::extractValue(
      mlir::OpBuilder& builder,
      mlir::Value structValue,
      mlir::Type type,
      unsigned int position) const
  {
    mlir::Location loc = structValue.getLoc();

    assert(structValue.getType().isa<mlir::LLVM::LLVMStructType>() && "Not an LLVM struct");
    auto structType = structValue.getType().cast<mlir::LLVM::LLVMStructType>();
    auto structTypes = structType.getBody();
    assert (position < structTypes.size() && "LLVM struct: index is out of bounds");

    mlir::Value var = builder.create<mlir::LLVM::ExtractValueOp>(loc, structTypes[position], structValue, builder.getIndexArrayAttr(position));
    return typeConverter->materializeSourceConversion(builder, loc, type, var);
  }

  mlir::Value ModelConverter::convertMember(mlir::OpBuilder& builder, MemberCreateOp op) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);

    using LoadReplacer = std::function<void(MemberLoadOp)>;
    using StoreReplacer = std::function<void(MemberStoreOp)>;

    mlir::Location loc = op->getLoc();

    auto memberType = op.resultType().cast<MemberType>();
    auto arrayType = memberType.toArrayType();
    assert(arrayType.getAllocationScope() == BufferAllocationScope::heap);

    // Create the memory buffer for the variable
    builder.setInsertionPoint(op);

    mlir::Value reference = builder.create<AllocOp>(
        loc, arrayType.getElementType(), arrayType.getShape(), op.dynamicDimensions(), false);

    // Replace loads and stores with appropriate instructions operating on the new memory buffer.
    // The way such replacements are executed depend on the nature of the variable.

    auto replacers = [&]() {
      if (arrayType.isScalar()) {
        assert(op.dynamicDimensions().empty());

        auto loadReplacer = [&builder, reference](MemberLoadOp loadOp) -> void {
          mlir::OpBuilder::InsertionGuard guard(builder);
          builder.setInsertionPoint(loadOp);
          loadOp.replaceAllUsesWith(reference);
          loadOp.erase();
        };

        auto storeReplacer = [&builder, reference](MemberStoreOp storeOp) -> void {
          mlir::OpBuilder::InsertionGuard guard(builder);
          builder.setInsertionPoint(storeOp);
          auto assignment = builder.create<AssignmentOp>(storeOp.getLoc(), storeOp.value(), reference);
          storeOp->replaceAllUsesWith(assignment);
          storeOp.erase();
        };

        return std::make_pair<LoadReplacer, StoreReplacer>(loadReplacer, storeReplacer);
      }

      auto loadReplacer = [&builder, reference](MemberLoadOp loadOp) -> void {
        mlir::OpBuilder::InsertionGuard guard(builder);
        builder.setInsertionPoint(loadOp);
        loadOp.replaceAllUsesWith(reference);
        loadOp.erase();
      };

      auto storeReplacer = [&builder, reference](MemberStoreOp storeOp) -> void {
        mlir::OpBuilder::InsertionGuard guard(builder);
        builder.setInsertionPoint(storeOp);
        auto assignment = builder.create<AssignmentOp>(storeOp.getLoc(), storeOp.value(), reference);
        storeOp->replaceAllUsesWith(assignment);
        storeOp.erase();
      };

      return std::make_pair<LoadReplacer, StoreReplacer>(loadReplacer, storeReplacer);
    };

    LoadReplacer loadReplacer;
    StoreReplacer storeReplacer;
    std::tie(loadReplacer, storeReplacer) = replacers();

    for (auto* user : op->getUsers()) {
      if (auto loadOp = mlir::dyn_cast<MemberLoadOp>(user)) {
        loadReplacer(loadOp);
      } else if (auto storeOp = mlir::dyn_cast<MemberStoreOp>(user)) {
        storeReplacer(storeOp);
      }
    }

    op.replaceAllUsesWith(reference);
    op.erase();
    return reference;
  }

  mlir::LogicalResult ModelConverter::createInitFunction(
      mlir::OpBuilder& builder,
      const Model<ScheduledEquationsBlock>& model,
      ExternalSolvers& externalSolvers,
      const mlir::BlockAndValueMapping& derivatives) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    auto modelOp = model.getOperation();
    mlir::Location loc = modelOp.getLoc();
    auto module = modelOp->getParentOfType<mlir::ModuleOp>();

    // Create the function inside the parent module
    builder.setInsertionPointToEnd(module.getBody());

    auto functionType = builder.getFunctionType(llvm::None, getVoidPtrType());
    auto function = builder.create<mlir::FuncOp>(loc, initFunctionName, functionType);

    mlir::BlockAndValueMapping mapping;

    // Move the initialization instructions into the new function
    modelOp.init().cloneInto(&function.getBody(), mapping);
    builder.setInsertionPointToStart(&function.getBody().front());

    // Convert the 'time' from a member type to a scalar value and set the start time
    auto terminator = mlir::cast<YieldOp>(function.getBody().back().getTerminator());
    auto timeMember = terminator.values()[0].getDefiningOp<MemberCreateOp>();
    auto startTime = builder.create<ConstantOp>(loc, modelOp.startTime());
    timeMember->replaceAllUsesWith(startTime);
    timeMember.erase();

    // Initialize the IDA solver
    if (auto res = externalSolvers.ida->init(builder, function); mlir::failed(res)) {
      return res;
    }

    // Process the variables that are handled by IDA
    if (auto res = externalSolvers.ida->processVariables(builder, function, derivatives); mlir::failed(res)) {
      return res;
    }

    // Process the equations that are handled by IDA
    if (auto res = externalSolvers.ida->processEquations(builder, model, function, modelOp.body().getArgumentTypes(), derivatives); mlir::failed(res)) {
      return res;
    }

    // The values to be packed into the structure to be passed around the simulation functions
    llvm::SmallVector<mlir::Value, 3> values;

    builder.setInsertionPointAfter(terminator);

    auto removeAllocationScopeFn = [&](mlir::Value value) -> mlir::Value {
      mlir::Type type = value.getType();

      if (auto arrayType = type.dyn_cast<ArrayType>()) {
        return builder.create<ArrayCastOp>(
            loc, value,
            value.getType().cast<ArrayType>().toUnknownAllocationScope());
      } else {
        return value;
      }
    };

    // Add variables to the struct to be passed around (i.e. to the step and
    // print functions).

    for (const auto& var : terminator.values()) {
      auto memberCreateOp = var.getDefiningOp<MemberCreateOp>();

      if (memberCreateOp != nullptr) {
        mlir::Value array = convertMember(builder, memberCreateOp);
        builder.setInsertionPointAfterValue(array);
        values.push_back(removeAllocationScopeFn(array));
      } else {
        builder.setInsertionPointAfterValue(var);
        values.push_back(removeAllocationScopeFn(var));
      }
    }

    builder.setInsertionPointAfter(terminator);

    // Determine the types composing the runtime data structure.
    std::vector<mlir::Type> runtimeDataStructTypes;

    // Determine the types of the external solvers structure.
    std::vector<mlir::Type> externalSolversStructTypes;

    externalSolversStructTypes.push_back(typeConverter->convertType(externalSolvers.ida->getSolverInstanceType(builder.getContext())));
    mlir::Type externalSolversStructType = mlir::LLVM::LLVMStructType::getLiteral(builder.getContext(), externalSolversStructTypes);
    mlir::Type externalSolversStructPtrType = mlir::LLVM::LLVMPointerType::get(externalSolversStructType);

    // Add the pointer to the external solvers structure
    runtimeDataStructTypes.push_back(externalSolversStructPtrType);

    // Add the types of the variables.
    // Notice that the 'time' variable is already present among the model
    // body, so its type is automatically added to the struct.
    for (const auto& type : modelOp.body().getArgumentTypes()) {
      auto convertedType = typeConverter->convertType(type);
      runtimeDataStructTypes.push_back(convertedType);
    }

    auto runtimeDataStructType = mlir::LLVM::LLVMStructType::getLiteral(modelOp.getContext(), runtimeDataStructTypes);

    mlir::Value structValue = builder.create<mlir::LLVM::UndefOp>(loc, runtimeDataStructType);

    for (const auto& var : llvm::enumerate(values)) {
      mlir::Type convertedType = typeConverter->convertType(var.value().getType());
      mlir::Value convertedVar = typeConverter->materializeTargetConversion(builder, loc, convertedType, var.value());
      structValue = builder.create<mlir::LLVM::InsertValueOp>(loc, structValue, convertedVar, builder.getIndexArrayAttr(var.index()));
    }

    // The data structure must be stored on the heap in order to escape
    // from the function.

    // Add the "malloc" function to the module
    auto heapAllocFunc = lookupOrCreateHeapAllocFn(builder, module);

    // Determine the size (in bytes) of the memory to be allocated
    mlir::Type structPtrType = mlir::LLVM::LLVMPointerType::get(runtimeDataStructType);
    mlir::Value nullPtr = builder.create<mlir::LLVM::NullOp>(loc, structPtrType);

    mlir::Value one = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(1));
    one = typeConverter->materializeTargetConversion(builder, loc, typeConverter->getIndexType(), one);

    mlir::Value gepPtr = builder.create<mlir::LLVM::GEPOp>(loc, structPtrType, llvm::ArrayRef<mlir::Value>{nullPtr, one});
    mlir::Value sizeBytes = builder.create<mlir::LLVM::PtrToIntOp>(loc, typeConverter->getIndexType(), gepPtr);
    mlir::Value resultOpaquePtr = createLLVMCall(builder, loc, heapAllocFunc, sizeBytes, getVoidPtrType())[0];

    // Store the struct into the heap memory
    mlir::Value resultCastedPtr = builder.create<mlir::LLVM::BitcastOp>(loc, structPtrType, resultOpaquePtr);
    builder.create<mlir::LLVM::StoreOp>(loc, structValue, resultCastedPtr);

    builder.create<mlir::ReturnOp>(loc, resultOpaquePtr);
    terminator->erase();

    return mlir::success();
  }

  mlir::LogicalResult ModelConverter::createDeinitFunction(mlir::OpBuilder& builder, ModelOp modelOp) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    mlir::Location loc = modelOp.getLoc();
    auto module = modelOp->getParentOfType<mlir::ModuleOp>();

    // Create the function inside the parent module
    builder.setInsertionPointToEnd(module.getBody());

    auto function = builder.create<mlir::FuncOp>(
        loc, deinitFunctionName,
        builder.getFunctionType(getVoidPtrType(), llvm::None));

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    // Extract the data from the struct
    mlir::TypeRange varTypes = modelOp.body().getArgumentTypes();
    mlir::Value structValue = loadDataFromOpaquePtr(builder, function.getArgument(0), varTypes);

    // Deallocate the arrays
    for (const auto& type : llvm::enumerate(varTypes)) {
      if (auto arrayType = type.value().dyn_cast<ArrayType>()) {
        mlir::Value var = extractValue(builder, structValue, varTypes[type.index()], type.index());
        var = builder.create<ArrayCastOp>(loc, var, arrayType.toAllocationScope(BufferAllocationScope::heap));
        builder.create<FreeOp>(loc, var);
      }
    }

    // Add "free" function to the module
    auto freeFunc = lookupOrCreateHeapFreeFn(builder, module);

    // Deallocate the data structure
    builder.create<mlir::LLVM::CallOp>(
        loc, llvm::None, builder.getSymbolRefAttr(freeFunc), function.getArgument(0));

    builder.create<mlir::ReturnOp>(loc);
    return mlir::success();
  }

  mlir::FuncOp ModelConverter::createEquationFunction(
      mlir::OpBuilder& builder,
      const ScheduledEquation& equation,
      llvm::StringRef equationFunctionName,
      mlir::FuncOp templateFunction,
      std::multimap<mlir::FuncOp, mlir::CallOp>& equationTemplateCalls,
      mlir::TypeRange varsTypes) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    mlir::Location loc = equation.getOperation().getLoc();

    auto module = equation.getOperation()->getParentOfType<mlir::ModuleOp>();
    builder.setInsertionPointToEnd(module.getBody());

    auto functionType = builder.getFunctionType(varsTypes, llvm::None);
    auto function = builder.create<mlir::FuncOp>(loc, equationFunctionName, functionType);

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    auto valuesFn = [&](marco::modeling::scheduling::Direction iterationDirection, Range range) -> std::tuple<mlir::Value, mlir::Value, mlir::Value> {
      assert(iterationDirection == marco::modeling::scheduling::Direction::Forward ||
          iterationDirection == marco::modeling::scheduling::Direction::Backward);

      if (iterationDirection == marco::modeling::scheduling::Direction::Forward) {
        mlir::Value begin = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.getBegin()));
        mlir::Value end = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.getEnd()));
        mlir::Value step = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(1));

        return std::make_tuple(begin, end, step);
      }

      mlir::Value begin = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.getEnd() - 1));
      mlir::Value end = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.getBegin() - 1));
      mlir::Value step = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(1));

      return std::make_tuple(begin, end, step);
    };

    std::vector<mlir::Value> args;
    auto iterationRanges = equation.getIterationRanges();

    for (size_t i = 0, e = equation.getNumOfIterationVars(); i < e; ++i) {
      auto values = valuesFn(equation.getSchedulingDirection(), iterationRanges[i]);

      args.push_back(std::get<0>(values));
      args.push_back(std::get<1>(values));
      args.push_back(std::get<2>(values));
    }

    mlir::ValueRange vars = function.getArguments();
    args.insert(args.end(), vars.begin(), vars.end());

    // Call the equation template function
    auto templateFunctionCall = builder.create<mlir::CallOp>(loc, templateFunction, args);
    equationTemplateCalls.emplace(templateFunction, templateFunctionCall);

    builder.create<mlir::ReturnOp>(loc);
    return function;
  }

  mlir::LogicalResult ModelConverter::createUpdateNonStateVariablesFunction(
      mlir::OpBuilder& builder,
      const Model<ScheduledEquationsBlock>& model,
      const ConversionInfo& conversionInfo,
      ExternalSolvers& externalSolvers) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    ModelOp modelOp = model.getOperation();
    mlir::Location loc = modelOp.getLoc();

    // Create the function inside the parent module
    builder.setInsertionPointToEnd(modelOp->getParentOfType<mlir::ModuleOp>().getBody());

    auto function = builder.create<mlir::FuncOp>(
        loc, updateNonStateVariablesFunctionName,
        builder.getFunctionType(getVoidPtrType(), llvm::None));

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    // Extract the data from the struct
    mlir::TypeRange varTypes = modelOp.body().getArgumentTypes();
    mlir::Value structValue = loadDataFromOpaquePtr(builder, function.getArgument(0), varTypes);

    llvm::SmallVector<mlir::Value, 3> vars;

    for (const auto& varType : llvm::enumerate(varTypes)) {
      vars.push_back(extractValue(builder, structValue, varTypes[varType.index()], varType.index()));
    }

    // Convert the equations into algorithmic code
    size_t equationTemplateCounter = 0;
    size_t equationCounter = 0;

    std::map<EquationTemplateInfo, mlir::FuncOp> equationTemplatesMap;
    std::set<mlir::FuncOp> equationTemplateFunctions;
    std::multimap<mlir::FuncOp, mlir::CallOp> equationTemplateCalls;
    std::set<mlir::FuncOp> equationFunctions;
    std::multimap<mlir::FuncOp, mlir::CallOp> equationCalls;

    // Get or create the template equation function for a scheduled equation
    auto getEquationTemplateFn = [&](const ScheduledEquation* equation) -> mlir::FuncOp {
      EquationTemplateInfo requestedTemplate(equation->getOperation(), equation->getSchedulingDirection());
      auto it = equationTemplatesMap.find(requestedTemplate);

      if (it != equationTemplatesMap.end()) {
        return it->second;
      }

      std::string templateFunctionName = "eq_template_" + std::to_string(equationTemplateCounter);
      ++equationTemplateCounter;

      auto explicitEquation = llvm::find_if(conversionInfo.explicitEquationsMap, [&](const auto& equationPtr) {
        return equationPtr.first == equation;
      });

      assert(explicitEquation != conversionInfo.explicitEquationsMap.end());

      // Create the equation template function
      auto templateFunction = explicitEquation->second->createTemplateFunction(
          builder, templateFunctionName, modelOp.body().getArguments(), equation->getSchedulingDirection());

      return templateFunction;
    };

    for (const auto& scheduledBlock : model.getScheduledBlocks()) {
      for (const auto& equation : *scheduledBlock) {
        if (externalSolvers.ida->hasEquation(equation.get())) {
          // Let IDA process the equation
          // TODO
          return mlir::failure();

        } else {
          // The equation is handled by MARCO
          auto templateFunction = getEquationTemplateFn(equation.get());

          if (templateFunction == nullptr) {
            return mlir::failure();
          }

          equationTemplateFunctions.insert(templateFunction);

          // Create the function that calls the template.
          // This function dictates the indices the template will work with.
          std::string equationFunctionName = "eq_" + std::to_string(equationCounter);
          ++equationCounter;

          auto equationFunction = createEquationFunction(
              builder, *equation, equationFunctionName, templateFunction,
              equationTemplateCalls,
              modelOp.body().getArgumentTypes());

          equationFunctions.insert(equationFunction);

          // Create the call to the instantiated template function
          auto equationCall = builder.create<mlir::CallOp>(loc, equationFunction, vars);
          equationCalls.emplace(equationFunction, equationCall);
        }
      }
    }

    builder.create<mlir::ReturnOp>(loc);

    // Remove the unused function arguments
    for (const auto& equationTemplateFunction : equationTemplateFunctions) {
      auto calls = equationTemplateCalls.equal_range(equationTemplateFunction);
      removeUnusedArguments(equationTemplateFunction, calls.first, calls.second);
    }

    for (const auto& equationFunction : equationFunctions) {
      auto calls = equationCalls.equal_range(equationFunction);
      removeUnusedArguments(equationFunction, calls.first, calls.second);
    }

    return mlir::success();
  }

  mlir::LogicalResult ModelConverter::createUpdateStateVariablesFunction(
      mlir::OpBuilder& builder, ModelOp modelOp, const DerivativesPositionsMap& derivatives) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    mlir::Location loc = modelOp.getLoc();

    // Create the function inside the parent module
    builder.setInsertionPointToEnd(modelOp->getParentOfType<mlir::ModuleOp>().getBody());

    auto function = builder.create<mlir::FuncOp>(
        loc, updateStateVariablesFunctionName,
        builder.getFunctionType(getVoidPtrType(), llvm::None));

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    // Extract the state variables from the opaque pointer
    mlir::TypeRange varTypes = modelOp.body().getArgumentTypes();
    mlir::Value structValue = loadDataFromOpaquePtr(builder, function.getArgument(0), varTypes);

    // Update the state variables by applying the forward Euler method
    mlir::Value timeStep = builder.create<ConstantOp>(loc, modelOp.timeStep());

    std::vector<std::pair<mlir::Value, mlir::Value>> varsAndDers;

    for (const auto& variable : modelOp.body().getArguments()) {
      size_t index = variable.getArgNumber();
      auto it = derivatives.find(variable.getArgNumber());

      if (it != derivatives.end()) {
        mlir::Value var = extractValue(builder, structValue, varTypes[index], index);
        mlir::Value der = extractValue(builder, structValue, varTypes[it->second], it->second);
        varsAndDers.emplace_back(var, der);
      }
    }

    for (const auto& [var, der] : varsAndDers) {
      mlir::Value nextValue = builder.create<MulOp>(loc, der.getType(), der, timeStep);
      nextValue = builder.create<AddOp>(loc, var.getType(), nextValue, var);
      builder.create<AssignmentOp>(loc, nextValue, var);
    }

    builder.create<mlir::ReturnOp>(loc);
    return mlir::success();
  }

  mlir::LogicalResult ModelConverter::createIncrementTimeFunction(
      mlir::OpBuilder& builder,
      const Model<ScheduledEquationsBlock>& model) const
  {
    mlir::OpBuilder::InsertionGuard guard(builder);
    ModelOp modelOp = model.getOperation();
    mlir::Location loc = modelOp.getLoc();

    // Create the function inside the parent module
    builder.setInsertionPointToEnd(modelOp->getParentOfType<mlir::ModuleOp>().getBody());

    auto function = builder.create<mlir::FuncOp>(
        loc, incrementTimeFunctionName,
        builder.getFunctionType(getVoidPtrType(), builder.getI1Type()));

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    // Extract the data from the struct
    mlir::TypeRange varTypes = modelOp.body().getArgumentTypes();
    mlir::Value structValue = loadDataFromOpaquePtr(builder, function.getArgument(0), varTypes);

    mlir::Value time = extractValue(builder, structValue, varTypes[timeVariablePosition], 0);

    // Increment the time
    mlir::Value timeStep = builder.create<ConstantOp>(loc, modelOp.timeStep());
    mlir::Value currentTime = builder.create<LoadOp>(loc, time);
    mlir::Value increasedTime = builder.create<AddOp>(loc, currentTime.getType(), currentTime, timeStep);
    builder.create<StoreOp>(loc, increasedTime, time);

    // Check if the current time is less than the end time
    mlir::Value endTime = builder.create<ConstantOp>(loc, modelOp.endTime());

    mlir::Value condition = builder.create<LteOp>(loc, BooleanType::get(modelOp->getContext()), increasedTime, endTime);
    condition = typeConverter->materializeTargetConversion(builder, condition.getLoc(), builder.getI1Type(), condition);
    builder.create<mlir::ReturnOp>(loc, condition);

    return mlir::success();
  }

  void ModelConverter::printSeparator(mlir::OpBuilder& builder, mlir::Value separator) const
  {
    auto module = separator.getParentRegion()->getParentOfType<mlir::ModuleOp>();
    auto printfRef = getOrInsertPrintf(builder, module);
    builder.create<mlir::LLVM::CallOp>(separator.getLoc(), printfRef, separator);
  }

  mlir::Value ModelConverter::getOrCreateGlobalString(
      mlir::Location loc,
      mlir::OpBuilder& builder,
      mlir::StringRef name,
      mlir::StringRef value,
      mlir::ModuleOp module) const
  {
    // Create the global at the entry of the module
    mlir::LLVM::GlobalOp global;

    if (!(global = module.lookupSymbol<mlir::LLVM::GlobalOp>(name))) {
      mlir::OpBuilder::InsertionGuard insertGuard(builder);
      builder.setInsertionPointToStart(module.getBody());
      auto type = mlir::LLVM::LLVMArrayType::get(mlir::IntegerType::get(builder.getContext(), 8), value.size());
      global = builder.create<mlir::LLVM::GlobalOp>(loc, type, true, mlir::LLVM::Linkage::Internal, name, builder.getStringAttr(value));
    }

    // Get the pointer to the first character in the global string
    mlir::Value globalPtr = builder.create<mlir::LLVM::AddressOfOp>(loc, global);

    mlir::Value cst0 = builder.create<mlir::LLVM::ConstantOp>(
        loc,
        mlir::IntegerType::get(builder.getContext(), 64),
        builder.getIntegerAttr(builder.getIndexType(), 0));

    return builder.create<mlir::LLVM::GEPOp>(
        loc,
        mlir::LLVM::LLVMPointerType::get(mlir::IntegerType::get(builder.getContext(), 8)),
        globalPtr, llvm::ArrayRef<mlir::Value>({cst0, cst0}));
  }

  mlir::Value ModelConverter::getSeparatorString(mlir::Location loc, mlir::OpBuilder& builder, mlir::ModuleOp module) const
  {
    return getOrCreateGlobalString(loc, builder, "semicolon", mlir::StringRef(";\0", 2), module);
  }

  mlir::Value ModelConverter::getNewlineString(mlir::Location loc, mlir::OpBuilder& builder, mlir::ModuleOp module) const
  {
    return getOrCreateGlobalString(loc, builder, "newline", mlir::StringRef("\n\0", 2), module);
  }

  mlir::LLVM::LLVMFuncOp ModelConverter::getOrInsertPrintf(mlir::OpBuilder& builder, mlir::ModuleOp module) const
  {
    auto *context = module.getContext();

    // Create a function declaration for printf, the signature is:
    //   * `i32 (i8*, ...)`
    auto llvmI32Ty = mlir::IntegerType::get(context, 32);
    auto llvmI8PtrTy = mlir::LLVM::LLVMPointerType::get(mlir::IntegerType::get(context, 8));
    auto llvmFnType = mlir::LLVM::LLVMFunctionType::get(llvmI32Ty, llvmI8PtrTy, true);

    // Insert the printf function into the body of the parent module
    return getOrInsertFunction(builder, module, "printf", llvmFnType);
  }

  void ModelConverter::printVariableName(
      mlir::OpBuilder& builder,
      mlir::Value name,
      mlir::Type type,
      VariableFilter::Filter filter,
      std::function<mlir::Value()> structValue,
      unsigned int position,
      mlir::ModuleOp module,
      mlir::Value separator,
      bool shouldPreprendSeparator) const
  {
    if (auto arrayType = type.dyn_cast<ArrayType>()) {
      if (arrayType.getRank() == 0) {
        printScalarVariableName(builder, name, module, separator, shouldPreprendSeparator);
      } else {
        printArrayVariableName(
            builder, name, type, filter, structValue, position, module, separator, shouldPreprendSeparator);
      }
    } else {
      printScalarVariableName(builder, name, module, separator, shouldPreprendSeparator);
    }
  }

  void ModelConverter::printScalarVariableName(
      mlir::OpBuilder& builder,
      mlir::Value name,
      mlir::ModuleOp module,
      mlir::Value separator,
      bool shouldPrependSeparator) const
  {
    if (shouldPrependSeparator) {
      printSeparator(builder, separator);
    }

    mlir::Location loc = name.getLoc();
    mlir::Value formatSpecifier = getOrCreateGlobalString(loc, builder, "frmt_spec_str", mlir::StringRef("%s\0", 3), module);
    auto printfRef = getOrInsertPrintf(builder, module);
    builder.create<mlir::LLVM::CallOp>(loc, printfRef, mlir::ValueRange({ formatSpecifier, name }));
  }

  void ModelConverter::printArrayVariableName(
      mlir::OpBuilder& builder,
      mlir::Value name,
      mlir::Type type,
      VariableFilter::Filter filter,
      std::function<mlir::Value()> structValue,
      unsigned int position,
      mlir::ModuleOp module,
      mlir::Value separator,
      bool shouldPrependSeparator) const
  {
    mlir::Location loc = name.getLoc();
    assert(type.isa<ArrayType>());

    // Get a reference to the printf function
    auto printfRef = getOrInsertPrintf(builder, module);

    // Create the brackets and comma strings
    mlir::Value lSquare = getOrCreateGlobalString(loc, builder, "lsquare", llvm::StringRef("[\0", 2), module);
    mlir::Value rSquare = getOrCreateGlobalString(loc, builder, "rsquare", llvm::StringRef("]\0", 2), module);
    mlir::Value comma = getOrCreateGlobalString(loc, builder, "comma", llvm::StringRef(",\0", 2), module);

    // Create the format strings
    mlir::Value stringFormatSpecifier = getOrCreateGlobalString(loc, builder, "frmt_spec_str", mlir::StringRef("%s\0", 3), module);
    mlir::Value integerFormatSpecifier = getOrCreateGlobalString(loc, builder, "frmt_spec_int", mlir::StringRef("%ld\0", 4), module);

    // Allow for the variable to lazily extracted if one of its dimension size
    // must be determined.
    bool valueLoaded = false;
    mlir::Value extractedValue = nullptr;
    auto insertionPoint = builder.saveInsertionPoint();

    auto var = [&]() -> mlir::Value {
      if (!valueLoaded) {
        mlir::OpBuilder::InsertionGuard guard(builder);
        builder.restoreInsertionPoint(insertionPoint);
        extractedValue = extractValue(builder, structValue(), type, position);
        valueLoaded = true;
      }

      return extractedValue;
    };

    // Create the lower and upper bounds
    auto ranges = filter.getRanges();
    auto arrayType = type.cast<ArrayType>();
    assert(arrayType.getRank() == ranges.size());

    llvm::SmallVector<mlir::Value, 3> lowerBounds;
    llvm::SmallVector<mlir::Value, 3> upperBounds;

    mlir::Value one = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(1));
    llvm::SmallVector<mlir::Value, 3> steps(arrayType.getRank(), one);

    for (const auto& range : llvm::enumerate(ranges)) {
      // In Modelica, arrays are 1-based. If present, we need to lower by 1
      // the value given by the variable filter.

      unsigned int lowerBound = range.value().hasLowerBound() ? range.value().getLowerBound() - 1 : 0;
      lowerBounds.push_back(builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(lowerBound)));

      // The upper bound is not lowered because the SCF's for operation assumes
      // them as excluded.

      if (range.value().hasUpperBound()) {
        mlir::Value upperBound = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.value().getUpperBound()));
        upperBounds.push_back(upperBound);
      } else {
        mlir::Value dim = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.index()));
        mlir::Value upperBound = builder.create<DimOp>(loc, var(), dim);
        upperBounds.push_back(upperBound);
      }
    }

    bool shouldPrintSeparator = false;

    // Create nested loops in order to iterate on each dimension of the array
    mlir::scf::buildLoopNest(
        builder, loc, lowerBounds, upperBounds, steps,
        [&](mlir::OpBuilder& nestedBuilder, mlir::Location loc, mlir::ValueRange indexes) {
          // Print the separator, the variable name and the left square bracket
          printSeparator(builder, separator);
          builder.create<mlir::LLVM::CallOp>(loc, printfRef, mlir::ValueRange({ stringFormatSpecifier, name }));
          builder.create<mlir::LLVM::CallOp>(loc, printfRef, lSquare);

          for (mlir::Value index : indexes) {
            if (shouldPrintSeparator)
              builder.create<mlir::LLVM::CallOp>(loc, printfRef, comma);

            shouldPrintSeparator = true;

            mlir::Type convertedType = typeConverter->convertType(index.getType());
            index = typeConverter->materializeTargetConversion(builder, loc, convertedType, index);

            // Arrays are 1-based in Modelica, so we add 1 in order to print
            // indexes that are coherent with the model source.
            mlir::Value increment = builder.create<mlir::ConstantOp>(loc, builder.getIntegerAttr(index.getType(), 1));
            index = builder.create<mlir::AddIOp>(loc, index.getType(), index, increment);

            builder.create<mlir::LLVM::CallOp>(loc, printfRef, mlir::ValueRange({ integerFormatSpecifier, index }));
          }

          // Print the right square bracket
          builder.create<mlir::LLVM::CallOp>(loc, printfRef, rSquare);
        });
  }

  mlir::LogicalResult ModelConverter::createPrintHeaderFunction(
      mlir::OpBuilder& builder,
      ModelOp op,
      DerivativesPositionsMap& derivativesPositions) const
  {
    mlir::TypeRange varTypes = op.body().getArgumentTypes();

    auto callback = [&](std::function<mlir::Value()> structValue, llvm::StringRef name, unsigned int position, VariableFilter::Filter filter, mlir::Value separator) -> mlir::LogicalResult {
      mlir::Location loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      std::string symbolName = "var" + std::to_string(position);
      llvm::SmallString<10> terminatedName(name);
      terminatedName.append("\0");
      mlir::Value symbol = getOrCreateGlobalString(loc, builder, symbolName, llvm::StringRef(terminatedName.c_str(), terminatedName.size() + 1), module);

      bool shouldPrintSeparator = position != 0;
      printVariableName(builder, symbol, varTypes[position], filter, structValue, position, module, separator, shouldPrintSeparator);
      return mlir::success();
    };

    return createPrintFunctionBody(builder, op, varTypes, derivativesPositions, printHeaderFunctionName, callback);
  }

  void ModelConverter::printVariable(
      mlir::OpBuilder& builder,
      mlir::Value var,
      VariableFilter::Filter filter,
      mlir::Value separator,
      bool shouldPreprendSeparator) const
  {
    if (auto arrayType = var.getType().dyn_cast<ArrayType>()) {
      if (arrayType.getRank() == 0) {
        mlir::Value value = builder.create<LoadOp>(var.getLoc(), var);
        printScalarVariable(builder, value, separator, shouldPreprendSeparator);
      } else {
        printArrayVariable(builder, var, filter, separator, shouldPreprendSeparator);
      }
    } else {
      printScalarVariable(builder, var, separator, shouldPreprendSeparator);
    }
  }

  void ModelConverter::printScalarVariable(
      mlir::OpBuilder& builder, mlir::Value var, mlir::Value separator, bool shouldPreprendSeparator) const
  {
    if (shouldPreprendSeparator) {
      printSeparator(builder, separator);
    }

    printElement(builder, var);
  }

  void ModelConverter::printArrayVariable(
      mlir::OpBuilder& builder,
      mlir::Value var,
      VariableFilter::Filter filter,
      mlir::Value separator,
      bool shouldPreprendSeparator) const
  {
    mlir::Location loc = var.getLoc();
    assert(var.getType().isa<ArrayType>());

    auto ranges = filter.getRanges();
    auto arrayType = var.getType().cast<ArrayType>();
    assert(arrayType.getRank() == ranges.size());

    llvm::SmallVector<mlir::Value, 3> lowerBounds;
    llvm::SmallVector<mlir::Value, 3> upperBounds;

    mlir::Value one = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(1));
    llvm::SmallVector<mlir::Value, 3> steps(arrayType.getRank(), one);

    for (const auto& range : llvm::enumerate(ranges)) {
      // In Modelica, arrays are 1-based. If present, we need to lower by 1
      // the value given by the variable filter.

      unsigned int lowerBound = range.value().hasLowerBound() ? range.value().getLowerBound() - 1 : 0;
      lowerBounds.push_back(builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(lowerBound)));

      // The upper bound is not lowered because the SCF's for operation assumes
      // them as excluded.

      if (range.value().hasUpperBound()) {
        mlir::Value upperBound = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.value().getUpperBound()));
        upperBounds.push_back(upperBound);
      } else {
        mlir::Value dim = builder.create<mlir::ConstantOp>(loc, builder.getIndexAttr(range.index()));
        mlir::Value upperBound = builder.create<DimOp>(loc, var, dim);
        upperBounds.push_back(upperBound);
      }
    }

    // Create nested loops in order to iterate on each dimension of the array
    mlir::scf::buildLoopNest(
        builder, loc, lowerBounds, upperBounds, steps,
        [&](mlir::OpBuilder& nestedBuilder, mlir::Location loc, mlir::ValueRange position) {
          mlir::Value value = nestedBuilder.create<LoadOp>(loc, var, position);

          printSeparator(nestedBuilder, separator);
          printElement(nestedBuilder, value);
        });
  }

  void ModelConverter::printElement(mlir::OpBuilder& builder, mlir::Value value) const
  {
    mlir::Location loc = value.getLoc();
    auto module = value.getParentRegion()->getParentOfType<mlir::ModuleOp>();
    auto printfRef = getOrInsertPrintf(builder, module);

    mlir::Type convertedType = typeConverter->convertType(value.getType());
    value = typeConverter->materializeTargetConversion(builder, loc, convertedType, value);
    mlir::Type type = value.getType();

    mlir::Value formatSpecifier;

    if (type.isa<mlir::IntegerType>()) {
      formatSpecifier = getOrCreateGlobalString(
          loc, builder, "frmt_spec_int", mlir::StringRef("%ld\0", 4), module);
    } else if (type.isa<mlir::FloatType>()) {
      formatSpecifier = getOrCreateGlobalString(
          loc, builder, "frmt_spec_float", mlir::StringRef("%.12f\0", 6), module);
    } else {
      assert(false && "Unknown type");
    }

    builder.create<mlir::LLVM::CallOp>(value.getLoc(), printfRef, mlir::ValueRange({ formatSpecifier, value }));
  }

  mlir::LogicalResult ModelConverter::createPrintFunction(
      mlir::OpBuilder& builder,
      ModelOp op,
      DerivativesPositionsMap& derivativesPositions) const
  {
    mlir::TypeRange varTypes = op.body().getArgumentTypes();

    auto callback = [&](std::function<mlir::Value()> structValue, llvm::StringRef name, unsigned int position, VariableFilter::Filter filter, mlir::Value separator) -> mlir::LogicalResult {
      mlir::Value var = extractValue(builder, structValue(), varTypes[position], position);
      bool shouldPrintSeparator = position != 0;
      printVariable(builder, var, filter, separator, shouldPrintSeparator);
      return mlir::success();
    };

    return createPrintFunctionBody(builder, op, varTypes, derivativesPositions, printFunctionName, callback);
  }

  mlir::LogicalResult ModelConverter::createPrintFunctionBody(
      mlir::OpBuilder& builder,
      ModelOp op,
      mlir::TypeRange varTypes,
      DerivativesPositionsMap& derivativesPositions,
      llvm::StringRef functionName,
      std::function<mlir::LogicalResult(std::function<mlir::Value()>, llvm::StringRef, unsigned int, VariableFilter::Filter, mlir::Value)> elementCallback) const
  {
    mlir::Location loc = op.getLoc();
    mlir::OpBuilder::InsertionGuard guard(builder);
    auto module = op->getParentOfType<mlir::ModuleOp>();

    // Create the function inside the parent module
    builder.setInsertionPointToEnd(module.getBody());

    auto function = builder.create<mlir::FuncOp>(
        loc, functionName,
        builder.getFunctionType(getVoidPtrType(), llvm::None));

    auto* entryBlock = function.addEntryBlock();
    builder.setInsertionPointToStart(entryBlock);

    // Create the separator and newline global strings
    mlir::Value separator = getSeparatorString(loc, builder, module);
    mlir::Value newline = getNewlineString(loc, builder, module);

    // Create the callback to load the data structure whenever needed
    bool structValueLoaded = false;
    mlir::Value structValue = nullptr;
    auto structValueInsertionPoint = builder.saveInsertionPoint();

    auto structValueCallback = [&]() -> mlir::Value {
      if (!structValueLoaded) {
        mlir::OpBuilder::InsertionGuard guard(builder);
        builder.restoreInsertionPoint(structValueInsertionPoint);
        structValue = loadDataFromOpaquePtr(builder, function.getArgument(0), varTypes);
      }

      return structValue;
    };

    // Get the names of the variables
    llvm::SmallVector<llvm::StringRef, 8> variableNames =
        llvm::to_vector<8>(op.variableNames().getAsValueRange<mlir::StringAttr>());

    // Map each variable to its position inside the data structure.
    // It must be noted that the data structure also contains derivative (if
    // existent), so its size can be greater than the number of names.

    assert(op.variableNames().size() <= varTypes.size());
    llvm::StringMap<size_t> variablePositionByName;

    for (const auto& var : llvm::enumerate(variableNames))
      variablePositionByName[var.value()] = var.index() + 1; // + 1 to skip the "time" variable

    // The positions have been saved, so we can now sort the names
    llvm::sort(variableNames, [](llvm::StringRef x, llvm::StringRef y) -> bool {
      return x.compare_insensitive(y) < 0;
    });

    if (auto status = elementCallback(
          structValueCallback, "time", timeVariablePosition, VariableFilter::Filter::visibleScalar(), separator); mlir::failed(status)) {
      return status;
    }

    // Print the other variables
    for (const auto& name : variableNames) {
      assert(variablePositionByName.count(name) != 0);
      size_t position = variablePositionByName[name];

      unsigned int rank = 0;

      if (auto arrayType = varTypes[position].dyn_cast<ArrayType>()) {
        rank = arrayType.getRank();
      }

      auto filter = options.variableFilter->getVariableInfo(name, rank);

      if (!filter.isVisible()) {
        continue;
      }

      if (auto status = elementCallback(
            structValueCallback, name, position, filter, separator); mlir::failed(status)) {
        return status;
      }
    }

    // Print the derivatives
    for (const auto& name : variableNames) {
      size_t varPosition = variablePositionByName[name];

      if (derivativesPositions.count(varPosition) == 0) {
        // The variable has no derivative
        continue;
      }

      size_t derivedVarPosition = derivativesPositions[varPosition];

      unsigned int rank = 0;

      if (auto arrayType = varTypes[derivedVarPosition].dyn_cast<ArrayType>()) {
        rank = arrayType.getRank();
      }

      auto filter = options.variableFilter->getVariableDerInfo(name, rank);

      if (!filter.isVisible()) {
        continue;
      }

      llvm::SmallString<15> derName;
      derName.append("der(");
      derName.append(name);
      derName.append(")");

      if (auto status = elementCallback(
            structValueCallback, derName, derivedVarPosition, filter, separator); mlir::failed(status)) {
        return status;
      }
    }

    // Print a newline character after all the variables have been processed
    builder.create<mlir::LLVM::CallOp>(loc, getOrInsertPrintf(builder, module), newline);

    builder.create<mlir::ReturnOp>(loc);
    return mlir::success();
  }
}
