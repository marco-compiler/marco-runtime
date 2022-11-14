#include "marco/Codegen/Conversion/KINSOLToLLVM/KINSOLToLLVM.h"
#include "marco/Codegen/Conversion/KINSOLToLLVM/LLVMTypeConverter.h"
#include "marco/Codegen/Runtime.h"
#include "marco/Dialect/Modelica/ModelicaDialect.h"
#include "mlir/Dialect/Arith/IR/Arith.h"
#include "mlir/Dialect/LLVMIR/LLVMDialect.h"
#include "mlir/Dialect/LLVMIR/FunctionCallUtils.h"
#include "mlir/Transforms/DialectConversion.h"
#include "mlir/Conversion/LLVMCommon/Pattern.h"

namespace mlir
{
#define GEN_PASS_DEF_KINSOLTOLLVMCONVERSIONPASS
#include "marco/Codegen/Conversion/Passes.h.inc"
}

using namespace ::marco;
using namespace ::marco::codegen;
using namespace ::mlir::kinsol;

static mlir::LLVM::LLVMFuncOp getOrDeclareFunction(
    mlir::OpBuilder& builder,
    mlir::ModuleOp module,
    mlir::Location loc,
    llvm::StringRef name,
    mlir::Type result,
    llvm::ArrayRef<mlir::Type> args)
{
  if (auto funcOp = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(name)) {
    return funcOp;
  }

  mlir::PatternRewriter::InsertionGuard guard(builder);
  builder.setInsertionPointToStart(module.getBody());

  return builder.create<mlir::LLVM::LLVMFuncOp>(loc, name, mlir::LLVM::LLVMFunctionType::get(result, args));
}

static mlir::LLVM::LLVMFuncOp getOrDeclareFunction(
    mlir::OpBuilder& builder,
    mlir::ModuleOp module,
    mlir::Location loc,
    llvm::StringRef name,
    mlir::Type result,
    mlir::ValueRange args)
{
  llvm::SmallVector<mlir::Type, 3> argsTypes;

  for (auto type : args.getTypes()) {
    argsTypes.push_back(type);
  }

  return getOrDeclareFunction(builder, module, loc, name, result, argsTypes);
}

namespace
{
  /// Generic conversion pattern that provides some utility functions.
  template<typename Op>
  class KINSOLOpConversion : public mlir::ConvertOpToLLVMPattern<Op>
  {
    public:
      KINSOLOpConversion(mlir::LLVMTypeConverter& typeConverter, unsigned int bitWidth)
          : mlir::ConvertOpToLLVMPattern<Op>(typeConverter),
            bitWidth(bitWidth)
      {
      }

      mlir::kinsol::LLVMTypeConverter& typeConverter() const
      {
        return *static_cast<mlir::kinsol::LLVMTypeConverter*>(this->getTypeConverter());
      }

      mlir::Type convertType(mlir::Type type) const
      {
        return typeConverter().convertType(type);
      }

      mlir::Value materializeTargetConversion(mlir::OpBuilder& builder, mlir::Value value) const
      {
        mlir::Type type = this->getTypeConverter()->convertType(value.getType());
        return this->getTypeConverter()->materializeTargetConversion(builder, value.getLoc(), type, value);
      }

      void materializeTargetConversion(mlir::OpBuilder& builder, llvm::SmallVectorImpl<mlir::Value>& values) const
      {
        for (auto& value : values) {
          value = materializeTargetConversion(builder, value);
        }
      }

    protected:
      unsigned int bitWidth;
  };

  struct CreateOpLowering : public KINSOLOpConversion<CreateOp>
  {
    using KINSOLOpConversion<CreateOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(CreateOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 2> newOperands;
      llvm::SmallVector<std::string, 2> mangledArgsTypes;

      // Scalar equations amount
      mlir::Value scalarEquationsAmount = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(op.getScalarEquations()));
      newOperands.push_back(scalarEquationsAmount);
      mangledArgsTypes.push_back(mangling.getIntegerType(scalarEquationsAmount.getType().getIntOrFloatBitWidth()));

      // Data bit-width
      mlir::Value dataBitWidth = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(bitWidth));
      newOperands.push_back(dataBitWidth);
      mangledArgsTypes.push_back(mangling.getIntegerType(dataBitWidth.getType().getIntOrFloatBitWidth()));

      // Create the call to the runtime library
      auto resultType = getVoidPtrType();
      auto mangledResultType = mangling.getVoidPointerType();
      auto functionName = mangling.getMangledFunction("kinsolCreate", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, newOperands);

      return mlir::success();
    }
  };

  struct AddEquationOpLowering : public KINSOLOpConversion<AddEquationOp>
  {
    using KINSOLOpConversion<AddEquationOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(AddEquationOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 3> newOperands;
      llvm::SmallVector<std::string, 3> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the array with the equation ranges
      mlir::Type dimensionSizeType = getTypeConverter()->convertType(rewriter.getI64Type());
      mlir::Value numOfElements = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), op.getEquationRanges().size() * 2));
      mlir::Type elementPtrType = mlir::LLVM::LLVMPointerType::get(dimensionSizeType);
      mlir::Value nullPtr = rewriter.create<mlir::LLVM::NullOp>(loc, elementPtrType);
      mlir::Value gepPtr = rewriter.create<mlir::LLVM::GEPOp>(loc, elementPtrType, nullPtr, numOfElements);
      mlir::Value sizeBytes = rewriter.create<mlir::LLVM::PtrToIntOp>(loc, getIndexType(), gepPtr);

      auto heapAllocFn = lookupOrCreateHeapAllocFn(module, getIndexType());
      mlir::Value equationRangesOpaquePtr = rewriter.create<mlir::LLVM::CallOp>(loc, heapAllocFn, sizeBytes).getResult();
      mlir::Value equationRangesPtr = rewriter.create<mlir::LLVM::BitcastOp>(loc, elementPtrType, equationRangesOpaquePtr);

      newOperands.push_back(equationRangesPtr);
      mangledArgsTypes.push_back(mangling.getPointerType(mangling.getIntegerType(dimensionSizeType.getIntOrFloatBitWidth())));

      // Populate the equation ranges
      for (const auto& range : llvm::enumerate(op.getEquationRanges())) {
        auto rangeAttr = range.value().cast<mlir::ArrayAttr>();
        assert(rangeAttr.size() == 2);

        for (const auto& index : llvm::enumerate(rangeAttr)) {
          auto indexAttr = index.value().cast<mlir::IntegerAttr>();
          mlir::Value offset = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), range.index() * 2 + index.index()));
          mlir::Value indexValue = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(indexAttr.getInt()));
          mlir::Value ptr = rewriter.create<mlir::LLVM::GEPOp>(loc, equationRangesPtr.getType(), equationRangesPtr, offset);
          rewriter.create<mlir::LLVM::StoreOp>(loc, indexValue, ptr);
        }
      }

      // Rank
      mlir::Value rank = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(op.getEquationRanges().size()));
      newOperands.push_back(rank);
      mangledArgsTypes.push_back(mangling.getIntegerType(rank.getType().getIntOrFloatBitWidth()));

      // Create the call to the runtime library
      auto resultType = getTypeConverter()->convertType(op.getResult().getType());
      auto mangledResultType = mangling.getIntegerType(resultType.getIntOrFloatBitWidth());
      auto functionName = mangling.getMangledFunction("kinsolAddEquation", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, newOperands);

      // Deallocate the ranges array
      auto heapFreeFn = lookupOrCreateHeapFreeFn(module);
      rewriter.create<mlir::LLVM::CallOp>(loc, heapFreeFn, equationRangesOpaquePtr);

      return mlir::success();
    }
  };

  struct AddAlgebraicVariableOpLowering : public KINSOLOpConversion<AddAlgebraicVariableOp>
  {
    using KINSOLOpConversion<AddAlgebraicVariableOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(AddAlgebraicVariableOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      auto heapAllocFn = lookupOrCreateHeapAllocFn(module, getIndexType());

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 5> newOperands;
      llvm::SmallVector<std::string, 3> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Variable
      mlir::Type variableType = adaptor.getVariable().getType();
      mlir::Type variablePtrType = mlir::LLVM::LLVMPointerType::get(variableType);
      mlir::Value variableNullPtr = rewriter.create<mlir::LLVM::NullOp>(loc, variablePtrType);
      mlir::Value one = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), 1));
      mlir::Value variableGepPtr = rewriter.create<mlir::LLVM::GEPOp>(loc, variablePtrType, variableNullPtr, one);
      mlir::Value variableSizeBytes = rewriter.create<mlir::LLVM::PtrToIntOp>(loc, getIndexType(), variableGepPtr);

      mlir::Value variableOpaquePtr = rewriter.create<mlir::LLVM::CallOp>(loc, heapAllocFn, variableSizeBytes).getResult();
      mlir::Value variablePtr = rewriter.create<mlir::LLVM::BitcastOp>(loc, variablePtrType, variableOpaquePtr);
      rewriter.create<mlir::LLVM::StoreOp>(loc, adaptor.getVariable(), variablePtr);

      newOperands.push_back(variableOpaquePtr);
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the array with the variable dimensions
      mlir::Type dimensionSizeType = getTypeConverter()->convertType(rewriter.getI64Type());
      mlir::Value numOfElements = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), op.getDimensions().size()));
      mlir::Type elementPtrType = mlir::LLVM::LLVMPointerType::get(dimensionSizeType);
      mlir::Value nullPtr = rewriter.create<mlir::LLVM::NullOp>(loc, elementPtrType);
      mlir::Value gepPtr = rewriter.create<mlir::LLVM::GEPOp>(loc, elementPtrType, nullPtr, numOfElements);
      mlir::Value sizeBytes = rewriter.create<mlir::LLVM::PtrToIntOp>(loc, getIndexType(), gepPtr);

      mlir::Value arrayDimensionsOpaquePtr = rewriter.create<mlir::LLVM::CallOp>(loc, heapAllocFn, sizeBytes).getResult();
      mlir::Value arrayDimensionsPtr = rewriter.create<mlir::LLVM::BitcastOp>(loc, elementPtrType, arrayDimensionsOpaquePtr);

      newOperands.push_back(arrayDimensionsPtr);
      mangledArgsTypes.push_back(mangling.getPointerType(mangling.getIntegerType(dimensionSizeType.getIntOrFloatBitWidth())));

      // Populate the dimensions list
      for (const auto& sizeAttr : llvm::enumerate(op.getDimensions())) {
        mlir::Value offset = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), sizeAttr.index()));
        mlir::Value size = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(sizeAttr.value().cast<mlir::IntegerAttr>().getInt()));
        mlir::Value ptr = rewriter.create<mlir::LLVM::GEPOp>(loc, arrayDimensionsPtr.getType(), arrayDimensionsPtr, offset);
        rewriter.create<mlir::LLVM::StoreOp>(loc, size, ptr);
      }

      // Rank
      mlir::Value rank = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(op.getDimensions().size()));
      newOperands.push_back(rank);
      mangledArgsTypes.push_back(mangling.getIntegerType(rank.getType().getIntOrFloatBitWidth()));

      // Variable getter function address
      auto getter = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(op.getGetter());
      mlir::Value getterAddress = rewriter.create<mlir::LLVM::AddressOfOp>(loc, getter);
      getterAddress = rewriter.create<mlir::LLVM::BitcastOp>(loc, getVoidPtrType(), getterAddress);

      newOperands.push_back(getterAddress);
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Variable setter function address
      auto setter = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(op.getSetter());
      mlir::Value setterAddress = rewriter.create<mlir::LLVM::AddressOfOp>(loc, setter);
      setterAddress = rewriter.create<mlir::LLVM::BitcastOp>(loc, getVoidPtrType(), setterAddress);

      newOperands.push_back(setterAddress);
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the call to the runtime library
      auto resultType = getTypeConverter()->convertType(op.getResult().getType());
      auto mangledResultType = mangling.getIntegerType(resultType.getIntOrFloatBitWidth());
      auto functionName = mangling.getMangledFunction("kinsolAddAlgebraicVariable", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, newOperands);

      // Deallocate the dimensions array
      auto heapFreeFn = lookupOrCreateHeapFreeFn(module);
      rewriter.create<mlir::LLVM::CallOp>(loc, heapFreeFn, arrayDimensionsOpaquePtr);

      return mlir::success();
    }
  };

  struct VariableGetterOpLowering : public KINSOLOpConversion<VariableGetterOp>
  {
    using KINSOLOpConversion<VariableGetterOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(VariableGetterOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();

      mlir::SmallVector<mlir::Type, 3> argsTypes;

      argsTypes.push_back(getVoidPtrType());
      argsTypes.push_back(mlir::LLVM::LLVMPointerType::get(getTypeConverter()->getIndexType()));

      auto functionType = mlir::LLVM::LLVMFunctionType::get(getTypeConverter()->convertType(op.getFunctionType().getResult(0)), argsTypes);

      auto newOp = rewriter.replaceOpWithNewOp<mlir::LLVM::LLVMFuncOp>(op, op.getSymName(), functionType);
      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      // The lowered function will receive a void pointer to the array descriptor of the variable
      mlir::Value variableOpaquePtr = newOp.getArgument(0);
      mlir::Value variablePtr = rewriter.create<mlir::LLVM::BitcastOp>(loc, mlir::LLVM::LLVMPointerType::get(op.getVariable().getType()), variableOpaquePtr);
      mlir::Value variable = rewriter.create<mlir::LLVM::LoadOp>(loc, variablePtr);
      mapping.map(op.getVariable(), variable);

      // The equation indices are also passed through an array
      mlir::Value variableIndicesPtr = newOp.getArgument(1);

      for (auto variableIndex : llvm::enumerate(op.getVariableIndices())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(variableIndicesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), variableIndex.index()));
        mlir::Value variableIndexPtr = rewriter.create<mlir::LLVM::GEPOp>(variableIndicesPtr.getLoc(), variableIndicesPtr.getType(), variableIndicesPtr, index);
        mlir::Value mappedVariableIndex = rewriter.create<mlir::LLVM::LoadOp>(variableIndexPtr.getLoc(), variableIndexPtr);
        mappedVariableIndex = getTypeConverter()->materializeSourceConversion(rewriter, mappedVariableIndex.getLoc(), rewriter.getIndexType(), mappedVariableIndex);
        mapping.map(variableIndex.value(), mappedVariableIndex);
      }

      // Clone the blocks structure
      for (auto& block : llvm::enumerate(op.getBodyRegion().getBlocks())) {
        if (block.index() != 0) {
          std::vector<mlir::Location> argLocations;

          for (const auto& arg : block.value().getArguments()) {
            argLocations.push_back(arg.getLoc());
          }

          mlir::Block* clonedBlock = rewriter.createBlock(
              &newOp.getBody(),
              newOp.getBody().end(),
              block.value().getArgumentTypes(),
              argLocations);

          mapping.map(&block.value(), clonedBlock);

          for (const auto& [original, cloned] : llvm::zip(block.value().getArguments(), clonedBlock->getArguments())) {
            mapping.map(original, cloned);
          }
        }
      }

      // Clone the original operations
      for (auto& block : llvm::enumerate(op.getBodyRegion())) {
        if (block.index() == 0) {
          rewriter.setInsertionPointToEnd(entryBlock);
        } else {
          rewriter.setInsertionPointToStart(mapping.lookup(&block.value()));
        }

        for (auto& bodyOp : block.value().getOperations()) {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };

  struct VariableSetterOpLowering : public KINSOLOpConversion<VariableSetterOp>
  {
    using KINSOLOpConversion<VariableSetterOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(VariableSetterOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();

      mlir::SmallVector<mlir::Type, 3> argsTypes;

      argsTypes.push_back(getVoidPtrType());
      argsTypes.push_back(op.getValue().getType());
      argsTypes.push_back(mlir::LLVM::LLVMPointerType::get(getTypeConverter()->getIndexType()));

      auto functionType = mlir::LLVM::LLVMFunctionType::get(mlir::LLVM::LLVMVoidType::get(rewriter.getContext()), argsTypes);

      auto newOp = rewriter.replaceOpWithNewOp<mlir::LLVM::LLVMFuncOp>(op, op.getSymName(), functionType);
      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      // The lowered function will receive a void pointer to the array descriptor of the variable
      mlir::Value variableOpaquePtr = newOp.getArgument(0);
      mlir::Value variablePtr = rewriter.create<mlir::LLVM::BitcastOp>(loc, mlir::LLVM::LLVMPointerType::get(op.getVariable().getType()), variableOpaquePtr);
      mlir::Value variable = rewriter.create<mlir::LLVM::LoadOp>(loc, variablePtr);
      mapping.map(op.getVariable(), variable);

      // Map the value
      mapping.map(op.getValue(), newOp.getArgument(1));

      // The equation indices are also passed through an array
      mlir::Value variableIndicesPtr = newOp.getArgument(2);

      for (auto variableIndex : llvm::enumerate(op.getVariableIndices())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(variableIndicesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), variableIndex.index()));
        mlir::Value variableIndexPtr = rewriter.create<mlir::LLVM::GEPOp>(variableIndicesPtr.getLoc(), variableIndicesPtr.getType(), variableIndicesPtr, index);
        mlir::Value mappedVariableIndex = rewriter.create<mlir::LLVM::LoadOp>(variableIndexPtr.getLoc(), variableIndexPtr);
        mappedVariableIndex = getTypeConverter()->materializeSourceConversion(rewriter, mappedVariableIndex.getLoc(), rewriter.getIndexType(), mappedVariableIndex);
        mapping.map(variableIndex.value(), mappedVariableIndex);
      }

      // Clone the blocks structure
      for (auto& block : llvm::enumerate(op.getBodyRegion().getBlocks())) {
        if (block.index() != 0) {
          std::vector<mlir::Location> argLocations;

          for (const auto& arg : block.value().getArguments()) {
            argLocations.push_back(arg.getLoc());
          }

          mlir::Block* clonedBlock = rewriter.createBlock(
              &newOp.getBody(),
              newOp.getBody().end(),
              block.value().getArgumentTypes(),
              argLocations);

          mapping.map(&block.value(), clonedBlock);

          for (const auto& [original, cloned] : llvm::zip(block.value().getArguments(), clonedBlock->getArguments())) {
            mapping.map(original, cloned);
          }
        }
      }

      // Clone the original operations
      for (auto& block : llvm::enumerate(op.getBodyRegion())) {
        if (block.index() == 0) {
          rewriter.setInsertionPointToEnd(entryBlock);
        } else {
          rewriter.setInsertionPointToStart(mapping.lookup(&block.value()));
        }

        for (auto& bodyOp : block.value().getOperations()) {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };

  struct AddVariableAccessOpLowering : public KINSOLOpConversion<AddVariableAccessOp>
  {
    using KINSOLOpConversion<AddVariableAccessOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(AddVariableAccessOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 5> newOperands;
      llvm::SmallVector<std::string, 3> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Equation
      newOperands.push_back(adaptor.getEquation());
      mangledArgsTypes.push_back(mangling.getIntegerType(adaptor.getEquation().getType().getIntOrFloatBitWidth()));

      // Variable
      newOperands.push_back(adaptor.getVariable());
      mangledArgsTypes.push_back(mangling.getIntegerType(adaptor.getVariable().getType().getIntOrFloatBitWidth()));

      // Create the array with the variable accesses
      auto dimensionAccesses = op.getAccess().getResults();
      llvm::SmallVector<mlir::Value, 6> accessValues;

      for (const auto dimensionAccess : dimensionAccesses) {
        if (dimensionAccess.isa<mlir::AffineConstantExpr>()) {
          auto constantAccess = dimensionAccess.cast<mlir::AffineConstantExpr>();
          accessValues.push_back(rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(-1)));
          accessValues.push_back(rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(constantAccess.getValue())));
        } else if (dimensionAccess.isa<mlir::AffineDimExpr>()) {
          auto dimension = dimensionAccess.cast<mlir::AffineDimExpr>();
          accessValues.push_back(rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(dimension.getPosition())));
          accessValues.push_back(rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(0)));
        } else {
          auto dynamicAccess = dimensionAccess.cast<mlir::AffineBinaryOpExpr>();
          auto dimension = dynamicAccess.getLHS().cast<mlir::AffineDimExpr>();
          auto offset = dynamicAccess.getRHS().cast<mlir::AffineConstantExpr>();
          accessValues.push_back(rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(dimension.getPosition())));
          accessValues.push_back(rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(offset.getValue())));
        }
      }

      mlir::Type dimensionSizeType = getTypeConverter()->convertType(rewriter.getI64Type());
      mlir::Value numOfElements = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), accessValues.size()));
      mlir::Type elementPtrType = mlir::LLVM::LLVMPointerType::get(dimensionSizeType);
      mlir::Value nullPtr = rewriter.create<mlir::LLVM::NullOp>(loc, elementPtrType);
      mlir::Value gepPtr = rewriter.create<mlir::LLVM::GEPOp>(loc, elementPtrType, nullPtr, numOfElements);
      mlir::Value sizeBytes = rewriter.create<mlir::LLVM::PtrToIntOp>(loc, getIndexType(), gepPtr);

      auto heapAllocFn = lookupOrCreateHeapAllocFn(op->getParentOfType<mlir::ModuleOp>(), getIndexType());
      mlir::Value accessesOpaquePtr = rewriter.create<mlir::LLVM::CallOp>(loc, heapAllocFn, sizeBytes).getResult();
      mlir::Value accessesPtr = rewriter.create<mlir::LLVM::BitcastOp>(loc, elementPtrType, accessesOpaquePtr);

      newOperands.push_back(accessesPtr);
      mangledArgsTypes.push_back(mangling.getPointerType(mangling.getIntegerType(dimensionSizeType.getIntOrFloatBitWidth())));

      // Populate the equation ranges
      for (const auto& accessValue : llvm::enumerate(accessValues)) {
        mlir::Value offset = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getIntegerAttr(getTypeConverter()->getIndexType(), accessValue.index()));
        mlir::Value ptr = rewriter.create<mlir::LLVM::GEPOp>(loc, accessesPtr.getType(), accessesPtr, offset);
        rewriter.create<mlir::LLVM::StoreOp>(loc, accessValue.value(), ptr);
      }

      // Rank
      mlir::Value rank = rewriter.create<mlir::arith::ConstantOp>(loc, rewriter.getI64IntegerAttr(op.getAccess().getResults().size()));
      newOperands.push_back(rank);
      mangledArgsTypes.push_back(mangling.getIntegerType(rank.getType().getIntOrFloatBitWidth()));

      // Create the call to the runtime library
      auto resultType = getVoidType();
      auto mangledResultType = mangling.getVoidType();
      auto functionName = mangling.getMangledFunction("kinsolAddVariableAccess", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, newOperands);

      // Deallocate the accesses array
      auto heapFreeFn = lookupOrCreateHeapFreeFn(op->getParentOfType<mlir::ModuleOp>());
      rewriter.create<mlir::LLVM::CallOp>(loc, heapFreeFn, accessesOpaquePtr);

      return mlir::success();
    }
  };

  struct ResidualFunctionOpLowering : public KINSOLOpConversion<ResidualFunctionOp>
  {
    using KINSOLOpConversion<ResidualFunctionOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(ResidualFunctionOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      mlir::SmallVector<mlir::Type, 3> argsTypes;

      argsTypes.push_back(op.getTime().getType());
      argsTypes.push_back(getVoidPtrType());
      argsTypes.push_back(mlir::LLVM::LLVMPointerType::get(getTypeConverter()->getIndexType()));

      auto functionType = mlir::LLVM::LLVMFunctionType::get(op.getFunctionType().getResult(0), argsTypes);

      auto newOp = rewriter.replaceOpWithNewOp<mlir::LLVM::LLVMFuncOp>(op, op.getSymName(), functionType);
      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      // Map the time variable
      mapping.map(op.getTime(), newOp.getArgument(0));

      // The lowered function will receive a pointer to the array of variables.
      assert(!op.getVariables().empty());

      mlir::Value variablesPtr = newOp.getArgument(1);
      variablesPtr = rewriter.create<mlir::LLVM::BitcastOp>(variablesPtr.getLoc(), mlir::LLVM::LLVMPointerType::get(getVoidPtrType()), variablesPtr);

      for (auto variable : llvm::enumerate(op.getVariables())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(variablesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), variable.index()));
        mlir::Value variablePtr = rewriter.create<mlir::LLVM::GEPOp>(variablesPtr.getLoc(), variablesPtr.getType(), variablesPtr, index);
        mlir::Value mappedVariable = rewriter.create<mlir::LLVM::LoadOp>(variablePtr.getLoc(), variablePtr);
        mappedVariable = rewriter.create<mlir::LLVM::BitcastOp>(mappedVariable.getLoc(), mlir::LLVM::LLVMPointerType::get(op.getVariables()[variable.index()].getType()), mappedVariable);
        mappedVariable = rewriter.create<mlir::LLVM::LoadOp>(mappedVariable.getLoc(), mappedVariable);

        mapping.map(variable.value(), mappedVariable);
      }

      // The equation indices are also passed through an array
      mlir::Value equationIndicesPtr = newOp.getArgument(2);

      for (auto equationIndex : llvm::enumerate(op.getEquationIndices())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(equationIndicesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), equationIndex.index()));
        mlir::Value equationIndexPtr = rewriter.create<mlir::LLVM::GEPOp>(equationIndicesPtr.getLoc(), equationIndicesPtr.getType(), equationIndicesPtr, index);
        mlir::Value mappedEquationIndex = rewriter.create<mlir::LLVM::LoadOp>(equationIndexPtr.getLoc(), equationIndexPtr);
        mappedEquationIndex = getTypeConverter()->materializeSourceConversion(rewriter, mappedEquationIndex.getLoc(), rewriter.getIndexType(), mappedEquationIndex);
        mapping.map(equationIndex.value(), mappedEquationIndex);
      }

      // Clone the blocks structure
      for (auto& block : llvm::enumerate(op.getBodyRegion().getBlocks())) {
        if (block.index() != 0) {
          std::vector<mlir::Location> argLocations;

          for (const auto& arg : block.value().getArguments()) {
            argLocations.push_back(arg.getLoc());
          }

          mlir::Block* clonedBlock = rewriter.createBlock(
              &newOp.getBody(),
              newOp.getBody().end(),
              block.value().getArgumentTypes(),
              argLocations);

          mapping.map(&block.value(), clonedBlock);

          for (const auto& [original, cloned] : llvm::zip(block.value().getArguments(), clonedBlock->getArguments())) {
            mapping.map(original, cloned);
          }
        }
      }

      // Clone the original operations
      for (auto& block : llvm::enumerate(op.getBodyRegion())) {
        if (block.index() == 0) {
          rewriter.setInsertionPointToEnd(entryBlock);
        } else {
          rewriter.setInsertionPointToStart(mapping.lookup(&block.value()));
        }

        for (auto& bodyOp : block.value().getOperations()) {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };

  struct JacobianFunctionOpLowering : public KINSOLOpConversion<JacobianFunctionOp>
  {
    using KINSOLOpConversion<JacobianFunctionOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(JacobianFunctionOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      mlir::SmallVector<mlir::Type, 5> argsTypes;

      argsTypes.push_back(op.getTime().getType());
      argsTypes.push_back(getVoidPtrType());
      argsTypes.push_back(mlir::LLVM::LLVMPointerType::get(getTypeConverter()->getIndexType()));
      argsTypes.push_back(mlir::LLVM::LLVMPointerType::get(getTypeConverter()->getIndexType()));
      argsTypes.push_back(op.getAlpha().getType());

      auto functionType = mlir::LLVM::LLVMFunctionType::get(op.getFunctionType().getResult(0), argsTypes);

      auto newOp = rewriter.replaceOpWithNewOp<mlir::LLVM::LLVMFuncOp>(op, op.getSymName(), functionType);
      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      // Map the "time" variable
      mapping.map(op.getTime(), newOp.getArgument(0));

      // The lowered function will receive a pointer to the array of variables.
      assert(!op.getVariables().empty());

      mlir::Value variablesPtr = newOp.getArgument(1);
      variablesPtr = rewriter.create<mlir::LLVM::BitcastOp>(variablesPtr.getLoc(), mlir::LLVM::LLVMPointerType::get(getVoidPtrType()), variablesPtr);

      for (auto variable : llvm::enumerate(op.getVariables())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(variablesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), variable.index()));
        mlir::Value variablePtr = rewriter.create<mlir::LLVM::GEPOp>(variablesPtr.getLoc(), variablesPtr.getType(), variablesPtr, index);
        mlir::Value mappedVariable = rewriter.create<mlir::LLVM::LoadOp>(variablePtr.getLoc(), variablePtr);
        mappedVariable = rewriter.create<mlir::LLVM::BitcastOp>(mappedVariable.getLoc(), mlir::LLVM::LLVMPointerType::get(op.getVariables()[variable.index()].getType()), mappedVariable);
        mappedVariable = rewriter.create<mlir::LLVM::LoadOp>(mappedVariable.getLoc(), mappedVariable);

        mapping.map(variable.value(), mappedVariable);
      }

      // The equation indices are also passed through an array
      mlir::Value equationIndicesPtr = newOp.getArgument(2);

      for (auto equationIndex : llvm::enumerate(op.getEquationIndices())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(equationIndicesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), equationIndex.index()));
        mlir::Value equationIndexPtr = rewriter.create<mlir::LLVM::GEPOp>(equationIndicesPtr.getLoc(), equationIndicesPtr.getType(), equationIndicesPtr, index);
        mlir::Value mappedEquationIndex = rewriter.create<mlir::LLVM::LoadOp>(equationIndexPtr.getLoc(), equationIndexPtr);
        mappedEquationIndex = getTypeConverter()->materializeSourceConversion(rewriter, mappedEquationIndex.getLoc(), rewriter.getIndexType(), mappedEquationIndex);
        mapping.map(equationIndex.value(), mappedEquationIndex);
      }

      // The variable indices are also passed through an array
      mlir::Value variableIndicesPtr = newOp.getArgument(3);

      for (auto variableIndex : llvm::enumerate(op.getVariableIndices())) {
        mlir::Value index = rewriter.create<mlir::arith::ConstantOp>(equationIndicesPtr.getLoc(), rewriter.getIntegerAttr(getIndexType(), variableIndex.index()));
        mlir::Value variableIndexPtr = rewriter.create<mlir::LLVM::GEPOp>(variableIndicesPtr.getLoc(), variableIndicesPtr.getType(), variableIndicesPtr, index);
        mlir::Value mappedVariableIndex = rewriter.create<mlir::LLVM::LoadOp>(variableIndexPtr.getLoc(), variableIndexPtr);
        mappedVariableIndex = getTypeConverter()->materializeSourceConversion(rewriter, mappedVariableIndex.getLoc(), rewriter.getIndexType(), mappedVariableIndex);
        mapping.map(variableIndex.value(), mappedVariableIndex);
      }

      // Add the "alpha" variable
      mapping.map(op.getAlpha(), newOp.getArgument(4));

      // Clone the blocks structure
      for (auto& block : llvm::enumerate(op.getBodyRegion().getBlocks())) {
        if (block.index() != 0) {
          std::vector<mlir::Location> argLocations;

          for (const auto& arg : block.value().getArguments()) {
            argLocations.push_back(arg.getLoc());
          }

          mlir::Block* clonedBlock = rewriter.createBlock(
              &newOp.getBody(),
              newOp.getBody().end(),
              block.value().getArgumentTypes(),
              argLocations);

          mapping.map(&block.value(), clonedBlock);

          for (const auto& [original, cloned] : llvm::zip(block.value().getArguments(), clonedBlock->getArguments())) {
            mapping.map(original, cloned);
          }
        }
      }

      // Clone the original operations
      for (auto& block : llvm::enumerate(op.getBodyRegion())) {
        if (block.index() == 0) {
          rewriter.setInsertionPointToEnd(entryBlock);
        } else {
          rewriter.setInsertionPointToStart(mapping.lookup(&block.value()));
        }

        for (auto& bodyOp : block.value().getOperations()) {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };

  struct ReturnOpLowering : public KINSOLOpConversion<ReturnOp>
  {
    using KINSOLOpConversion<ReturnOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(ReturnOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      rewriter.replaceOpWithNewOp<mlir::LLVM::ReturnOp>(op, adaptor.getOperands());
      return mlir::success();
    }
  };

  struct AddResidualOpLowering : public KINSOLOpConversion<AddResidualOp>
  {
    using KINSOLOpConversion<AddResidualOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(AddResidualOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 3> newOperands;
      llvm::SmallVector<std::string, 3> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Equation
      newOperands.push_back(adaptor.getEquation());
      mangledArgsTypes.push_back(mangling.getIntegerType(adaptor.getEquation().getType().getIntOrFloatBitWidth()));

      // Residual function address
      auto function = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(op.getFunction());
      mlir::Value functionAddress = rewriter.create<mlir::LLVM::AddressOfOp>(loc, function);
      functionAddress = rewriter.create<mlir::LLVM::BitcastOp>(loc, getVoidPtrType(), functionAddress);

      newOperands.push_back(functionAddress);
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the call to the runtime library
      auto resultType = getVoidType();
      auto mangledResultType = mangling.getVoidType();
      auto functionName = mangling.getMangledFunction("kinsolAddResidual", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, newOperands);

      return mlir::success();
    }
  };

  struct AddJacobianOpLowering : public KINSOLOpConversion<AddJacobianOp>
  {
    using KINSOLOpConversion<AddJacobianOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(AddJacobianOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 4> newOperands;
      llvm::SmallVector<std::string, 4> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Equation
      newOperands.push_back(adaptor.getEquation());
      mangledArgsTypes.push_back(mangling.getIntegerType(adaptor.getEquation().getType().getIntOrFloatBitWidth()));

      // Variable
      newOperands.push_back(adaptor.getVariable());
      mangledArgsTypes.push_back(mangling.getIntegerType(adaptor.getVariable().getType().getIntOrFloatBitWidth()));

      // Jacobian function address
      auto function = module.lookupSymbol<mlir::LLVM::LLVMFuncOp>(op.getFunction());
      mlir::Value functionAddress = rewriter.create<mlir::LLVM::AddressOfOp>(loc, function);
      functionAddress = rewriter.create<mlir::LLVM::BitcastOp>(loc, getVoidPtrType(), functionAddress);

      newOperands.push_back(functionAddress);
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the call to the runtime library
      auto resultType = getVoidType();
      auto mangledResultType = mangling.getVoidType();
      auto functionName = mangling.getMangledFunction("kinsolAddJacobian", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, newOperands);

      return mlir::success();
    }
  };

  struct InitOpLowering : public KINSOLOpConversion<InitOp>
  {
    using KINSOLOpConversion<InitOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(InitOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 3> newOperands;
      llvm::SmallVector<std::string, 1> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the call to the runtime library
      auto resultType = getVoidType();
      auto mangledResultType = mangling.getVoidType();
      auto functionName = mangling.getMangledFunction("kinsolInit", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, adaptor.getOperands());

      return mlir::success();
    }
  };

  struct FreeOpLowering : public KINSOLOpConversion<FreeOp>
  {
    using KINSOLOpConversion<FreeOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(FreeOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 1> newOperands;
      llvm::SmallVector<std::string, 1> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the call to the runtime library
      auto resultType = getVoidType();
      auto mangledResultType = mangling.getVoidType();
      auto functionName = mangling.getMangledFunction("kinsolFree", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, newOperands);
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, adaptor.getOperands());

      return mlir::success();
    }
  };

  struct PrintStatisticsOpLowering : public KINSOLOpConversion<PrintStatisticsOp>
  {
    using KINSOLOpConversion<PrintStatisticsOp>::KINSOLOpConversion;

    mlir::LogicalResult matchAndRewrite(PrintStatisticsOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto loc = op.getLoc();
      auto module = op->getParentOfType<mlir::ModuleOp>();

      RuntimeFunctionsMangling mangling;

      llvm::SmallVector<mlir::Value, 1> newOperands;
      llvm::SmallVector<std::string, 1> mangledArgsTypes;

      // KINSOL instance
      newOperands.push_back(adaptor.getInstance());
      mangledArgsTypes.push_back(mangling.getVoidPointerType());

      // Create the call to the runtime library
      auto resultType = getVoidType();
      auto mangledResultType = mangling.getVoidType();
      auto functionName = mangling.getMangledFunction("kinsolPrintStatistics", mangledResultType, mangledArgsTypes);
      auto callee = getOrDeclareFunction(rewriter, module, loc, functionName, resultType, adaptor.getOperands());
      rewriter.replaceOpWithNewOp<mlir::LLVM::CallOp>(op, callee, adaptor.getOperands());

      return mlir::success();
    }
  };
}

static void populateKINSOLFunctionLikeOpsConversionPatterns(
    mlir::RewritePatternSet& patterns,
    mlir::kinsol::LLVMTypeConverter& typeConverter,
    unsigned int bitWidth)
{
  patterns.insert<
      VariableGetterOpLowering,
      VariableSetterOpLowering,
      ResidualFunctionOpLowering,
      JacobianFunctionOpLowering,
      ReturnOpLowering>(typeConverter, bitWidth);
}

static void populateKINSOLConversionPatterns(
    mlir::RewritePatternSet& patterns,
    mlir::kinsol::LLVMTypeConverter& typeConverter,
    unsigned int bitWidth)
{
  patterns.insert<
      CreateOpLowering,
      AddEquationOpLowering,
      AddAlgebraicVariableOpLowering,
      AddVariableAccessOpLowering,
      AddResidualOpLowering,
      AddJacobianOpLowering,
      InitOpLowering,
      FreeOpLowering,
      PrintStatisticsOpLowering>(typeConverter, bitWidth);
}

namespace marco::codegen
{
  class KINSOLToLLVMConversionPass : public mlir::impl::KINSOLToLLVMConversionPassBase<KINSOLToLLVMConversionPass>
  {
    public:
      using KINSOLToLLVMConversionPassBase::KINSOLToLLVMConversionPassBase;

      void runOnOperation() override
      {
        // Convert the function-like operations.
        // This must be done first and as an independent step, as other operations within
        // the KINSOL dialect will need the addresses of such functions when being converted.

        if (mlir::failed(convertFunctionsLikeOps())) {
          mlir::emitError(getOperation().getLoc(), "Error in converting the KINSOL function-like operations");
          return signalPassFailure();
        }

        // Convert the rest of the KINSOL dialect
        if (mlir::failed(convertOperations())) {
          mlir::emitError(getOperation().getLoc(), "Error in converting the KINSOL operations");
          return signalPassFailure();
        }
      }

    private:
      mlir::LogicalResult convertFunctionsLikeOps()
      {
        auto module = getOperation();
        mlir::ConversionTarget target(getContext());

        target.addLegalDialect<mlir::LLVM::LLVMDialect>();

        target.addIllegalOp<
            VariableGetterOp,
            VariableSetterOp,
            ResidualFunctionOp,
            JacobianFunctionOp,
            ReturnOp>();

        target.markUnknownOpDynamicallyLegal([](mlir::Operation* op) {
          return true;
        });

        mlir::LowerToLLVMOptions llvmLoweringOptions(&getContext());
        LLVMTypeConverter typeConverter(&getContext(), llvmLoweringOptions);

        mlir::RewritePatternSet patterns(&getContext());
        populateKINSOLFunctionLikeOpsConversionPatterns(patterns, typeConverter, 64);

        if (auto status = applyPartialConversion(module, target, std::move(patterns)); mlir::failed(status)) {
          return status;
        }

        return mlir::success();
      }

      mlir::LogicalResult convertOperations()
      {
        auto module = getOperation();
        mlir::ConversionTarget target(getContext());

        target.addIllegalDialect<mlir::kinsol::KINSOLDialect>();
        target.addLegalDialect<mlir::LLVM::LLVMDialect>();

        target.markUnknownOpDynamicallyLegal([](mlir::Operation* op) {
          return true;
        });

        mlir::LowerToLLVMOptions llvmLoweringOptions(&getContext());
        llvmLoweringOptions.dataLayout.reset(dataLayout);

        LLVMTypeConverter typeConverter(&getContext(), llvmLoweringOptions);

        mlir::RewritePatternSet patterns(&getContext());
        populateKINSOLConversionPatterns(patterns, typeConverter, 64);

        if (auto status = applyPartialConversion(module, target, std::move(patterns)); mlir::failed(status)) {
          return status;
        }

        return mlir::success();
      }
  };
}

namespace
{
  struct AddAlgebraicVariableOpTypes : public mlir::OpConversionPattern<AddAlgebraicVariableOp>
  {
    using mlir::OpConversionPattern<AddAlgebraicVariableOp>::OpConversionPattern;

    mlir::LogicalResult matchAndRewrite(AddAlgebraicVariableOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      auto newOp = mlir::cast<AddAlgebraicVariableOp>(rewriter.cloneWithoutRegions(*op.getOperation()));
      newOp->setOperands(adaptor.getOperands());

      for (auto result : newOp->getResults()) {
        result.setType(getTypeConverter()->convertType(result.getType()));
      }

      rewriter.replaceOp(op, newOp->getResults());
      return mlir::success();
    }
  };

  struct VariableGetterOpTypes : public mlir::OpConversionPattern<VariableGetterOp>
  {
    using mlir::OpConversionPattern<VariableGetterOp>::OpConversionPattern;

    mlir::LogicalResult matchAndRewrite(VariableGetterOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      mlir::Type resultType = getTypeConverter()->convertType(op.getFunctionType().getResult(0));
      mlir::Type variableType = getTypeConverter()->convertType(op.getVariable().getType());

      auto newOp = rewriter.replaceOpWithNewOp<VariableGetterOp>(
          op, op.getSymName(), resultType, variableType, op.getVariableIndices().size());

      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      for (const auto& [original, cloned] : llvm::zip(op.getFunctionBody().getArguments(), newOp.getFunctionBody().getArguments())) {
        mapping.map(original, cloned);
      }

      auto castArgFn = [&](mlir::Value originalArg) {
        auto castedClonedArg = getTypeConverter()->materializeSourceConversion(
            rewriter, originalArg.getLoc(), originalArg.getType(), mapping.lookup(originalArg));

        mapping.map(originalArg, castedClonedArg);
      };

      castArgFn(op.getVariable());
      assert(op.getFunctionBody().getBlocks().size() == 1);

      for (auto& bodyOp : op.getFunctionBody().getOps()) {
        if (auto returnOp = mlir::dyn_cast<ReturnOp>(bodyOp)) {
          std::vector<mlir::Value> returnValues;

          for (auto returnValue : returnOp.operands()) {
            returnValues.push_back(getTypeConverter()->materializeTargetConversion(
                rewriter, returnOp.getLoc(),
                getTypeConverter()->convertType(returnValue.getType()),
                mapping.lookup(returnValue)));
          }

          rewriter.create<ReturnOp>(returnOp.getLoc(), returnValues);
        } else {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };

  struct VariableSetterOpTypes : public mlir::OpConversionPattern<VariableSetterOp>
  {
    using mlir::OpConversionPattern<VariableSetterOp>::OpConversionPattern;

    mlir::LogicalResult matchAndRewrite(VariableSetterOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      mlir::Type variableType = getTypeConverter()->convertType(op.getVariable().getType());
      mlir::Type valueType = getTypeConverter()->convertType(op.getValue().getType());

      auto newOp = rewriter.replaceOpWithNewOp<VariableSetterOp>(
          op, op.getSymName(), variableType, valueType, op.getVariableIndices().size());

      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      for (const auto& [original, cloned] : llvm::zip(op.getFunctionBody().getArguments(), newOp.getFunctionBody().getArguments())) {
        mapping.map(original, cloned);
      }

      auto castArgFn = [&](mlir::Value originalArg) {
        auto castedClonedArg = getTypeConverter()->materializeSourceConversion(
            rewriter, originalArg.getLoc(), originalArg.getType(), mapping.lookup(originalArg));

        mapping.map(originalArg, castedClonedArg);
      };

      castArgFn(op.getVariable());
      castArgFn(op.getValue());

      assert(op.getFunctionBody().getBlocks().size() == 1);

      for (auto& bodyOp : op.getFunctionBody().getOps()) {
        rewriter.clone(bodyOp, mapping);
      }

      return mlir::success();
    }
  };

  struct ResidualFunctionOpTypes : public mlir::OpConversionPattern<ResidualFunctionOp>
  {
    using mlir::OpConversionPattern<ResidualFunctionOp>::OpConversionPattern;

    mlir::LogicalResult matchAndRewrite(ResidualFunctionOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      mlir::Type timeType = getTypeConverter()->convertType(op.getTime().getType());

      llvm::SmallVector<mlir::Type> variablesTypes;

      for (auto variable : op.getVariables()) {
        variablesTypes.push_back(getTypeConverter()->convertType(variable.getType()));
      }

      mlir::Type differenceType = getTypeConverter()->convertType(op.getFunctionType().getResult(0));

      auto newOp = rewriter.replaceOpWithNewOp<ResidualFunctionOp>(
          op, op.getSymName(), timeType, variablesTypes, op.getEquationRank().getSExtValue(), differenceType);

      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      for (const auto& [original, cloned] : llvm::zip(op.getFunctionBody().getArguments(), newOp.getFunctionBody().getArguments())) {
        mapping.map(original, cloned);
      }

      auto castArgFn = [&](mlir::Value originalArg) {
        auto castedClonedArg = getTypeConverter()->materializeSourceConversion(
            rewriter, originalArg.getLoc(), originalArg.getType(), mapping.lookup(originalArg));

        mapping.map(originalArg, castedClonedArg);
      };

      castArgFn(op.getTime());

      for (auto variable : op.getVariables()) {
        castArgFn(variable);
      }

      assert(op.getFunctionBody().getBlocks().size() == 1);

      for (auto& bodyOp : op.getFunctionBody().getOps()) {
        if (auto returnOp = mlir::dyn_cast<ReturnOp>(bodyOp)) {
          std::vector<mlir::Value> returnValues;

          for (auto returnValue : returnOp.operands()) {
            returnValues.push_back(getTypeConverter()->materializeTargetConversion(
                rewriter, returnOp.getLoc(), differenceType, mapping.lookup(returnValue)));
          }

          rewriter.create<ReturnOp>(returnOp.getLoc(), returnValues);
        } else {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };

  struct JacobianFunctionOpTypes : public mlir::OpConversionPattern<JacobianFunctionOp>
  {
    using mlir::OpConversionPattern<JacobianFunctionOp>::OpConversionPattern;

    mlir::LogicalResult matchAndRewrite(JacobianFunctionOp op, OpAdaptor adaptor, mlir::ConversionPatternRewriter& rewriter) const override
    {
      mlir::Type timeType = getTypeConverter()->convertType(op.getTime().getType());

      llvm::SmallVector<mlir::Type> variablesTypes;

      for (auto variable : op.getVariables()) {
        variablesTypes.push_back(getTypeConverter()->convertType(variable.getType()));
      }

      mlir::Type alphaType = getTypeConverter()->convertType(op.getAlpha().getType());
      mlir::Type resultType = getTypeConverter()->convertType(op.getFunctionType().getResult(0));

      auto newOp = rewriter.replaceOpWithNewOp<JacobianFunctionOp>(
          op, op.getSymName(), timeType, variablesTypes,
          op.getEquationRank().getSExtValue(),
          op.getVariableRank().getSExtValue(),
          alphaType, resultType);

      mlir::Block* entryBlock = newOp.addEntryBlock();
      rewriter.setInsertionPointToStart(entryBlock);

      mlir::BlockAndValueMapping mapping;

      for (const auto& [original, cloned] : llvm::zip(op.getFunctionBody().getArguments(), newOp.getFunctionBody().getArguments())) {
        mapping.map(original, cloned);
      }

      auto castArgFn = [&](mlir::Value originalArg) {
        auto castedClonedArg = getTypeConverter()->materializeSourceConversion(
            rewriter, originalArg.getLoc(), originalArg.getType(), mapping.lookup(originalArg));

        mapping.map(originalArg, castedClonedArg);
      };

      castArgFn(op.getTime());

      for (auto variable : op.getVariables()) {
        castArgFn(variable);
      }

      castArgFn(op.getAlpha());

      assert(op.getFunctionBody().getBlocks().size() == 1);

      for (auto& bodyOp : op.getFunctionBody().getOps()) {
        if (auto returnOp = mlir::dyn_cast<ReturnOp>(bodyOp)) {
          std::vector<mlir::Value> returnValues;

          for (auto returnValue : returnOp.operands()) {
            returnValues.push_back(getTypeConverter()->materializeTargetConversion(
                rewriter, returnOp.getLoc(), resultType, mapping.lookup(returnValue)));
          }

          rewriter.create<ReturnOp>(returnOp.getLoc(), returnValues);
        } else {
          rewriter.clone(bodyOp, mapping);
        }
      }

      return mlir::success();
    }
  };
}

namespace mlir
{
  void populateKINSOLStructuralTypeConversionsAndLegality(
      mlir::LLVMTypeConverter& typeConverter,
      mlir::RewritePatternSet& patterns,
      mlir::ConversionTarget& target)
  {
    typeConverter.addConversion([&](InstanceType type) {
     return type;
    });

    typeConverter.addConversion([&](EquationType type) {
      return type;
    });

    typeConverter.addConversion([&](VariableType type) {
      return type;
    });

    patterns.add<
        AddAlgebraicVariableOpTypes>(typeConverter, patterns.getContext());

    patterns.add<
        VariableGetterOpTypes,
        VariableSetterOpTypes,
        ResidualFunctionOpTypes,
        JacobianFunctionOpTypes>(typeConverter, patterns.getContext());

    target.addDynamicallyLegalOp<AddAlgebraicVariableOp>([&](mlir::Operation* op) {
      return typeConverter.isLegal(op);
    });

    target.addDynamicallyLegalOp<VariableGetterOp>([&](VariableGetterOp op) {
      if (!typeConverter.isLegal(op.getFunctionType().getResult(0))) {
        return false;
      }

      if (!typeConverter.isLegal(op.getVariable().getType())) {
        return false;
      }

      return true;
    });

    target.addDynamicallyLegalOp<VariableSetterOp>([&](VariableSetterOp op) {
      if (!typeConverter.isLegal(op.getVariable().getType())) {
        return false;
      }

      if (!typeConverter.isLegal(op.getValue().getType())) {
        return false;
      }

      return true;
    });

    target.addDynamicallyLegalOp<ResidualFunctionOp>([&](ResidualFunctionOp op) {
      if (!typeConverter.isLegal(op.getTime().getType())) {
        return false;
      }

      for (auto variable : op.getVariables()) {
        if (!typeConverter.isLegal(variable.getType())) {
          return false;
        }
      }

      if (!typeConverter.isLegal(op.getFunctionType().getResult(0))) {
        return false;
      }

      return true;
    });

    target.addDynamicallyLegalOp<JacobianFunctionOp>([&](JacobianFunctionOp op) {
      if (!typeConverter.isLegal(op.getTime().getType())) {
        return false;
      }

      for (auto variable : op.getVariables()) {
        if (!typeConverter.isLegal(variable.getType())) {
          return false;
        }
      }

      if (!typeConverter.isLegal(op.getAlpha().getType())) {
        return false;
      }

      if (!typeConverter.isLegal(op.getFunctionType().getResult(0))) {
        return false;
      }

      return true;
    });
  }

  std::unique_ptr<mlir::Pass> createKINSOLToLLVMConversionPass()
  {
    return std::make_unique<KINSOLToLLVMConversionPass>();
  }

  std::unique_ptr<mlir::Pass> createKINSOLToLLVMConversionPass(const KINSOLToLLVMConversionPassOptions& options)
  {
    return std::make_unique<KINSOLToLLVMConversionPass>(options);
  }
}
