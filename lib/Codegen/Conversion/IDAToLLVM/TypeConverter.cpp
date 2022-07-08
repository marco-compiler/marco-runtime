#include "marco/Codegen/Conversion/IDAToLLVM/TypeConverter.h"

namespace mlir::ida
{
  TypeConverter::TypeConverter(mlir::MLIRContext* context, mlir::LowerToLLVMOptions options)
      : mlir::LLVMTypeConverter(context, options)
  {
    addConversion([&](InstanceType type) {
      return convertInstanceType(type);
    });

    addConversion([&](VariableType type) {
      return convertVariableType(type);
    });

    addConversion([&](EquationType type) {
      return convertEquationType(type);
    });

    addTargetMaterialization([&](mlir::OpBuilder& builder, mlir::LLVM::LLVMPointerType resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
      return opaquePointerTypeTargetMaterialization(builder, resultType, inputs, loc);
    });

    addTargetMaterialization([&](mlir::OpBuilder& builder, mlir::IntegerType resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
      return integerTypeTargetMaterialization(builder, resultType, inputs, loc);
    });

    addSourceMaterialization([&](mlir::OpBuilder& builder, InstanceType resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
      return instanceTypeSourceMaterialization(builder, resultType, inputs, loc);
    });

    addSourceMaterialization([&](mlir::OpBuilder& builder, VariableType resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
      return variableTypeSourceMaterialization(builder, resultType, inputs, loc);
    });

    addSourceMaterialization([&](mlir::OpBuilder& builder, EquationType resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
      return equationTypeSourceMaterialization(builder, resultType, inputs, loc);
    });
  }

  mlir::Type TypeConverter::convertInstanceType(InstanceType type)
  {
    auto integerType = convertType(mlir::IntegerType::get(type.getContext(), 8));
    return mlir::LLVM::LLVMPointerType::get(integerType);
  }

  mlir::Type TypeConverter::convertVariableType(VariableType type)
  {
    return mlir::IntegerType::get(type.getContext(), 64);
  }

  mlir::Type TypeConverter::convertEquationType(EquationType type)
  {
    return mlir::IntegerType::get(type.getContext(), 64);
  }

  llvm::Optional<mlir::Value> TypeConverter::opaquePointerTypeTargetMaterialization(
      mlir::OpBuilder& builder, mlir::LLVM::LLVMPointerType resultType, mlir::ValueRange inputs, mlir::Location loc) const
  {
    if (inputs.size() != 1) {
      return llvm::None;
    }

    if (!inputs[0].getType().isa<InstanceType>()) {
      return llvm::None;
    }

    auto elementType = resultType.getElementType().dyn_cast<mlir::IntegerType>();

    if (!elementType || elementType.getIntOrFloatBitWidth() != 8) {
      return llvm::None;
    }

    return builder.create<mlir::UnrealizedConversionCastOp>(loc, resultType, inputs[0]).getResult(0);
  }

  llvm::Optional<mlir::Value> TypeConverter::integerTypeTargetMaterialization(
      mlir::OpBuilder& builder, mlir::IntegerType resultType, mlir::ValueRange inputs, mlir::Location loc) const
  {
    if (inputs.size() != 1) {
      return llvm::None;
    }

    if (!inputs[0].getType().isa<VariableType, EquationType>()) {
      return llvm::None;
    }

    return builder.create<mlir::UnrealizedConversionCastOp>(loc, resultType, inputs[0]).getResult(0);
  }

  llvm::Optional<mlir::Value> TypeConverter::instanceTypeSourceMaterialization(
      mlir::OpBuilder& builder, InstanceType resultType, mlir::ValueRange inputs, mlir::Location loc) const
  {
    if (inputs.size() != 1) {
      return llvm::None;
    }

    auto pointerType = inputs[0].getType().dyn_cast<mlir::LLVM::LLVMPointerType>();

    if (!pointerType) {
      return llvm::None;
    }

    auto elementType = pointerType.getElementType().dyn_cast<mlir::IntegerType>();

    if (!elementType || elementType.getIntOrFloatBitWidth() != 8) {
      return llvm::None;
    }

    return builder.create<mlir::UnrealizedConversionCastOp>(loc, resultType, inputs[0]).getResult(0);
  }

  llvm::Optional<mlir::Value> TypeConverter::variableTypeSourceMaterialization(
      mlir::OpBuilder& builder, VariableType resultType, mlir::ValueRange inputs, mlir::Location loc) const
  {
    if (inputs.size() != 1) {
      return llvm::None;
    }

    if (!inputs[0].getType().isa<mlir::IntegerType>()) {
      return llvm::None;
    }

    return builder.create<mlir::UnrealizedConversionCastOp>(loc, resultType, inputs[0]).getResult(0);
  }

  llvm::Optional<mlir::Value> TypeConverter::equationTypeSourceMaterialization(
      mlir::OpBuilder& builder, EquationType resultType, mlir::ValueRange inputs, mlir::Location loc) const
  {
    if (inputs.size() != 1) {
      return llvm::None;
    }

    if (!inputs[0].getType().isa<mlir::IntegerType>()) {
      return llvm::None;
    }

    return builder.create<mlir::UnrealizedConversionCastOp>(loc, resultType, inputs[0]).getResult(0);
  }
}