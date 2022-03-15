#ifndef MARCO_CODEGEN_CONVERSION_MODELICA_TYPECONVERTER_H
#define MARCO_CODEGEN_CONVERSION_MODELICA_TYPECONVERTER_H

#include "marco/Codegen/dialects/modelica/Type.h"
#include "marco/Dialect/IDA/IDADialect.h"
#include "mlir/Conversion/LLVMCommon/TypeConverter.h"
#include "mlir/Dialect/LLVMIR/LLVMDialect.h"
#include "mlir/Dialect/StandardOps/IR/Ops.h"
#include "mlir/IR/BuiltinDialect.h"
#include "mlir/IR/MLIRContext.h"
#include "mlir/Transforms/DialectConversion.h"

namespace mlir::modelica
{
  // We inherit from the LLVMTypeConverter in order to retrieve the converted MLIR index type.
	class TypeConverter : public mlir::LLVMTypeConverter
  {
		public:
      TypeConverter(mlir::MLIRContext* context, mlir::LowerToLLVMOptions options, unsigned int bitWidth);

      mlir::Type convertBooleanType(marco::codegen::modelica::BooleanType type);
      mlir::Type convertIntegerType(marco::codegen::modelica::IntegerType type);
      mlir::Type convertRealType(marco::codegen::modelica::RealType type);
      mlir::Type convertArrayType(marco::codegen::modelica::ArrayType type);
      mlir::Type convertUnsizedArrayType(marco::codegen::modelica::UnsizedArrayType type);

      llvm::Optional<mlir::Value> integerTypeTargetMaterialization(
          mlir::OpBuilder& builder, mlir::IntegerType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> floatTypeTargetMaterialization(
          mlir::OpBuilder& builder, mlir::FloatType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> llvmStructTypeTargetMaterialization(
          mlir::OpBuilder& builder, mlir::LLVM::LLVMStructType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> booleanTypeSourceMaterialization(
          mlir::OpBuilder& builder, marco::codegen::modelica::BooleanType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> integerTypeSourceMaterialization(
          mlir::OpBuilder& builder, marco::codegen::modelica::IntegerType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> realTypeSourceMaterialization(
          mlir::OpBuilder& builder, marco::codegen::modelica::RealType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> arrayTypeSourceMaterialization(
          mlir::OpBuilder& builder, marco::codegen::modelica::ArrayType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

      llvm::Optional<mlir::Value> unsizedArrayTypeSourceMaterialization(
          mlir::OpBuilder& builder, marco::codegen::modelica::UnsizedArrayType resultType, mlir::ValueRange inputs, mlir::Location loc) const;

    private:
      llvm::SmallVector<mlir::Type, 3> getArrayDescriptorFields(marco::codegen::modelica::ArrayType type);
      llvm::SmallVector<mlir::Type, 3> getUnsizedArrayDescriptorFields(marco::codegen::modelica::UnsizedArrayType type);

    private:
		  unsigned int bitWidth;
	};
}

#endif // MARCO_CODEGEN_CONVERSION_MODELICA_TYPECONVERTER_H
