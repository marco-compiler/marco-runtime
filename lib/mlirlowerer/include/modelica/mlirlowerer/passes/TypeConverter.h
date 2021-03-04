#pragma once

#include <mlir/Conversion/StandardToLLVM/ConvertStandardToLLVM.h>
#include <mlir/Dialect/LLVMIR/LLVMDialect.h>
#include <mlir/IR/MLIRContext.h>
#include <modelica/mlirlowerer/Type.h>

namespace modelica
{
	class TypeConverter : public mlir::LLVMTypeConverter
	{
		public:
		TypeConverter(mlir::MLIRContext* context, mlir::LowerToLLVMOptions options) : mlir::LLVMTypeConverter(context, options)
		{
			addConversion([&](BooleanType type) { return convertBooleanType(type); });
			addConversion([&](IntegerType type) { return convertIntegerType(type); });
			addConversion([&](RealType type) { return convertRealType(type); });
			addConversion([&](PointerType type) { return convertPointerType(type); });

			/*
			addSourceMaterialization(
					[&](mlir::OpBuilder &builder, BooleanType resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
						if (inputs.size() != 1)
							return llvm::None;

						if (inputs[0].getType().isa<mlir::IntegerType>() && inputs[0].getType().getIntOrFloatBitWidth() == 1)
						{
							inputs[0].setType(resultType);
							return inputs[0];
						}

						return llvm::None;
			});
			 */

			/*
			addTargetMaterialization(
					[&](mlir::OpBuilder &builder, mlir::Type resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
						if (inputs.size() != 1)
							return llvm::None;

						return inputs[0];
					});

			addSourceMaterialization(
					[&](mlir::OpBuilder &builder, mlir::Type resultType, mlir::ValueRange inputs, mlir::Location loc) -> llvm::Optional<mlir::Value> {
						if (inputs.size() != 1)
							return llvm::None;

						return inputs[0];
					});
					*/
		}

		[[nodiscard]] mlir::Type indexType();

		[[nodiscard]] mlir::Type voidPtrType();

		private:
		mlir::Type convertBooleanType(BooleanType type);
		mlir::Type convertIntegerType(IntegerType type);
		mlir::Type convertRealType(RealType type);
		mlir::Type convertPointerType(PointerType type);

		llvm::SmallVector<mlir::Type, 3> getPointerDescriptorFields(PointerType type);
	};
}
