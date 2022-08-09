#include "marco/Codegen/Conversion/ModelicaCommon/Utils.h"
#include "marco/Dialect/Modelica/ModelicaDialect.h"

using namespace ::marco;
using namespace ::marco::codegen;
using namespace ::mlir::modelica;

namespace marco::codegen
{
  bool isNumeric(mlir::Type type)
  {
    return type.isa<BooleanType, IntegerType, RealType, mlir::IndexType>();
  }

  bool isNumeric(mlir::Value value)
  {
    return isNumeric(value.getType());
  }

  mlir::Type castToMostGenericType(
      mlir::OpBuilder& builder,
      mlir::ValueRange values,
      llvm::SmallVectorImpl<mlir::Value>& castedValues)
  {
    mlir::Type resultType = nullptr;
    mlir::Type resultBaseType = nullptr;

    for (const auto& value : values) {
      mlir::Type type = value.getType();
      mlir::Type baseType = type;

      if (resultType == nullptr) {
        resultType = type;
        resultBaseType = type;

        while (resultBaseType.isa<ArrayType>()) {
          resultBaseType = resultBaseType.cast<ArrayType>().getElementType();
        }

        continue;
      }

      if (type.isa<ArrayType>()) {
        while (baseType.isa<ArrayType>()) {
          baseType = baseType.cast<ArrayType>().getElementType();
        }
      }

      if (resultBaseType.isa<mlir::IndexType>() || baseType.isa<RealType>()) {
        resultType = type;
        resultBaseType = baseType;
      }
    }

    llvm::SmallVector<mlir::Type, 3> types;

    for (const auto& value : values) {
      mlir::Type type = value.getType();

      if (type.isa<ArrayType>()) {
        auto arrayType = type.cast<ArrayType>();
        auto shape = arrayType.getShape();
        types.emplace_back(ArrayType::get(shape, resultBaseType));
      }
      else {
        types.emplace_back(resultBaseType);
      }
    }

    for (const auto& [value, type] : llvm::zip(values, types)) {
      if (value.getType() != type) {
        mlir::Value castedValue = builder.create<CastOp>(value.getLoc(), type, value);
        castedValues.push_back(castedValue);
      } else {
        castedValues.push_back(value);
      }
    }

    return types[0];
  }

  std::vector<mlir::Value> getArrayDynamicDimensions(mlir::OpBuilder& builder, mlir::Location loc, mlir::Value array)
  {
    std::vector<mlir::Value> result;

    assert(array.getType().isa<ArrayType>());
    auto arrayType = array.getType().cast<ArrayType>();

    for (const auto& dimension : llvm::enumerate(arrayType.getShape())) {
      if (dimension.value() == ArrayType::kDynamicSize) {
        mlir::Value dim = builder.create<ConstantOp>(loc, builder.getIndexAttr(dimension.index()));
        result.push_back(builder.create<DimOp>(loc, array, dim));
      }
    }

    return result;
  }
}
