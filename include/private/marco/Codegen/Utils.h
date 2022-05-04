#ifndef MARCO_CODEGEN_UTILS_H
#define MARCO_CODEGEN_UTILS_H

#include "mlir/IR/Builders.h"
#include "mlir/IR/BuiltinOps.h"
#include <functional>

namespace marco::codegen
{
  mlir::Type getMostGenericType(mlir::Value x, mlir::Value y);
  mlir::Type getMostGenericType(mlir::Type x, mlir::Type y);

  std::string getUniqueSymbolName(mlir::ModuleOp module, std::function<std::string()> tryFn);

  void copyArray(mlir::OpBuilder& builder, mlir::Location loc, mlir::Value source, mlir::Value destination);
}

#endif // MARCO_CODEGEN_UTILS_H
