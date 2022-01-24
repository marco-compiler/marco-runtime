#pragma once

#include "marco/codegen/dialects/modelica/Attribute.h"
#include "marco/codegen/dialects/modelica/Ops.h"
#include "marco/codegen/dialects/modelica/Type.h"
#include "mlir/IR/Dialect.h"

namespace marco::codegen::modelica
{
	class ModelicaDialect : public mlir::Dialect
	{
		public:
		explicit ModelicaDialect(mlir::MLIRContext* context);

		static llvm::StringRef getDialectNamespace();

		mlir::Type parseType(mlir::DialectAsmParser& parser) const override;
		void printType(mlir::Type type, mlir::DialectAsmPrinter& printer) const override;

		mlir::Attribute parseAttribute(mlir::DialectAsmParser& parser, mlir::Type type) const override;
		void printAttribute(mlir::Attribute attribute, mlir::DialectAsmPrinter& printer) const override;

		mlir::Operation* materializeConstant(mlir::OpBuilder& builder, mlir::Attribute value, mlir::Type type, mlir::Location loc) override;
	};
}
