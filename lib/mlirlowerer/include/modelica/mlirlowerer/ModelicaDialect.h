#pragma once

#include <mlir/IR/Dialect.h>
#include <modelica/mlirlowerer/Ops.h>
#include <modelica/mlirlowerer/Type.h>

namespace modelica
{
	class ModelicaDialect : public mlir::Dialect
	{
		public:
		explicit ModelicaDialect(mlir::MLIRContext* context);

		/**
		 * Get the dialect namespace.
		 *
		 * @return dialect namespace
		 */
		static llvm::StringRef getDialectNamespace();

		void printType(mlir::Type type, mlir::DialectAsmPrinter& printer) const override;
	};
}