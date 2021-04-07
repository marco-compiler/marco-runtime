#pragma once

#include <memory>
#include <mlir/IR/Value.h>

namespace modelica::codegen::model
{
	class Constant
	{
		public:
		explicit Constant(mlir::Value value);

		[[nodiscard]] mlir::Value getValue() const;

		private:
		mlir::Value value;
	};
}
