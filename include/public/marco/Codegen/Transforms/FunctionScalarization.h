#ifndef MARCO_CODEN_TRANSFORMS_FUNCTIONSCALARIZATION_H
#define MARCO_CODEN_TRANSFORMS_FUNCTIONSCALARIZATION_H

#include "mlir/Pass/Pass.h"

namespace marco::codegen
{
	struct FunctionScalarizationOptions
	{
		bool assertions = true;

		static const FunctionScalarizationOptions& getDefaultOptions() {
			static FunctionScalarizationOptions options;
			return options;
		}
	};

	std::unique_ptr<mlir::Pass> createFunctionScalarizationPass(
      FunctionScalarizationOptions options = FunctionScalarizationOptions::getDefaultOptions());

	inline void registerFunctionScalarizationPass()
	{
		mlir::registerPass(
        "scalarize", "Convert vectorized functions in loops with scalar calls",
        []() -> std::unique_ptr<::mlir::Pass> {
          return createFunctionScalarizationPass();
        });
	}
}

#endif // MARCO_CODEN_TRANSFORMS_FUNCTIONSCALARIZATION_H
