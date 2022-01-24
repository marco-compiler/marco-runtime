#pragma once

#include "mlir/Pass/Pass.h"

namespace marco::codegen
{
	struct FunctionsVectorizationOptions
	{
		bool assertions = true;

		static const FunctionsVectorizationOptions& getDefaultOptions() {
			static FunctionsVectorizationOptions options;
			return options;
		}
	};

	std::unique_ptr<mlir::Pass> createFunctionsVectorizationPass(
			FunctionsVectorizationOptions options = FunctionsVectorizationOptions::getDefaultOptions());

	inline void registerFunctionsVectorizationPass()
	{
		mlir::registerPass("vectorize-functions", "Convert vectorized functions in loops with scalar calls",
											 []() -> std::unique_ptr<::mlir::Pass> {
												 return createFunctionsVectorizationPass();
											 });
	}
}
