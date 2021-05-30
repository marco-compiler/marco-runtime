#pragma once

// Just a convenience header file to include the Modelica passes

#include "passes/AutomaticDifferentiation.h"
#include "passes/BufferDeallocation.h"
#include "passes/ExplicitCastInsertion.h"
#include "passes/LowerToLLVM.h"
#include "passes/ModelicaConversion.h"
#include "passes/ResultBuffersToArgs.h"
#include "passes/SolveModel.h"

namespace modelica::codegen
{
	inline void registerModelicaPasses() {
		registerAutomaticDifferentiationPass();
		registerBufferDeallocationPass();
		registerExplicitCastInsertionPass();
		registerLLVMLoweringPass();
		registerFunctionConversionPass();
		registerModelicaConversionPass();
		registerLowerToCFGPass();
		registerResultBuffersToArgsPass();
		registerSolveModelPass();
	}
}
