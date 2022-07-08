#ifndef MARCO_CODEGEN_CONVERSION_MODELICA_LOWERTOCFG_H
#define MARCO_CODEGEN_CONVERSION_MODELICA_LOWERTOCFG_H

#include "mlir/Pass/Pass.h"

namespace marco::codegen
{
  struct ModelicaToCFOptions
  {
    unsigned int bitWidth = 64;
    bool outputArraysPromotion = true;
    bool inlining = true;

    static const ModelicaToCFOptions& getDefaultOptions();
  };

	std::unique_ptr<mlir::Pass> createModelicaToCFPass(
      ModelicaToCFOptions options = ModelicaToCFOptions::getDefaultOptions());
}

#endif // MARCO_CODEGEN_CONVERSION_MODELICA_LOWERTOCFG_H