#ifndef MARCO_RUNTIME_BUILTINFUNCTIONS_H
#define MARCO_RUNTIME_BUILTINFUNCTIONS_H

#include "marco/Runtime/Support/Mangling.h"
#include "mlir/ExecutionEngine/CRunnerUtils.h"
#include <cstdint>

RUNTIME_FUNC_DECL(abs, bool, bool)
RUNTIME_FUNC_DECL(abs, int32_t , int32_t)
RUNTIME_FUNC_DECL(abs, int64_t, int64_t)
RUNTIME_FUNC_DECL(abs, float, float)
RUNTIME_FUNC_DECL(abs, double, double)

RUNTIME_FUNC_DECL(acos, float, float)
RUNTIME_FUNC_DECL(acos, double, double)

RUNTIME_FUNC_DECL(asin, float, float)
RUNTIME_FUNC_DECL(asin, double, double)

RUNTIME_FUNC_DECL(atan, float, float)
RUNTIME_FUNC_DECL(atan, double, double)

RUNTIME_FUNC_DECL(atan2, float, float, float)
RUNTIME_FUNC_DECL(atan2, double, double, double)

RUNTIME_FUNC_DECL(ceil, bool, bool)
RUNTIME_FUNC_DECL(ceil, int32_t, int32_t)
RUNTIME_FUNC_DECL(ceil, int64_t, int64_t)
RUNTIME_FUNC_DECL(ceil, float, float)
RUNTIME_FUNC_DECL(ceil, double, double)

RUNTIME_FUNC_DECL(cos, float, float)
RUNTIME_FUNC_DECL(cos, double, double)

RUNTIME_FUNC_DECL(cosh, float, float)
RUNTIME_FUNC_DECL(cosh, double, double)

RUNTIME_FUNC_DECL(diagonal, void, ARRAY(bool), ARRAY(bool))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(bool), ARRAY(int32_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(bool), ARRAY(int64_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(bool), ARRAY(float))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(bool), ARRAY(double))

RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int32_t), ARRAY(bool))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int32_t), ARRAY(int32_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int32_t), ARRAY(int64_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int32_t), ARRAY(float))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int32_t), ARRAY(double))

RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int64_t), ARRAY(bool))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int64_t), ARRAY(int32_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int64_t), ARRAY(int64_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int64_t), ARRAY(float))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(int64_t), ARRAY(double))

RUNTIME_FUNC_DECL(diagonal, void, ARRAY(float), ARRAY(bool))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(float), ARRAY(int32_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(float), ARRAY(int64_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(float), ARRAY(float))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(float), ARRAY(double))

RUNTIME_FUNC_DECL(diagonal, void, ARRAY(double), ARRAY(bool))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(double), ARRAY(int32_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(double), ARRAY(int64_t))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(double), ARRAY(float))
RUNTIME_FUNC_DECL(diagonal, void, ARRAY(double), ARRAY(double))

RUNTIME_FUNC_DECL(div, bool, bool, bool)
RUNTIME_FUNC_DECL(div, int32_t, int32_t, int32_t)
RUNTIME_FUNC_DECL(div, int64_t, int64_t, int64_t)
RUNTIME_FUNC_DECL(div, float, float, float)
RUNTIME_FUNC_DECL(div, double, double, double)

RUNTIME_FUNC_DECL(exp, float, float)
RUNTIME_FUNC_DECL(exp, double, double)

RUNTIME_FUNC_DECL(floor, bool, bool)
RUNTIME_FUNC_DECL(floor, int32_t, int32_t)
RUNTIME_FUNC_DECL(floor, int64_t, int64_t)
RUNTIME_FUNC_DECL(floor, float, float)
RUNTIME_FUNC_DECL(floor, double, double)

RUNTIME_FUNC_DECL(identity, void, ARRAY(bool))
RUNTIME_FUNC_DECL(identity, void, ARRAY(int32_t))
RUNTIME_FUNC_DECL(identity, void, ARRAY(int64_t))
RUNTIME_FUNC_DECL(identity, void, ARRAY(float))
RUNTIME_FUNC_DECL(identity, void, ARRAY(double))

RUNTIME_FUNC_DECL(integer, bool, bool)
RUNTIME_FUNC_DECL(integer, int32_t, int32_t)
RUNTIME_FUNC_DECL(integer, int64_t, int64_t)
RUNTIME_FUNC_DECL(integer, float, float)
RUNTIME_FUNC_DECL(integer, double, double)

RUNTIME_FUNC_DECL(linspace, void, ARRAY(bool), float, float)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(bool), double, double)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(int32_t), float, float)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(int32_t), double, double)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(int64_t), float, float)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(int64_t), double, double)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(float), float, float)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(float), double, double)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(double), float, float)
RUNTIME_FUNC_DECL(linspace, void, ARRAY(double), double, double)

RUNTIME_FUNC_DECL(log, float, float)
RUNTIME_FUNC_DECL(log, double, double)

RUNTIME_FUNC_DECL(log10, float, float)
RUNTIME_FUNC_DECL(log10, double, double)

RUNTIME_FUNC_DECL(maxArray, bool, ARRAY(bool))
RUNTIME_FUNC_DECL(maxArray, int32_t, ARRAY(int32_t))
RUNTIME_FUNC_DECL(maxArray, int64_t, ARRAY(int64_t))
RUNTIME_FUNC_DECL(maxArray, float, ARRAY(float))
RUNTIME_FUNC_DECL(maxArray, double, ARRAY(double))

RUNTIME_FUNC_DECL(maxScalars, bool, bool, bool)
RUNTIME_FUNC_DECL(maxScalars, int32_t, int32_t, int32_t)
RUNTIME_FUNC_DECL(maxScalars, int64_t, int64_t, int64_t)
RUNTIME_FUNC_DECL(maxScalars, float, float, float)
RUNTIME_FUNC_DECL(maxScalars, double, double, double)

RUNTIME_FUNC_DECL(minArray, bool, ARRAY(bool))
RUNTIME_FUNC_DECL(minArray, int32_t, ARRAY(int32_t))
RUNTIME_FUNC_DECL(minArray, int64_t, ARRAY(int64_t))
RUNTIME_FUNC_DECL(minArray, float, ARRAY(float))
RUNTIME_FUNC_DECL(minArray, double, ARRAY(double))

RUNTIME_FUNC_DECL(minScalars, bool, bool, bool)
RUNTIME_FUNC_DECL(minScalars, int32_t, int32_t, int32_t)
RUNTIME_FUNC_DECL(minScalars, int64_t, int64_t, int64_t)
RUNTIME_FUNC_DECL(minScalars, float, float, float)
RUNTIME_FUNC_DECL(minScalars, double, double, double)

RUNTIME_FUNC_DECL(mod, bool, bool, bool)
RUNTIME_FUNC_DECL(mod, int32_t, int32_t, int32_t)
RUNTIME_FUNC_DECL(mod, int64_t, int64_t, int64_t)
RUNTIME_FUNC_DECL(mod, float, float, float)
RUNTIME_FUNC_DECL(mod, double, double, double)

RUNTIME_FUNC_DECL(ones, void, ARRAY(bool))
RUNTIME_FUNC_DECL(ones, void, ARRAY(int32_t))
RUNTIME_FUNC_DECL(ones, void, ARRAY(int64_t))
RUNTIME_FUNC_DECL(ones, void, ARRAY(float))
RUNTIME_FUNC_DECL(ones, void, ARRAY(double))

RUNTIME_FUNC_DECL(product, bool, ARRAY(bool))
RUNTIME_FUNC_DECL(product, int32_t, ARRAY(int32_t))
RUNTIME_FUNC_DECL(product, int64_t, ARRAY(int64_t))
RUNTIME_FUNC_DECL(product, float, ARRAY(float))
RUNTIME_FUNC_DECL(product, double, ARRAY(double))

RUNTIME_FUNC_DECL(rem, bool, bool, bool)
RUNTIME_FUNC_DECL(rem, int32_t, int32_t, int32_t)
RUNTIME_FUNC_DECL(rem, int64_t, int64_t, int64_t)
RUNTIME_FUNC_DECL(rem, float, float, float)
RUNTIME_FUNC_DECL(rem, double, double, double)

RUNTIME_FUNC_DECL(sign, int32_t, bool)
RUNTIME_FUNC_DECL(sign, int32_t, int32_t)
RUNTIME_FUNC_DECL(sign, int32_t, int64_t)
RUNTIME_FUNC_DECL(sign, int32_t, float)
RUNTIME_FUNC_DECL(sign, int32_t, double)

RUNTIME_FUNC_DECL(sign, int64_t, bool)
RUNTIME_FUNC_DECL(sign, int64_t, int32_t)
RUNTIME_FUNC_DECL(sign, int64_t, int64_t)
RUNTIME_FUNC_DECL(sign, int64_t, float)
RUNTIME_FUNC_DECL(sign, int64_t, double)

RUNTIME_FUNC_DECL(sin, float, float)
RUNTIME_FUNC_DECL(sin, double, double)

RUNTIME_FUNC_DECL(sinh, float, float)
RUNTIME_FUNC_DECL(sinh, double, double)

RUNTIME_FUNC_DECL(sqrt, float, float)
RUNTIME_FUNC_DECL(sqrt, double, double)

RUNTIME_FUNC_DECL(sum, bool, ARRAY(bool))
RUNTIME_FUNC_DECL(sum, int32_t, ARRAY(int32_t))
RUNTIME_FUNC_DECL(sum, int64_t, ARRAY(int64_t))
RUNTIME_FUNC_DECL(sum, float, ARRAY(float))
RUNTIME_FUNC_DECL(sum, double, ARRAY(double))

RUNTIME_FUNC_DECL(symmetric, void, ARRAY(bool), ARRAY(bool))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(bool), ARRAY(int32_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(bool), ARRAY(int64_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(bool), ARRAY(float))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(bool), ARRAY(double))

RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int32_t), ARRAY(bool))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int32_t), ARRAY(int32_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int32_t), ARRAY(int64_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int32_t), ARRAY(float))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int32_t), ARRAY(double))

RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int64_t), ARRAY(bool))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int64_t), ARRAY(int32_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int64_t), ARRAY(int64_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int64_t), ARRAY(float))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(int64_t), ARRAY(double))

RUNTIME_FUNC_DECL(symmetric, void, ARRAY(float), ARRAY(bool))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(float), ARRAY(int32_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(float), ARRAY(int64_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(float), ARRAY(float))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(float), ARRAY(double))

RUNTIME_FUNC_DECL(symmetric, void, ARRAY(double), ARRAY(bool))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(double), ARRAY(int32_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(double), ARRAY(int64_t))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(double), ARRAY(float))
RUNTIME_FUNC_DECL(symmetric, void, ARRAY(double), ARRAY(double))

RUNTIME_FUNC_DECL(tan, float, float)
RUNTIME_FUNC_DECL(tan, double, double)

RUNTIME_FUNC_DECL(tanh, float, float)
RUNTIME_FUNC_DECL(tanh, double, double)

RUNTIME_FUNC_DECL(transpose, void, ARRAY(bool), ARRAY(bool))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(bool), ARRAY(int32_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(bool), ARRAY(int64_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(bool), ARRAY(float))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(bool), ARRAY(double))

RUNTIME_FUNC_DECL(transpose, void, ARRAY(int32_t), ARRAY(bool))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int32_t), ARRAY(int32_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int32_t), ARRAY(int64_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int32_t), ARRAY(float))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int32_t), ARRAY(double))

RUNTIME_FUNC_DECL(transpose, void, ARRAY(int64_t), ARRAY(bool))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int64_t), ARRAY(int32_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int64_t), ARRAY(int64_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int64_t), ARRAY(float))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(int64_t), ARRAY(double))

RUNTIME_FUNC_DECL(transpose, void, ARRAY(float), ARRAY(bool))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(float), ARRAY(int32_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(float), ARRAY(int64_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(float), ARRAY(float))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(float), ARRAY(double))

RUNTIME_FUNC_DECL(transpose, void, ARRAY(double), ARRAY(bool))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(double), ARRAY(int32_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(double), ARRAY(int64_t))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(double), ARRAY(float))
RUNTIME_FUNC_DECL(transpose, void, ARRAY(double), ARRAY(double))

RUNTIME_FUNC_DECL(zeros, void, ARRAY(bool))
RUNTIME_FUNC_DECL(zeros, void, ARRAY(int32_t))
RUNTIME_FUNC_DECL(zeros, void, ARRAY(int64_t))
RUNTIME_FUNC_DECL(zeros, void, ARRAY(float))
RUNTIME_FUNC_DECL(zeros, void, ARRAY(double))

#endif	// MARCO_RUNTIME_BUILTINFUNCTIONS_H
