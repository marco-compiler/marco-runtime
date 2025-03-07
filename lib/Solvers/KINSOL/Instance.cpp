#ifdef SUNDIALS_ENABLE

#include "marco/Runtime/Solvers/KINSOL/Instance.h"
#include "kinsol/kinsol.h"
#include "marco/Runtime/Simulation/Options.h"
#include "marco/Runtime/Solvers/KINSOL/Options.h"
#include "marco/Runtime/Solvers/KINSOL/Profiler.h"
#include "marco/Runtime/Support/MemoryManagement.h"
#include <algorithm>
#include <atomic>
#include <cassert>
#include <functional>
#include <iomanip>
#include <iostream>
#include <set>

using namespace marco::runtime::sundials;
using namespace marco::runtime::sundials::kinsol;

//===---------------------------------------------------------------------===//
// Solver
//===---------------------------------------------------------------------===//

namespace marco::runtime::sundials::kinsol {
KINSOLInstance::KINSOLInstance() {
  // Initially there is are no variables or equations in the instance.
  variableOffsets.push_back(0);
  equationOffsets.push_back(0);

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Instance created" << std::endl;
  }
}

KINSOLInstance::~KINSOLInstance() {
  if (getNumOfScalarEquations() != 0) {
    N_VDestroy(variablesVector);
    N_VDestroy(tolerancesVector);

    KINFree(&kinsolMemory);
    SUNLinSolFree(linearSolver);
    SUNMatDestroy(sparseMatrix);
  }

#if SUNDIALS_VERSION_MAJOR >= 6
  if (ctx != nullptr) {
    SUNContext_Free(&ctx);
  }
#endif

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Instance destroyed" << std::endl;
  }
}

Variable KINSOLInstance::addVariable(uint64_t rank, const uint64_t *dimensions,
                                     VariableGetter getterFunction,
                                     VariableSetter setterFunction,
                                     const char *name) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Adding algebraic variable";

    if (name != nullptr) {
      std::cerr << " \"" << name << "\"";
    }

    std::cerr << std::endl;
  }

  // Add variable offset and dimensions.
  assert(variableOffsets.size() == variablesDimensions.size() + 1);

  VariableDimensions varDimension(rank);
  uint64_t flatSize = 1;

  for (uint64_t i = 0; i < rank; ++i) {
    flatSize *= dimensions[i];
    varDimension[i] = dimensions[i];
  }

  variablesDimensions.push_back(std::move(varDimension));

  size_t offset = variableOffsets.back();
  variableOffsets.push_back(offset + flatSize);

  // Store the getter and setter functions.
  variableGetters.push_back(getterFunction);
  variableSetters.push_back(setterFunction);

  // Return the index of the variable.
  Variable id = getNumOfArrayVariables() - 1;

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "  - ID: " << id << std::endl;
    std::cerr << "  - Rank: " << rank << std::endl;
    std::cerr << "  - Dimensions: [";

    for (uint64_t i = 0; i < rank; ++i) {
      if (i != 0) {
        std::cerr << ",";
      }

      std::cerr << dimensions[i];
    }

    std::cerr << "]" << std::endl;
    std::cerr << "  - Getter function address: "
              << reinterpret_cast<void *>(getterFunction) << std::endl;
    std::cerr << "  - Setter function address: "
              << reinterpret_cast<void *>(setterFunction) << std::endl;
  }

  return id;
}

Equation KINSOLInstance::addEquation(const int64_t *ranges,
                                     uint64_t equationRank,
                                     const char *stringRepresentation) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Adding equation";

    if (stringRepresentation != nullptr) {
      std::cerr << " \"" << stringRepresentation << "\"";
    }

    std::cerr << std::endl;
  }

  // Add the start and end dimensions of the current equation.
  MultidimensionalRange eqRanges = {};
  uint64_t flatSize = 1;

  for (size_t i = 0, e = equationRank * 2; i < e; i += 2) {
    int64_t begin = ranges[i];
    int64_t end = ranges[i + 1];
    eqRanges.push_back({begin, end});
    flatSize *= end - begin;
  }

  equationRanges.push_back(eqRanges);
  size_t offset = equationOffsets.back();
  equationOffsets.push_back(offset + flatSize);

  // Return the index of the equation.
  Equation id = getNumOfVectorizedEquations() - 1;

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "  - ID: " << id << std::endl;
    std::cerr << "  - Rank: " << equationRank << std::endl;
    std::cerr << "  - Ranges: [";

    for (uint64_t i = 0; i < equationRank; ++i) {
      if (i != 0) {
        std::cerr << ",";
      }

      std::cerr << "[" << ranges[i * 2] << "," << (ranges[i * 2 + 1] - 1)
                << "]";
    }

    std::cerr << "]" << std::endl;
  }

  return id;
}

void KINSOLInstance::addVariableAccess(Equation equation, Variable variable,
                                       AccessFunction accessFunction) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Adding access information" << std::endl;
    std::cerr << "  - Equation: " << equation << std::endl;
    std::cerr << "  - Variable: " << variable << std::endl;
    std::cerr << "  - Access function address: "
              << reinterpret_cast<void *>(accessFunction) << std::endl;
  }

  assert(equation < getNumOfVectorizedEquations());
  assert(variable < getNumOfArrayVariables());

  precomputedAccesses = true;

  if (variableAccesses.size() <= (size_t)equation) {
    variableAccesses.resize(equation + 1);
  }

  auto &varAccessList = variableAccesses[equation];
  varAccessList.emplace_back(variable, accessFunction);
}

void KINSOLInstance::setResidualFunction(Equation equation,
                                         ResidualFunction residualFunction) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Setting residual function for equation " << equation
              << ". Address: " << reinterpret_cast<void *>(residualFunction)
              << std::endl;
  }

  if (residualFunctions.size() <= equation) {
    residualFunctions.resize(equation + 1, nullptr);
  }

  residualFunctions[equation] = residualFunction;
}

void KINSOLInstance::addJacobianFunction(Equation equation, Variable variable,
                                         JacobianFunction jacobianFunction,
                                         uint64_t numOfSeeds,
                                         uint64_t *seedSizes) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Setting jacobian function for equation " << equation
              << " and variable " << variable << "\n"
              << "  - Address: " << reinterpret_cast<void *>(jacobianFunction)
              << "\n"
              << "  - Number of seeds : " << numOfSeeds;

    if (numOfSeeds != 0) {
      std::cerr << "\n" << "  - Seed sizes: ";

      for (uint64_t i = 0; i < numOfSeeds; ++i) {
        if (i != 0) {
          std::cerr << ", ";
        }

        std::cerr << seedSizes[i];
      }
    }

    std::cerr << std::endl;
  }

  if (jacobianFunctions.size() <= equation) {
    jacobianFunctions.resize(equation + 1, {});
  }

  if (jacobianFunctions[equation].size() <= variable) {
    jacobianFunctions[equation].resize(
        variable + 1, std::make_pair(nullptr, std::vector<uint64_t>{}));
  }

  jacobianFunctions[equation][variable].first = jacobianFunction;
  jacobianFunctions[equation][variable].second.resize(numOfSeeds);

  for (uint64_t i = 0; i < numOfSeeds; ++i) {
    assert(seedSizes[i] != 0);
    jacobianFunctions[equation][variable].second[i] = seedSizes[i];
  }
}

bool KINSOLInstance::initialize() {
  assert(!initialized && "KINSOL instance has already been initialized");

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Performing initialization" << std::endl;
  }

  memoryPoolId = MemoryPoolManager::getInstance().create();

  // Compute the number of scalar variables.
  scalarVariablesNumber = 0;

  for (Variable var = 0, e = getNumOfArrayVariables(); var < e; ++var) {
    scalarVariablesNumber += getVariableFlatSize(var);
  }

  // Compute the number of scalar equations.
  scalarEquationsNumber = 0;

  for (Equation eq = 0, e = getNumOfVectorizedEquations(); eq < e; ++eq) {
    scalarEquationsNumber += getEquationFlatSize(eq);
  }

  assert(getNumOfScalarVariables() == getNumOfScalarEquations() &&
         "Unbalanced system");

  if (scalarEquationsNumber == 0) {
    // KINSOL has nothing to solve.
    initialized = true;
    return true;
  }

#if SUNDIALS_VERSION_MAJOR >= 6
  // Create the SUNDIALS context.
  if (SUNContext_Create(nullptr, &ctx) != 0) {
    return false;
  }
#endif

  // Create and initialize the variables vector.
#if SUNDIALS_VERSION_MAJOR >= 6
  variablesVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber), ctx);
#else
  variablesVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber));
#endif

  assert(
      checkAllocation(static_cast<void *>(variablesVector), "N_VNew_Serial"));

  for (uint64_t i = 0; i < scalarVariablesNumber; ++i) {
    N_VGetArrayPointer(variablesVector)[i] = 0;
  }

  // Create and initialize the tolerances vector.
#if SUNDIALS_VERSION_MAJOR >= 6
  tolerancesVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber), ctx);
#else
  tolerancesVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber));
#endif

  assert(
      checkAllocation(static_cast<void *>(tolerancesVector), "N_VNew_Serial"));

  for (Variable var = 0; var < getNumOfArrayVariables(); ++var) {
    uint64_t arrayOffset = variableOffsets[var];
    uint64_t flatSize = getVariableFlatSize(var);

    for (uint64_t scalarOffset = 0; scalarOffset < flatSize; ++scalarOffset) {
      uint64_t offset = arrayOffset + scalarOffset;

      N_VGetArrayPointer(tolerancesVector)[offset] =
          std::min(getOptions().maxAlgebraicAbsoluteTolerance,
                   getOptions().absoluteTolerance);
    }
  }

#if SUNDIALS_VERSION_MAJOR >= 6
  variableScaleVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber), ctx);
#else
  variableScaleVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber));
#endif

#if SUNDIALS_VERSION_MAJOR >= 6
  residualScaleVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber), ctx);
#else
  residualScaleVector =
      N_VNew_Serial(static_cast<sunindextype>(scalarVariablesNumber));
#endif

  for (uint64_t i = 0; i < scalarVariablesNumber; ++i) {
    N_VGetArrayPointer(variableScaleVector)[i] = 1;
    N_VGetArrayPointer(residualScaleVector)[i] = 1;
  }

  // Check that all the residual functions have been set.
  assert(residualFunctions.size() == getNumOfVectorizedEquations());

  assert(std::all_of(
      residualFunctions.begin(), residualFunctions.end(),
      [](const ResidualFunction &function) { return function != nullptr; }));

  // Check if the KINSOL instance is not informed about the accesses that all
  // the jacobian functions have been set.
  assert(precomputedAccesses ||
         jacobianFunctions.size() == getNumOfVectorizedEquations());

  assert(precomputedAccesses ||
         std::all_of(jacobianFunctions.begin(), jacobianFunctions.end(),
                     [&](std::vector<JacobianFunctionDescriptor> functions) {
                       if (functions.size() != variableGetters.size()) {
                         return false;
                       }

                       return std::all_of(
                           functions.begin(), functions.end(),
                           [](const JacobianFunctionDescriptor &function) {
                             return function.first != nullptr;
                           });
                     }));

  // Check that all the getters and setters have been set.
  assert(
      std::none_of(variableGetters.begin(), variableGetters.end(),
                   [](VariableGetter getter) { return getter == nullptr; }) &&
      "Not all the variable getters have been set");

  assert(
      std::none_of(variableSetters.begin(), variableSetters.end(),
                   [](VariableSetter setter) { return setter == nullptr; }) &&
      "Not all the variable setters have been set");

  // Reserve the space for data of the jacobian matrix.
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Reserving space for the data of the Jacobian matrix"
              << std::endl;
  }

  jacobianMatrixData.resize(scalarEquationsNumber);

  uint64_t numOfVectorizedEquations = getNumOfVectorizedEquations();

  for (Equation eq = 0; eq < numOfVectorizedEquations; ++eq) {
    std::vector<int64_t> equationIndices;
    getEquationBeginIndices(eq, equationIndices);

    do {
      uint64_t equationArrayOffset = equationOffsets[eq];

      uint64_t equationScalarOffset =
          getEquationFlatIndex(equationIndices, equationRanges[eq]);

      uint64_t scalarEquationIndex = equationArrayOffset + equationScalarOffset;

      // Compute the column indexes that may be non-zeros.
      std::vector<JacobianColumn> jacobianColumns =
          computeJacobianColumns(eq, equationIndices.data());

      jacobianMatrixData[scalarEquationIndex].resize(jacobianColumns.size());

      if (marco::runtime::simulation::getOptions().debug) {
        std::cerr << "  - Equation " << eq << std::endl;
        std::cerr << "    Equation indices: ";
        printIndices(equationIndices);
        std::cerr << std::endl;

        std::cerr << "    Scalar equation index: " << scalarEquationIndex
                  << std::endl;

        std::cerr << "    Number of possibly non-zero columns: "
                  << jacobianColumns.size() << std::endl;
      }
    } while (advanceEquationIndices(equationIndices, equationRanges[eq]));
  }

  // Compute the total amount of non-zero values in the Jacobian Matrix.
  computeNNZ();

  // Compute the equation chunks for each thread.
  computeThreadChunks();

  // Create and initialize the memory for KINSOL.
#if SUNDIALS_VERSION_MAJOR >= 6
  kinsolMemory = KINCreate(ctx);
#else
  kinsolMemory = KINCreate();
#endif

  if (!checkAllocation(kinsolMemory, "KINCreate")) {
    return false;
  }

  if (!kinsolInit()) {
    return false;
  }

  // Create sparse SUNMatrix for use in linear solver.
#if SUNDIALS_VERSION_MAJOR >= 6
  sparseMatrix = SUNSparseMatrix(
      static_cast<sunindextype>(scalarEquationsNumber),
      static_cast<sunindextype>(scalarEquationsNumber),
      static_cast<sunindextype>(nonZeroValuesNumber), CSR_MAT, ctx);
#else
  sparseMatrix =
      SUNSparseMatrix(static_cast<sunindextype>(scalarEquationsNumber),
                      static_cast<sunindextype>(scalarEquationsNumber),
                      static_cast<sunindextype>(nonZeroValuesNumber), CSR_MAT);
#endif

  if (!checkAllocation(static_cast<void *>(sparseMatrix), "SUNSparseMatrix")) {
    return false;
  }

  // Create and attach a KLU SUNLinearSolver object.
#if SUNDIALS_VERSION_MAJOR >= 6
  linearSolver = SUNLinSol_KLU(variablesVector, sparseMatrix, ctx);
#else
  linearSolver = SUNLinSol_KLU(variablesVector, sparseMatrix);
#endif

  if (!checkAllocation(static_cast<void *>(linearSolver), "SUNLinSol_KLU")) {
    return false;
  }

  if (!kinsolSetLinearSolver()) {
    return false;
  }

  if (!kinsolSetUserData() || !kinsolSetJacobianFunction()) {
    return false;
  }

  initialized = true;

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Initialization completed" << std::endl;
  }

  return true;
}

bool KINSOLInstance::solve() {
  if (!initialized) {
    if (!initialize()) {
      return false;
    }
  }

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Computing solution" << std::endl;
  }

  if (getNumOfScalarEquations() == 0) {
    // KINSOL has nothing to solve.
    return true;
  }

  // Update the values of the variables living inside KINSOL.
  copyVariablesFromMARCO(variablesVector);

  auto solveRetVal = KINSol(kinsolMemory, variablesVector, KIN_LINESEARCH,
                            variableScaleVector, residualScaleVector);

  if (solveRetVal != KIN_SUCCESS) {
    // TODO handle errors
    return false;
  }

  copyVariablesIntoMARCO(variablesVector);

  return true;
}

int KINSOLInstance::residualFunction(N_Vector variables, N_Vector residuals,
                                     void *userData) {
  KINSOL_PROFILER_RESIDUALS_CALL_COUNTER_INCREMENT

  realtype *rval = N_VGetArrayPointer(residuals);
  auto *instance = static_cast<KINSOLInstance *>(userData);

  // Copy the values of the variables and derivatives provided by KINSOL into
  // the variables owned by MARCO, so that the residual functions operate on
  // the current iteration values.
  instance->copyVariablesIntoMARCO(variables);

  // For every vectorized equation, set the residual values of the variables
  // it writes into.
  KINSOL_PROFILER_RESIDUALS_START

  instance->equationsParallelIteration(
      [&](Equation eq, const std::vector<int64_t> &equationIndices,
          const JacobianSeedsMap &jacobianSeedsMap) {
        assert(equationIndices.size() == instance->getEquationRank(eq));

        uint64_t equationArrayOffset = instance->equationOffsets[eq];

        uint64_t equationScalarOffset =
            getEquationFlatIndex(equationIndices, instance->equationRanges[eq]);

        uint64_t offset = equationArrayOffset + equationScalarOffset;

        auto residualFn = instance->residualFunctions[eq];
        auto *eqIndicesPtr = equationIndices.data();

        auto residualFunctionResult = residualFn(eqIndicesPtr);
        *(rval + offset) = residualFunctionResult;
      });

  KINSOL_PROFILER_RESIDUALS_STOP

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Residuals function called" << std::endl;
    std::cerr << "Variables:" << std::endl;
    instance->printVariablesVector(variables);
    std::cerr << "Residuals vector:" << std::endl;
    instance->printResidualsVector(residuals);
  }

  return KIN_SUCCESS;
}

int KINSOLInstance::jacobianMatrix(N_Vector variables, N_Vector residuals,
                                   SUNMatrix jacobianMatrix, void *userData,
                                   N_Vector tempv1, N_Vector tempv2) {
  KINSOL_PROFILER_PARTIAL_DERIVATIVES_CALL_COUNTER_INCREMENT

  realtype *jacobian = SUNSparseMatrix_Data(jacobianMatrix);
  auto *instance = static_cast<KINSOLInstance *>(userData);

  // Copy the values of the variables and derivatives provided by KINSOL into
  // the variables owned by MARCO, so that the jacobian functions operate on
  // the current iteration values.
  instance->copyVariablesIntoMARCO(variables);

  KINSOL_PROFILER_PARTIAL_DERIVATIVES_START

  instance->equationsParallelIteration(
      [&](Equation eq, const std::vector<int64_t> &equationIndices,
          const JacobianSeedsMap &jacobianSeedsMap) {
        uint64_t equationArrayOffset = instance->equationOffsets[eq];

        uint64_t equationScalarOffset =
            getEquationFlatIndex(equationIndices, instance->equationRanges[eq]);

        uint64_t scalarEquationIndex =
            equationArrayOffset + equationScalarOffset;

        assert(scalarEquationIndex < instance->getNumOfScalarEquations());

        // Compute the column indexes that may be non-zeros.
        std::vector<JacobianColumn> jacobianColumns =
            instance->computeJacobianColumns(eq, equationIndices.data());

        // For every scalar variable with respect to which the equation must be
        // partially differentiated.
        for (size_t i = 0, e = jacobianColumns.size(); i < e; ++i) {
          const JacobianColumn &column = jacobianColumns[i];
          Variable variable = column.first;
          const auto &variableIndices = column.second;

          uint64_t variableArrayOffset = instance->variableOffsets[variable];

          uint64_t variableScalarOffset = getVariableFlatIndex(
              instance->variablesDimensions[variable], column.second);

          auto jacobianFunction =
              instance->jacobianFunctions[eq][variable].first;

          assert(jacobianFunction != nullptr);
          auto seedsMapIt = jacobianSeedsMap.find(jacobianFunction);
          const uint64_t *seedsPtr = nullptr;

          if (seedsMapIt != jacobianSeedsMap.end()) {
            seedsPtr = seedsMapIt->second.data();
          }

          auto jacobianFunctionResult =
              jacobianFunction(equationIndices.data(), variableIndices.data(),
                               instance->memoryPoolId, seedsPtr);

          instance->jacobianMatrixData[scalarEquationIndex][i].second =
              jacobianFunctionResult;

          auto index = static_cast<sunindextype>(variableArrayOffset +
                                                 variableScalarOffset);

          instance->jacobianMatrixData[scalarEquationIndex][i].first = index;
        }
      });

  // Move the partial derivatives into the SUNDIALS sparse matrix.
  sunindextype *rowPtrs = SUNSparseMatrix_IndexPointers(jacobianMatrix);
  sunindextype *columnIndices = SUNSparseMatrix_IndexValues(jacobianMatrix);

  sunindextype offset = 0;
  *rowPtrs++ = offset;

  for (const auto &row : instance->jacobianMatrixData) {
    offset += static_cast<sunindextype>(row.size());
    *rowPtrs++ = offset;

    for (const auto &column : row) {
      *columnIndices++ = column.first;
      *jacobian++ = column.second;
    }
  }

  assert(rowPtrs == SUNSparseMatrix_IndexPointers(jacobianMatrix) +
                        instance->getNumOfScalarEquations() + 1);
  assert(columnIndices == SUNSparseMatrix_IndexValues(jacobianMatrix) +
                              instance->nonZeroValuesNumber);
  assert(jacobian ==
         SUNSparseMatrix_Data(jacobianMatrix) + instance->nonZeroValuesNumber);

  KINSOL_PROFILER_PARTIAL_DERIVATIVES_STOP

  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Jacobian matrix function called" << std::endl;
    std::cerr << "Variables:" << std::endl;
    instance->printVariablesVector(variables);
    std::cerr << "Residuals vector:" << std::endl;
    instance->printResidualsVector(residuals);
    std::cerr << "Jacobian matrix:" << std::endl;
    instance->printJacobianMatrix(jacobianMatrix);
  }

  return KIN_SUCCESS;
}

uint64_t KINSOLInstance::getNumOfArrayVariables() const {
  return variablesDimensions.size();
}

uint64_t KINSOLInstance::getNumOfScalarVariables() const {
  return scalarVariablesNumber;
}

uint64_t KINSOLInstance::getVariableFlatSize(Variable variable) const {
  uint64_t result = 1;

  for (uint64_t dimension : variablesDimensions[variable]) {
    result *= dimension;
  }

  return result;
}

uint64_t KINSOLInstance::getNumOfVectorizedEquations() const {
  return equationRanges.size();
}

uint64_t KINSOLInstance::getNumOfScalarEquations() const {
  return scalarEquationsNumber;
}

uint64_t KINSOLInstance::getEquationRank(Equation equation) const {
  return equationRanges[equation].size();
}

uint64_t KINSOLInstance::getEquationFlatSize(Equation equation) const {
  assert(equation < getNumOfVectorizedEquations());
  uint64_t result = 1;

  for (const Range &range : equationRanges[equation]) {
    result *= range.end - range.begin;
  }

  return result;
}

uint64_t KINSOLInstance::getVariableRank(Variable variable) const {
  return variablesDimensions[variable].rank();
}

void KINSOLInstance::iterateAccessedArrayVariables(
    Equation equation, std::function<void(Variable)> callback) const {
  if (precomputedAccesses) {
    for (const auto &access : variableAccesses[equation]) {
      callback(access.first);
    }
  } else {
    uint64_t numOfArrayVariables = getNumOfArrayVariables();

    for (Variable variable = 0; variable < numOfArrayVariables; ++variable) {
      callback(variable);
    }
  }
}

/// Determine which of the columns of the current Jacobian row has to be
/// populated, and with respect to which variable the partial derivative has
/// to be performed. The row is determined by the indices of the equation.
std::vector<JacobianColumn>
KINSOLInstance::computeJacobianColumns(Equation eq,
                                       const int64_t *equationIndices) const {
  std::set<JacobianColumn> uniqueColumns;

  if (precomputedAccesses) {
    for (const auto &access : variableAccesses[eq]) {
      Variable variable = access.first;
      AccessFunction accessFunction = access.second;

      uint64_t variableRank = getVariableRank(variable);

      std::vector<uint64_t> variableIndices;
      variableIndices.resize(variableRank, 0);
      accessFunction(equationIndices, variableIndices.data());

      assert([&]() -> bool {
        for (uint64_t i = 0; i < variableRank; ++i) {
          if (variableIndices[i] >= variablesDimensions[variable][i]) {
            return false;
          }
        }

        return true;
      }() && "Access out of bounds");

      uniqueColumns.insert({variable, variableIndices});
    }
  } else {
    for (size_t variableIndex = 0, e = getNumOfArrayVariables();
         variableIndex < e; ++variableIndex) {
      const auto &dimensions = variablesDimensions[variableIndex];

      for (auto indices = dimensions.indicesBegin(),
                end = dimensions.indicesEnd();
           indices != end; ++indices) {
        JacobianColumn column(variableIndex, {});

        for (size_t dim = 0; dim < dimensions.rank(); ++dim) {
          column.second.push_back((*indices)[dim]);
        }

        uniqueColumns.insert(std::move(column));
      }
    }
  }

  std::vector<JacobianColumn> orderedColumns;

  for (const JacobianColumn &column : uniqueColumns) {
    orderedColumns.push_back(column);
  }

  std::sort(orderedColumns.begin(), orderedColumns.end(),
            [](const JacobianColumn &first, const JacobianColumn &second) {
              if (first.first != second.first) {
                return first.first < second.first;
              }

              assert(first.second.size() == second.second.size());

              for (size_t i = 0, e = first.second.size(); i < e; ++i) {
                if (first.second[i] < second.second[i]) {
                  return true;
                }
              }

              return false;
            });

  return orderedColumns;
}

/// Compute the number of non-zero values in the Jacobian Matrix. Also
/// compute the column indexes of all non-zero values in the Jacobian Matrix.
/// This allows to avoid the recomputation of such indexes during the
/// Jacobian evaluation.
void KINSOLInstance::computeNNZ() {
  nonZeroValuesNumber = 0;
  std::vector<int64_t> equationIndices;

  for (size_t eq = 0; eq < getNumOfVectorizedEquations(); ++eq) {
    // Initialize the multidimensional interval of the vector equation.
    uint64_t equationRank = equationRanges[eq].size();
    equationIndices.resize(equationRank);

    for (size_t i = 0; i < equationRank; ++i) {
      const auto &iterationRange = equationRanges[eq][i];
      int64_t beginIndex = iterationRange.begin;
      equationIndices[i] = beginIndex;
    }

    // For every scalar equation in the vector equation.
    do {
      // Compute the column indexes that may be non-zeros
      nonZeroValuesNumber +=
          computeJacobianColumns(eq, equationIndices.data()).size();

    } while (advanceEquationIndices(equationIndices, equationRanges[eq]));
  }
}

void KINSOLInstance::computeThreadChunks() {
  unsigned int numOfThreads = threadPool.getNumOfThreads();

  int64_t chunksFactor = getOptions().equationsChunksFactor;
  int64_t numOfChunks = numOfThreads * chunksFactor;

  uint64_t numOfVectorizedEquations = getNumOfVectorizedEquations();
  uint64_t numOfScalarEquations = getNumOfScalarEquations();

  size_t chunkSize = (numOfScalarEquations + numOfChunks - 1) / numOfChunks;

  // The number of vectorized equations whose indices have been completely
  // assigned.
  uint64_t processedEquations = 0;

  while (processedEquations < numOfVectorizedEquations) {
    Equation equation = processedEquations;
    uint64_t equationFlatSize = getEquationFlatSize(equation);
    uint64_t equationFlatIndex = 0;

    // Divide the ranges into chunks.
    while (equationFlatIndex < equationFlatSize) {
      uint64_t beginFlatIndex = equationFlatIndex;

      uint64_t endFlatIndex = std::min(
          beginFlatIndex + static_cast<uint64_t>(chunkSize), equationFlatSize);

      std::vector<int64_t> beginIndices;
      std::vector<int64_t> endIndices;

      getEquationIndicesFromFlatIndex(beginFlatIndex, beginIndices,
                                      equationRanges[equation]);

      if (endFlatIndex == equationFlatSize) {
        getEquationEndIndices(equation, endIndices);
      } else {
        getEquationIndicesFromFlatIndex(endFlatIndex, endIndices,
                                        equationRanges[equation]);
      }

      JacobianSeedsMap jacobianSeedsMap;

      iterateAccessedArrayVariables(equation, [&](Variable variable) {
        auto jacobianFunction = jacobianFunctions[equation][variable].first;
        const auto &seedSizes = jacobianFunctions[equation][variable].second;

        for (const auto &seedSize : seedSizes) {
          MemoryPool &memoryPool =
              MemoryPoolManager::getInstance().get(memoryPoolId);
          uint64_t seedId = memoryPool.create(seedSize);
          jacobianSeedsMap[jacobianFunction].push_back(seedId);
        }
      });

      threadEquationsChunks.emplace_back(equation, std::move(beginIndices),
                                         std::move(endIndices),
                                         std::move(jacobianSeedsMap));

      // Move to the next chunk.
      equationFlatIndex = endFlatIndex;
    }

    // Move to the next vectorized equation.
    ++processedEquations;
  }
}

void KINSOLInstance::copyVariablesFromMARCO(N_Vector variables) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Copying variables from MARCO" << std::endl;
  }

  KINSOL_PROFILER_COPY_VARS_FROM_MARCO_START

  realtype *varsPtr = N_VGetArrayPointer(variables);
  uint64_t numOfArrayVariables = getNumOfArrayVariables();

  for (Variable var = 0; var < numOfArrayVariables; ++var) {
    uint64_t variableArrayOffset = variableOffsets[var];
    const auto &dimensions = variablesDimensions[var];

    std::vector<uint64_t> varIndices;
    getVariableBeginIndices(var, varIndices);

    do {
      uint64_t variableScalarOffset =
          getVariableFlatIndex(dimensions, varIndices.data());

      uint64_t offset = variableArrayOffset + variableScalarOffset;

      // Get the variable.
      auto getterFn = variableGetters[var];
      auto value = static_cast<realtype>(getterFn(varIndices.data()));
      varsPtr[offset] = value;

      if (marco::runtime::simulation::getOptions().debug) {
        std::cerr << "Got var " << var << " ";
        printIndices(varIndices);
        std::cerr << " with value " << std::fixed << std::setprecision(9)
                  << value << std::endl;
      }
    } while (advanceVariableIndices(varIndices, variablesDimensions[var]));
  }

  KINSOL_PROFILER_COPY_VARS_FROM_MARCO_STOP
}

void KINSOLInstance::copyVariablesIntoMARCO(N_Vector variables) {
  if (marco::runtime::simulation::getOptions().debug) {
    std::cerr << "[KINSOL] Copying variables into MARCO" << std::endl;
  }

  KINSOL_PROFILER_COPY_VARS_INTO_MARCO_START

  realtype *varsPtr = N_VGetArrayPointer(variables);
  uint64_t numOfArrayVariables = getNumOfArrayVariables();

  for (Variable var = 0; var < numOfArrayVariables; ++var) {
    uint64_t variableArrayOffset = variableOffsets[var];
    const auto &dimensions = variablesDimensions[var];

    std::vector<uint64_t> varIndices;
    getVariableBeginIndices(var, varIndices);

    do {
      uint64_t variableScalarOffset =
          getVariableFlatIndex(dimensions, varIndices.data());

      uint64_t offset = variableArrayOffset + variableScalarOffset;

      // Set the variable.
      auto setterFn = variableSetters[var];
      auto value = static_cast<double>(varsPtr[offset]);

      if (marco::runtime::simulation::getOptions().debug) {
        std::cerr << "Setting var " << var << " ";
        printIndices(varIndices);
        std::cerr << " to " << value << std::endl;
      }

      setterFn(value, varIndices.data());

      assert([&]() -> bool {
        auto getterFn = variableGetters[var];
        return getterFn(varIndices.data()) == value;
      }() && "Variable value not set correctly");
    } while (advanceVariableIndices(varIndices, variablesDimensions[var]));
  }

  KINSOL_PROFILER_COPY_VARS_INTO_MARCO_STOP
}

void KINSOLInstance::equationsParallelIteration(
    std::function<void(Equation equation,
                       const std::vector<int64_t> &equationIndices,
                       const JacobianSeedsMap &jacobianSeedsMap)>
        processFn) {
  // Shard the work among multiple threads.
  unsigned int numOfThreads = threadPool.getNumOfThreads();
  std::atomic_size_t chunkIndex = 0;

  for (unsigned int thread = 0; thread < numOfThreads; ++thread) {
    threadPool.async([&]() {
      size_t assignedChunk;

      while ((assignedChunk = chunkIndex++) < threadEquationsChunks.size()) {
        const ThreadEquationsChunk &chunk =
            threadEquationsChunks[assignedChunk];

        Equation equation = std::get<0>(chunk);
        std::vector<int64_t> equationIndices = std::get<1>(chunk);

        do {
          assert([&]() -> bool {
            if (equationIndices.size() != equationRanges[equation].size()) {
              return false;
            }

            for (size_t i = 0, rank = equationIndices.size(); i < rank; ++i) {
              if (equationIndices[i] < equationRanges[equation][i].begin ||
                  equationIndices[i] >= equationRanges[equation][i].end) {
                return false;
              }
            }

            return true;
          }() && "Invalid equation indices");

          processFn(equation, equationIndices, std::get<3>(chunk));
        } while (advanceEquationIndicesUntil(
            equationIndices, equationRanges[equation], std::get<2>(chunk)));
      }
    });
  }

  threadPool.wait();
}

void KINSOLInstance::getVariableBeginIndices(
    Variable variable, std::vector<uint64_t> &indices) const {
  uint64_t variableRank = getVariableRank(variable);
  indices.resize(variableRank);

  for (uint64_t i = 0; i < variableRank; ++i) {
    indices[i] = 0;
  }
}

void KINSOLInstance::getVariableEndIndices(
    Variable variable, std::vector<uint64_t> &indices) const {
  uint64_t variableRank = getVariableRank(variable);
  indices.resize(variableRank);

  for (uint64_t i = 0; i < variableRank; ++i) {
    indices[i] = variablesDimensions[variable][i];
  }
}

void KINSOLInstance::getEquationBeginIndices(
    Equation equation, std::vector<int64_t> &indices) const {
  uint64_t equationRank = getEquationRank(equation);
  indices.resize(equationRank);

  for (uint64_t i = 0; i < equationRank; ++i) {
    indices[i] = equationRanges[equation][i].begin;
  }
}

void KINSOLInstance::getEquationEndIndices(
    Equation equation, std::vector<int64_t> &indices) const {
  uint64_t equationRank = getEquationRank(equation);
  indices.resize(equationRank);

  for (uint64_t i = 0; i < equationRank; ++i) {
    indices[i] = equationRanges[equation][i].end;
  }
}

bool KINSOLInstance::kinsolInit() {
  auto retVal = KINInit(kinsolMemory, residualFunction, variablesVector);

  if (retVal == KIN_MEM_NULL) {
    std::cerr << "KINInit - The kinsol_mem pointer is NULL" << std::endl;
    return false;
  }

  if (retVal == KIN_MEM_FAIL) {
    std::cerr << "KINInit - A memory allocation request has failed"
              << std::endl;
    return false;
  }

  if (retVal == KIN_ILL_INPUT) {
    std::cerr << "KINInit - An input argument to KINInit has an illegal value"
              << std::endl;
    return false;
  }

  return retVal == KIN_SUCCESS;
}

bool KINSOLInstance::kinsolFNTolerance() {
  auto retVal = KINSetFuncNormTol(kinsolMemory, getOptions().fnormtol);

  if (retVal == KIN_MEM_NULL) {
    std::cerr << "KINSVtolerances - The kinsol_mem pointer is NULL"
              << std::endl;
    return false;
  }

  if (retVal == KIN_ILL_INPUT) {
    std::cerr << "KINSVtolerances - The relative error tolerance was negative "
                 "or the absolute tolerance vector had a negative component"
              << std::endl;
    return false;
  }

  return retVal == KIN_SUCCESS;
}

bool KINSOLInstance::kinsolSSTolerance() {
  auto retVal = KINSetScaledStepTol(kinsolMemory, getOptions().scsteptol);

  if (retVal == KIN_MEM_NULL) {
    std::cerr << "KINSVtolerances - The kinsol_mem pointer is NULL"
              << std::endl;
    return false;
  }

  if (retVal == KIN_ILL_INPUT) {
    std::cerr << "KINSVtolerances - The relative error tolerance was negative "
                 "or the absolute tolerance vector had a negative component"
              << std::endl;
    return false;
  }

  return retVal == KIN_SUCCESS;
}

bool KINSOLInstance::kinsolSetLinearSolver() {
  auto retVal = KINSetLinearSolver(kinsolMemory, linearSolver, sparseMatrix);

  if (retVal == KINLS_MEM_NULL) {
    std::cerr << "KINSetLinearSolver - The kinsol_mem pointer is NULL"
              << std::endl;
    return false;
  }

  if (retVal == KINLS_ILL_INPUT) {
    std::cerr << "KINSetLinearSolver - The KINLS interface is not compatible "
                 "with the LS or J input objects or is incompatible with the "
                 "N_Vector object passed to KINInit"
              << std::endl;
    return false;
  }

  if (retVal == KINLS_SUNLS_FAIL) {
    std::cerr << "KINSetLinearSolver - A call to the LS object failed"
              << std::endl;
    return false;
  }

  if (retVal == KINLS_MEM_FAIL) {
    std::cerr << "KINSetLinearSolver - A memory allocation request failed"
              << std::endl;
    return false;
  }

  return retVal == KINLS_SUCCESS;
}

bool KINSOLInstance::kinsolSetUserData() {
  auto retVal = KINSetUserData(kinsolMemory, this);

  if (retVal == KIN_MEM_NULL) {
    std::cerr << "KINSetUserData - The kinsol_mem pointer is NULL" << std::endl;
    return false;
  }

  return retVal == KIN_SUCCESS;
}

bool KINSOLInstance::kinsolSetJacobianFunction() {
  auto retVal = KINSetJacFn(kinsolMemory, jacobianMatrix);

  if (retVal == KIN_MEM_NULL) {
    std::cerr << "KINSetJacFn - The kinsol_mem pointer is NULL" << std::endl;
    return false;
  }

  if (retVal == KINLS_LMEM_NULL) {
    std::cerr << "KINSetJacFn - The KINLS linear solver interface has not been "
                 "initialized"
              << std::endl;
    return false;
  }

  return retVal == KIN_SUCCESS;
}

void KINSOLInstance::printVariablesVector(N_Vector variables) const {
  realtype *data = N_VGetArrayPointer(variables);
  uint64_t numOfArrayVariables = getNumOfArrayVariables();

  for (Variable var = 0; var < numOfArrayVariables; ++var) {
    std::vector<uint64_t> indices;
    getVariableBeginIndices(var, indices);

    do {
      std::cerr << "var " << var << " ";
      printIndices(indices);
      std::cerr << "\t" << std::fixed << std::setprecision(9) << *data
                << std::endl;
      ++data;
    } while (advanceVariableIndices(indices, variablesDimensions[var]));
  }
}

void KINSOLInstance::printResidualsVector(N_Vector residuals) const {
  realtype *data = N_VGetArrayPointer(residuals);
  uint64_t numOfVectorizedEquations = getNumOfVectorizedEquations();

  for (Equation eq = 0; eq < numOfVectorizedEquations; ++eq) {
    std::vector<int64_t> equationIndices;
    getEquationBeginIndices(eq, equationIndices);

    do {
      std::cerr << "eq " << eq << " ";
      printIndices(equationIndices);
      std::cerr << "\t" << std::fixed << std::setprecision(9) << *data << "\n";
      ++data;
    } while (advanceEquationIndices(equationIndices, equationRanges[eq]));
  }
}

void KINSOLInstance::printJacobianMatrix(SUNMatrix jacobianMatrix) const {
  uint64_t numOfArrayVariables = getNumOfArrayVariables();

  // Print the heading row.
  for (Variable var = 0; var < numOfArrayVariables; ++var) {
    std::vector<uint64_t> variableIndices;
    getVariableBeginIndices(var, variableIndices);

    do {
      std::cerr << "\tvar " << var << " ";
      printIndices(variableIndices);
    } while (advanceVariableIndices(variableIndices, variablesDimensions[var]));
  }

  std::cerr << std::endl;

  // Print the matrix.
  uint64_t numOfVectorizedEquations = getNumOfVectorizedEquations();
  uint64_t rowFlatIndex = 0;

  for (Equation eq = 0; eq < numOfVectorizedEquations; ++eq) {
    std::vector<int64_t> equationIndices;
    getEquationBeginIndices(eq, equationIndices);

    do {
      std::cerr << "eq " << eq << " ";
      printIndices(equationIndices);

      uint64_t columnFlatIndex = 0;
      bool previousNegative = false;

      for (Variable var = 0; var < numOfArrayVariables; ++var) {
        std::vector<uint64_t> varIndices;
        getVariableBeginIndices(var, varIndices);

        do {
          auto value = getCellFromSparseMatrix(jacobianMatrix, rowFlatIndex,
                                               columnFlatIndex);

          if (!previousNegative) {
            std::cerr << " ";
          }

          std::cerr << "\t" << std::fixed << std::setprecision(9) << value;
          previousNegative = value < 0;

          columnFlatIndex++;
        } while (advanceVariableIndices(varIndices, variablesDimensions[var]));
      }

      std::cerr << std::endl;
      rowFlatIndex++;
    } while (advanceEquationIndices(equationIndices, equationRanges[eq]));
  }
}
} // namespace marco::runtime::sundials::kinsol

//===---------------------------------------------------------------------===//
// Exported functions
//===---------------------------------------------------------------------===//

//===---------------------------------------------------------------------===//
// kinsolCreate

static void *kinsolCreate_pvoid() {
  auto *instance = new KINSOLInstance();
  return static_cast<void *>(instance);
}

RUNTIME_FUNC_DEF(kinsolCreate, PTR(void))

//===---------------------------------------------------------------------===//
// kinsolSolve

static void kinsolSolve_void(void *instance) {
  [[maybe_unused]] bool result =
      static_cast<KINSOLInstance *>(instance)->solve();
  assert(result && "KINSOL solve failed");
}

RUNTIME_FUNC_DEF(kinsolSolve, void, PTR(void))

//===---------------------------------------------------------------------===//
// kinsolFree

static void kinsolFree_void(void *instance) {
  delete static_cast<KINSOLInstance *>(instance);
}

RUNTIME_FUNC_DEF(kinsolFree, void, PTR(void))

//===---------------------------------------------------------------------===//
// kinsolAddVariable

static uint64_t kinsolAddVariable_i64(void *instance, uint64_t rank,
                                      uint64_t *dimensions, void *getter,
                                      void *setter, void *name) {
  return static_cast<KINSOLInstance *>(instance)->addVariable(
      rank, dimensions, reinterpret_cast<VariableGetter>(getter),
      reinterpret_cast<VariableSetter>(setter),
      static_cast<const char *>(name));
}

RUNTIME_FUNC_DEF(kinsolAddVariable, uint64_t, PTR(void), uint64_t,
                 PTR(uint64_t), PTR(void), PTR(void), PTR(void))

//===---------------------------------------------------------------------===//
// kinsolAddVariableAccess

static void kinsolAddVariableAccess_void(void *instance, uint64_t equationIndex,
                                         uint64_t variableIndex,
                                         void *accessFunction) {
  static_cast<KINSOLInstance *>(instance)->addVariableAccess(
      equationIndex, variableIndex,
      reinterpret_cast<AccessFunction>(accessFunction));
}

RUNTIME_FUNC_DEF(kinsolAddVariableAccess, void, PTR(void), uint64_t, uint64_t,
                 PTR(void))

//===---------------------------------------------------------------------===//
// kinsolAddEquation

static uint64_t kinsolAddEquation_i64(void *instance, int64_t *ranges,
                                      uint64_t rank,
                                      void *stringRepresentation) {
  return static_cast<KINSOLInstance *>(instance)->addEquation(
      ranges, rank, static_cast<const char *>(stringRepresentation));
}

RUNTIME_FUNC_DEF(kinsolAddEquation, uint64_t, PTR(void), PTR(int64_t), uint64_t,
                 PTR(void))

//===---------------------------------------------------------------------===//
// kinsolSetResidual

static void kinsolSetResidual_void(void *instance, uint64_t equationIndex,
                                   void *residualFunction) {
  static_cast<KINSOLInstance *>(instance)->setResidualFunction(
      equationIndex, reinterpret_cast<ResidualFunction>(residualFunction));
}

RUNTIME_FUNC_DEF(kinsolSetResidual, void, PTR(void), uint64_t, PTR(void))

//===---------------------------------------------------------------------===//
// kinsolAddJacobian

static void kinsolAddJacobian_void(void *instance, uint64_t equationIndex,
                                   uint64_t variableIndex,
                                   void *jacobianFunction, uint64_t numOfSeeds,
                                   uint64_t *seedSizes) {
  static_cast<KINSOLInstance *>(instance)->addJacobianFunction(
      equationIndex, variableIndex,
      reinterpret_cast<JacobianFunction>(jacobianFunction), numOfSeeds,
      seedSizes);
}

RUNTIME_FUNC_DEF(kinsolAddJacobian, void, PTR(void), uint64_t, uint64_t,
                 PTR(void), uint64_t, PTR(uint64_t))

#endif // SUNDIALS_ENABLE
