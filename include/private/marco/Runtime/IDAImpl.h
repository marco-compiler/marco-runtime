#ifndef MARCO_RUNMTIME_IDAIMPL_H
#define MARCO_RUNMTIME_IDAIMPL_H

#include "ida/ida.h"
#include "nvector/nvector_serial.h"
#include "sundials/sundials_config.h"
#include "sundials/sundials_types.h"
#include "sunlinsol/sunlinsol_klu.h"
#include "sunmatrix/sunmatrix_sparse.h"
#include <set>
#include <vector>

namespace marco::runtime::ida
{
  struct Options
  {
    // Relative tolerance is intended as the difference between the values
    // computed through the n-th and the (n+1)-th order BDF method, divided
    // by the absolute value given by the (n+1)-th order BDF method.
    //
    // It is mandatory to set the parameter higher than the minimum
    // precision of the floating point unit roundoff (10^-15 for doubles).
    //
    // It is also highly suggested setting the parameter lower than 10^-3 in
    // order to avoid inaccurate results. IDA defaults to 10^-6.
    realtype relativeTolerance = 1e-06;

    // Absolute tolerance is intended as the maximum acceptable difference
    // between the values computed through the n-th and  the (n+1)-th order
    // BDF method.
    //
    // Absolute tolerance is used to substitute relative tolerance when the
    // value converges to zero. When this happens, in fact, the relative
    // error would tend to infinity, thus exceeding the set tolerance.
    //
    // It is mandatory to set the parameter higher than the minimum precision
    // of the floating point unit roundoff (10^-15 for doubles).
    //
    // It is also highly suggested setting the parameter lower than 10^-3 in
    // order to avoid inaccurate results. IDA defaults to 10^-6.
    realtype absoluteTolerance = 1e-06;

    // Arbitrary initial guess made in the 20/12/2021 Modelica Call
    realtype maxAlgebraicAbsoluteTolerance = 1e-12;

    // Arbitrary initial guess made in the 20/12/2021 Modelica Call
    realtype timeScalingFactorInit = 1e5;

    // Whether to print the Jacobian matrices while debugging
    bool printJacobian = false;

    // Maximum number of steps to reach the next output time
    long maxSteps = 1e4;

    // Initial step size
    realtype initialStepSize = 0;

    // Minimum absolute value of the step size
    realtype minStepSize = 0;

    // Maximum absolute value of the step size
    realtype maxStepSize = 0;

    // Maximum number of error test failures in attempting one step
    int maxErrTestFails = 10;

    // Whether to suppress algebraic variables in the local error test
    booleantype suppressAlg = SUNFALSE;

    // Maximum number of nonlinear solver iterations in one solve attempt
    int maxNonlinIters = 4;

    // Maximum number of nonlinear solver convergence failures in one step
    int maxConvFails = 10;

    // Safety factor in the nonlinear convergence test
    realtype nonlinConvCoef = 0.33;

    // Positive constant in the Newton iteration convergence test within
    // the initial condition calculation.
    realtype nonlinConvCoefIC = 0.0033;

    // Maximum number of steps allowed for IC
    long maxStepsIC = 5;

    // Maximum number of the approximate Jacobian or preconditioner
    // evaluations allowed when the Newton iteration appears to be slowly
    // converging.
    int maxNumJacsIC = 4;

    // Maximum number of Newton iterations allowed in any one attempt
    // to solve the initial conditions calculation problem.
    int maxNumItersIC = 10;

    // Whether to turn on or off the linesearch algorithm
    booleantype lineSearchOff = SUNFALSE;
  };

  Options& getOptions();
}

namespace marco::runtime::ida
{
  using Access = std::vector<std::pair<sunindextype, sunindextype>>;
  using VarAccessList = std::vector<std::pair<sunindextype, Access>>;

  using DerivativeVariable = std::pair<size_t, std::vector<size_t>>;

  class VariableIndicesIterator;

  /// The list of dimensions of an array variable.
  class VariableDimensions
  {
    private:
    using Container = std::vector<size_t>;

    public:
    using iterator = typename Container::iterator;
    using const_iterator = typename Container::const_iterator;

    VariableDimensions(size_t rank);

    size_t rank() const;

    size_t& operator[](size_t index);
    const size_t& operator[](size_t index) const;

    /// @name Dimensions iterators
    /// {

    const_iterator begin() const;
    const_iterator end() const;

    /// }
    /// @name Indices iterators
    /// {

    VariableIndicesIterator indicesBegin() const;
    VariableIndicesIterator indicesEnd() const;

    /// }

    private:
    /// Check that all the dimensions have been correctly initialized.
    [[maybe_unused]] bool isValid() const;

    private:
    Container dimensions;
  };

  /// This class is used to iterate on all the possible combination of indices of a variable.
  class VariableIndicesIterator
  {
    public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = size_t*;
    using difference_type = std::ptrdiff_t;
    using pointer = size_t**;
    using reference = size_t*&;

    ~VariableIndicesIterator();

    static VariableIndicesIterator begin(const VariableDimensions& dimensions);

    static VariableIndicesIterator end(const VariableDimensions& dimensions);

    bool operator==(const VariableIndicesIterator& it) const;

    bool operator!=(const VariableIndicesIterator& it) const;

    VariableIndicesIterator& operator++();
    VariableIndicesIterator operator++(int);

    size_t* operator*() const;

    private:
    VariableIndicesIterator(const VariableDimensions& dimensions);

    void fetchNext();

    private:
    size_t* indices;
    const VariableDimensions* dimensions;
  };

  using EqDimension = std::vector<std::pair<size_t, size_t>>;

  /// Signature of variable getter functions.
  /// The 1st argument is an opaque pointer to the variable descriptor.
  /// The 2nd argument is a pointer to the indices list.
  /// The result is the scalar value.
  template<typename FloatType>
  using VariableGetterFunction = FloatType(*)(void*, size_t*);

  /// Signature of variable setter functions.
  /// The 1st argument is an opaque pointer to the variable descriptor.
  /// The 2nd argument is the value to be set.
  /// The 3rd argument is a pointer to the indices list.
  template<typename FloatType>
  using VariableSetterFunction = void(*)(void*, FloatType, size_t*);

  /// Signature of residual functions.
  /// The 1st argument is the current time.
  /// The 2nd argument is an opaque pointer to the simulation data.
  /// The 3rd argument is a pointer to the list of equation indices.
  /// The result is the residual value.
  template<typename FloatType>
  using ResidualFunction = FloatType(*)(FloatType, void*, size_t*);

  /// Signature of Jacobian functions.
  /// The 1st argument is the current time.
  /// The 2nd argument is an opaque pointer to the simulation data.
  /// The 3rd argument is a pointer to the list of equation indices.
  /// The 4th argument is a pointer to the list of variable indices.
  /// The 5th argument is the 'alpha' value.
  /// The result is the Jacobian value.
  template<typename FloatType>
  using JacobianFunction = FloatType(*)(FloatType, void*, size_t*, size_t*, FloatType);

  class IDAInstance
  {
    public:
      /// Constant used to indicate that no fixed time step has been set.
      static constexpr realtype kUndefinedTimeStep = -1;

      IDAInstance(int64_t marcoBitWidth, int64_t scalarEquationsNumber);

      ~IDAInstance();

      void setStartTime(double time);
      void setEndTime(double time);
      void setTimeStep(double time);

      /// Add and initialize a new variable given its array.
      int64_t addAlgebraicVariable(void* variable, int64_t* dimensions, int64_t rank, void* getter, void* setter);

      int64_t addStateVariable(void* variable, int64_t* dimensions, int64_t rank, void* getter, void* setter);

      void setDerivative(int64_t stateVariable, void* derivative, void* getter, void* setter);

      /// Add the dimension of an equation to the IDA user data.
      int64_t addEquation(int64_t* ranges, int64_t rank);

      void addVariableAccess(int64_t equationIndex, int64_t variableIndex, int64_t* access, int64_t rank);

      /// Add the function pointer that computes the index-th residual function to the
      /// IDA user data.
      void addResidualFunction(int64_t equationIndex, void* residualFunction);

      /// Add the function pointer that computes the index-th jacobian row to the user
      /// data.
      void addJacobianFunction(int64_t equationIndex, int64_t variableIndex, void* jacobianFunction);

      /// Instantiate and initialize all the classes needed by IDA in order to solve
      /// the given system of equations. It also sets optional simulation parameters
      /// for IDA. It must be called before the first usage of idaStep() and after a
      /// call to idaAllocData(). It may fail in case of malformed model.
      bool initialize();

      /// Invoke IDA to perform one step of the computation. If a time step is given,
      /// the output will show the variables in an equidistant time grid based on the
      /// step time parameter. Otherwise, the output will show the variables at every
      /// step of the computation. Returns true if the computation was successful,
      /// false otherwise.
      bool step();

      /// Returns the time reached by the solver after the last step.
      realtype getCurrentTime() const;

      /// Prints statistics regarding the computation of the system.
      void printStatistics() const;

      /// IDAResFn user-defined residual function, passed to IDA through IDAInit.
      /// It contains how to compute the Residual Function of the system, starting
      /// from the provided UserData struct, iterating through every equation.
      static int residualFunction(
          realtype time,
          N_Vector variables, N_Vector derivatives, N_Vector residuals,
          void* userData);

      /// IDALsJacFn user-defined Jacobian approximation function, passed to IDA
      /// through IDASetJacFn. It contains how to compute the Jacobian Matrix of
      /// the system, starting from the provided UserData struct, iterating through
      /// every equation and variable. The matrix is represented in CSR format.
      static int jacobianMatrix(
          realtype time, realtype alpha,
          N_Vector variables, N_Vector derivatives, N_Vector residuals,
          SUNMatrix jacobianMatrix,
          void* userData,
          N_Vector tempv1, N_Vector tempv2, N_Vector tempv3);

    private:
      std::set<DerivativeVariable> computeIndexSet(size_t eq, size_t* eqIndexes) const;

      void computeNNZ();

      void copyVariablesFromMARCO(N_Vector values);

      void copyDerivativesFromMARCO(N_Vector values);

      void copyVariablesIntoMARCO(N_Vector values);

      void copyDerivativesIntoMARCO(N_Vector values);

      /// Prints the Jacobian incidence matrix of the system.
      void printIncidenceMatrix() const;

    private:
      /// @name Forwarded methods
      /// {

      bool idaInit();
      bool idaSVTolerances();
      bool idaSetLinearSolver();
      bool idaSetUserData();
      bool idaSetMaxNumSteps();
      bool idaSetInitialStepSize();
      bool idaSetMinStepSize();
      bool idaSetMaxStepSize();
      bool idaSetStopTime();
      bool idaSetMaxErrTestFails();
      bool idaSetSuppressAlg();
      bool idaSetId();
      bool idaSetJacobianFunction();
      bool idaSetMaxNonlinIters();
      bool idaSetMaxConvFails();
      bool idaSetNonlinConvCoef();
      bool idaSetNonlinConvCoefIC();
      bool idaSetMaxNumStepsIC();
      bool idaSetMaxNumJacsIC();
      bool idaSetMaxNumItersIC();
      bool idaSetLineSearchOffIC();

      /// }

    private:
      SUNContext ctx;
      bool initialized;

      int64_t marcoBitWidth;

      // Model size
      int64_t scalarEquationsNumber;
      int64_t nonZeroValuesNumber;

      // Equations data
      std::vector<EqDimension> equationDimensions;
      std::vector<void*> residuals;
      std::vector<std::vector<void*>> jacobians;
      std::vector<VarAccessList> variableAccesses;

      // The offset of each array variable inside the flattened variables vector
      std::vector<sunindextype> variableOffsets;

      // The dimensions list of each array variable
      std::vector<VariableDimensions> variablesDimensions;
      std::vector<VariableDimensions> derivativesDimensions;

      // Simulation times
      realtype startTime;
      realtype endTime;
      realtype timeStep;
      realtype currentTime;

      // Variables vectors and values
      N_Vector variablesVector;
      N_Vector derivativesVector;

      // The vector stores whether each scalar variable is an algebraic or a state one.
      // 0 = algebraic
      // 1 = state
      N_Vector idVector;

      // The tolerance for each scalar variable.
      N_Vector tolerancesVector;

      // IDA classes
      void* idaMemory;
      SUNMatrix sparseMatrix;
      SUNLinearSolver linearSolver;

      std::vector<void*> variables;
      std::vector<void*> derivatives;

      std::vector<void*> variablesGetters;
      std::vector<void*> derivativesGetters;

      std::vector<void*> variablesSetters;
      std::vector<void*> derivativesSetters;

      void** simulationData;
  };
}

#endif // MARCO_RUNMTIME_IDAIMPL_H