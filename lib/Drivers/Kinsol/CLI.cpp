#ifdef CLI_ENABLE
#ifdef SUNDIALS_ENABLE

#include "marco/Runtime/Drivers/KINSOL/CLI.h"
#include "marco/Runtime/Solvers/KINSOL/Options.h"

#include <iostream>

namespace marco::runtime::sundials::kinsol {
std::string CommandLineOptions::getTitle() const { return "KINSOL"; }

void CommandLineOptions::printCommandLineOptions(std::ostream &os) const {
  // clang-format off
  os << "  --kinsol-relative-tolerance=<value>     Set the relative tolerance. Defaults to " << getOptions().relativeTolerance << "." << std::endl;
  os << "  --kinsol-absolute-tolerance=<value>     Set the absolute tolerance. Defaults to " << getOptions().absoluteTolerance << "." << std::endl;

  os << "  --kinsol-max-steps=<value>              Set the maximum number of steps to be taken by the solver in its attempt to reach the next output time. Defaults to " << getOptions().maxSteps << "." << std::endl;
  os << "  --kinsol-initial-step-size=<value>      Set the initial step size. Defaults to " << getOptions().initialStepSize << "." << std::endl;
  os << "  --kinsol-min-step-size=<value>          Set the minimum absolute value of the step size. Defaults to " << getOptions().minStepSize << "." << std::endl;
  os << "  --kinsol-max-step-size=<value>          Set the maximum absolute value of the step size. Defaults to " << getOptions().maxStepSize << "." << std::endl;
  os << "  --kinsol-max-err-test-fails=<value>     Set the maximum number of error test failures in attempting one step. Defaults to " << getOptions().maxErrTestFails << "." << std::endl;
  os << "  --kinsol-suppress-alg-vars              Suppress algebraic variables in the local error test." << std::endl;
  os << "  --kinsol-max-nonlin-iters=<value>       Maximum number of nonlinear solver iterations in one solve attempt. Defaults to " << getOptions().maxNonlinIters << "." << std::endl;
  os << "  --kinsol-max-conv-fails=<value>         Maximum number of nonlinear solver convergence failures in one step. Defaults to " << getOptions().maxConvFails << "." << std::endl;
  os << "  --kinsol-nonlin-conv-coef=<value>       Safety factor in the nonlinear convergence test. Defaults to " << getOptions().nonlinConvCoef << "." << std::endl;
  os << "  --kinsol-nonlin-conv-coef-ic=<value>    Positive constant in the Newton iteration convergence test within the initial condition calculation. Defaults to " << getOptions().nonlinConvCoefIC << "." << std::endl;
  os << "  --kinsol-max-steps-ic=<value>           Maximum number of steps allowed for IC. Defaults to " << getOptions().maxStepsIC << "." << std::endl;
  os << "  --kinsol-max-jacs-ic=<value>            Maximum number of the approximate Jacobian or preconditioner evaluations allowed when the Newton iteration appears to be slowly converging. Defaults to " << getOptions().maxNumJacsIC << "." << std::endl;
  os << "  --kinsol-max-iters-ic=<value>           Maximum number of Newton iterations allowed in any one attempt to solve the initial conditions calculation problem. Defaults to " << getOptions().maxNumItersIC << "." << std::endl;
  os << "  --kinsol-line-search-off                Disable the linesearch algorithm" << std::endl;

  os << "  --kinsol-print-jacobian                 Whether to print the Jacobian matrices while debugging." << std::endl;
  // clang-format on
}

void CommandLineOptions::parseCommandLineOptions(
    const argh::parser &options) const {
  // clang-format off
  options("kinsol-relative-tolerance", getOptions().relativeTolerance) >> getOptions().relativeTolerance;
  options("kinsol-absolute-tolerance", getOptions().absoluteTolerance) >> getOptions().absoluteTolerance;
  options("kinsol-max-steps", getOptions().maxSteps) >> getOptions().maxSteps;
  options("kinsol-initial-step-size", getOptions().initialStepSize) >> getOptions().initialStepSize;
  options("kinsol-min-step-size", getOptions().minStepSize) >> getOptions().minStepSize;
  options("kinsol-max-step-size", getOptions().maxStepSize) >> getOptions().maxStepSize;
  options("kinsol-max-err-test-fails", getOptions().maxErrTestFails) >> getOptions().maxErrTestFails;
  getOptions().suppressAlg = options["kinsol-suppress-alg-vars"] ? SUNTRUE : SUNFALSE;
  options("kinsol-max-nonlin-iters", getOptions().maxNonlinIters) >> getOptions().maxNonlinIters;
  options("kinsol-max-conv-fails", getOptions().maxConvFails) >> getOptions().maxConvFails;
  options("kinsol-nonlin-conv-coef", getOptions().nonlinConvCoef) >> getOptions().nonlinConvCoef;
  options("kinsol-nonlin-conv-coef-ic", getOptions().nonlinConvCoefIC) >> getOptions().nonlinConvCoefIC;
  options("kinsol-max-steps-ic", getOptions().maxStepsIC) >> getOptions().maxStepsIC;
  options("kinsol-max-jacs-ic", getOptions().maxNumJacsIC) >> getOptions().maxNumJacsIC;
  options("kinsol-max-iters-ic", getOptions().maxNumItersIC) >> getOptions().maxNumItersIC;
  getOptions().lineSearchOff = options["kinsol-line-search-off"] ? SUNTRUE : SUNFALSE;
  getOptions().printJacobian = options["kinsol-print-jacobian"];
  // clang-format on
}
} // namespace marco::runtime::sundials::kinsol

#endif // SUNDIALS_ENABLE
#endif // CLI_ENABLE
