#ifndef MARCO_MATCHING_MODEL_H
#define MARCO_MATCHING_MODEL_H

#include "marco/codegen/dialects/modelica/ModelicaDialect.h"
#include "marco/codegen/passes/model/Equation.h"
#include "marco/codegen/passes/model/Path.h"
#include "marco/codegen/passes/model/Variable.h"

namespace marco::codegen
{
  template<typename EquationType = Equation>
  class Model
  {
    public:
      Model(modelica::ModelOp modelOp)
        : modelOp(modelOp.getOperation())
      {
      }

      modelica::ModelOp getOperation() const
      {
        return mlir::cast<modelica::ModelOp>(modelOp);
      }

      /// Get the variables that are managed by this model.
      Variables getVariables() const
      {
        return variables;
      }

      /// Set the variables the are managed by this model.
      void setVariables(Variables value)
      {
        this->variables = std::move(value);
      }

      /// Get the equations that are managed by this model.
      Equations<EquationType> getEquations() const
      {
        return equations;
      }

      /// Set the equations that are managed by this model.
      void setEquations(Equations<EquationType> value)
      {
        this->equations = std::move(value);
      }

    private:
      mlir::Operation* modelOp;
      Variables variables;
      Equations<EquationType> equations;
  };
}

#endif // MARCO_MATCHING_MODEL_H