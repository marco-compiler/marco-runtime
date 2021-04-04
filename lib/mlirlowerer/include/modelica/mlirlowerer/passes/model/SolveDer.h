#pragma once

#include <modelica/mlirlowerer/passes/model/Model.h>
#include <modelica/mlirlowerer/ModelicaDialect.h>

namespace modelica::codegen::model
{
	class Constant;
	class Expression;
	class Model;
	class Reference;
	class Operation;

	/**
	 * Search for the derivative calls, set the reference variable as a
	 * state one and for each new state allocate a new buffer that will
	 * hold the value that the derivative will have when simulating.
	 */
	class DerSolver
	{
		public:
		DerSolver(SimulationOp simulation, Model& model);

		void solve();

		private:
		void solve(Equation& equation);

		template<typename T>
		void solve(Expression& expression);

		SimulationOp simulation;
		Model& model;
	};

	template<>
	void DerSolver::solve<Expression>(Expression& expression);

	template<>
	void DerSolver::solve<Constant>(Expression& expression);

	template<>
	void DerSolver::solve<Reference>(Expression& expression);

	template<>
	void DerSolver::solve<Operation>(Expression& expression);
}