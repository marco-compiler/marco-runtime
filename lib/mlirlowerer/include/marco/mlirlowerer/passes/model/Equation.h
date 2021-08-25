#pragma once

#include <llvm/Support/Error.h>
#include <mlir/IR/Builders.h>
#include <mlir/IR/Operation.h>
#include <marco/mlirlowerer/dialects/modelica/ModelicaDialect.h>
#include <marco/utils/IndexSet.hpp>

#include "Expression.h"
#include "Path.h"

namespace marco::codegen::model
{
	class AccessToVar;
	class Expression;
	class ExpressionPath;
	class VectorAccess;

	class Equation
	{
		private:
		class Impl;

		std::shared_ptr<Impl> impl;

		public:
		Equation(mlir::Operation* op,
						 Expression left,
						 Expression right,
						 bool isForward = true,
						 std::optional<EquationPath> path = std::nullopt);

		static Equation build(mlir::Operation* op);
		static Equation build(modelica::EquationOp op);
		static Equation build(modelica::ForEquationOp op);

		bool operator==(const Equation& rhs) const;
		bool operator!=(const Equation& rhs) const;

		bool operator<(const Equation& rhs) const;
		bool operator>(const Equation& rhs) const;
		bool operator<=(const Equation& rhs) const;
		bool operator>=(const Equation& rhs) const;

		[[nodiscard]] modelica::EquationInterface getOp() const;

		[[nodiscard]] Expression lhs() const;
		[[nodiscard]] Expression rhs() const;

		[[nodiscard]] size_t amount() const;

		[[nodiscard]] MultiDimInterval getInductions() const;
		void setInductions(MultiDimInterval inductions);

		[[nodiscard]] size_t dimensions() const;

		[[nodiscard]] bool isForward() const;
		void setForward(bool isForward);

		[[nodiscard]] bool isMatched() const;
		[[nodiscard]] Expression getMatchedExp() const;
		void setMatchedExp(EquationPath path);

		[[nodiscard]] AccessToVar getDeterminedVariable() const;

		[[nodiscard]] ExpressionPath getMatchedExpressionPath() const;

		[[nodiscard]] mlir::LogicalResult normalize();

		mlir::LogicalResult explicitate(mlir::OpBuilder& builder, size_t argumentIndex, bool left);
		mlir::LogicalResult explicitate(const ExpressionPath& path);
		mlir::LogicalResult explicitate();

		[[nodiscard]] Equation clone() const;
		void update();

		[[nodiscard]] Equation composeAccess(const VectorAccess& transformation) const;

		template<typename Path>
		[[nodiscard]] Expression reachExp(Path& path) const
		{
			return path.isOnEquationLeftHand() ? path.reach(lhs()) : path.reach(rhs());
		}

		[[nodiscard]] modelica::EquationSidesOp getTerminator() const;

		private:
		void getEquationsAmount(mlir::ValueRange values, llvm::SmallVectorImpl<long>& amounts) const;
	};
}