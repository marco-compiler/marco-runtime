#include <llvm/ADT/SmallVector.h>
#include <llvm/Support/raw_ostream.h>
#include <mlir/IR/BlockAndValueMapping.h>
#include <modelica/mlirlowerer/passes/model/Equation.h>
#include <modelica/mlirlowerer/passes/model/Expression.h>
#include <modelica/mlirlowerer/passes/model/LinSolver.h>
#include <modelica/mlirlowerer/passes/model/Model.h>
#include <modelica/mlirlowerer/passes/model/ReferenceMatcher.h>
#include <modelica/mlirlowerer/passes/model/Variable.h>
#include <modelica/mlirlowerer/passes/model/VectorAccess.h>
#include <modelica/utils/IndexSet.hpp>
#include <queue>

using namespace modelica::codegen;
using namespace model;

namespace modelica::codegen::model
{
	static void replaceUses(mlir::OpBuilder& builder, const Equation& newEq, Equation& original)
	{
		mlir::OpBuilder::InsertionGuard guard(builder);

		auto var = newEq.getDeterminedVariable();
		ReferenceMatcher matcher(original);

		for (auto& access : matcher)
		{
			auto pathToVar = AccessToVar::fromExp(access.getExp());

			if (pathToVar.getVar() != var.getVar())
				continue;

			auto composed = newEq.composeAccess(pathToVar.getAccess());

			// Map the old induction values with the ones in the new equation
			mlir::BlockAndValueMapping mapper;

			for (auto [oldInduction, newInduction] :
					llvm::zip(composed.getOp().inductions(), original.getOp().inductions()))
				mapper.map(oldInduction, newInduction);

			// Copy all the operations from the explicitated equation into the
			// one whose member has to be replaced.
			builder.setInsertionPointToStart(original.getOp().body());
			EquationSidesOp clonedTerminator;

			for (auto& op : composed.getOp().body()->getOperations())
			{
				mlir::Operation* clonedOp = builder.clone(op, mapper);

				if (auto terminator = mlir::dyn_cast<EquationSidesOp>(clonedOp))
					clonedTerminator = terminator;
			}

			// Remove the cloned terminator. In fact, in the generic case we need
			// to preserve the original left-hand and right-hand sides of the
			// equations. If the member to be replaced is the same as a side of
			// the original equations, if will be automatically replaced inside
			// the remaining block terminator.
			clonedTerminator.erase();

			// Replace the uses of the value we want to replace.
			for (auto& use : original.reachExp(access).getOp()->getUses())
			{
				// We need to check if we are inside the equation body block. In fact,
				// if the value to be replaced is an array (and not a scalar or a
				// subscription), we would replace the array instantiation itself,
				// which is outside the simulation block and thus would impact also
				// other equations.
				if (!original.getOp()->isAncestor(use.getOwner()))
					continue;

				if (auto loadOp = mlir::dyn_cast<LoadOp>(use.getOwner()); loadOp.indexes().empty())
				{
					// If the value to be replaced is the declaration of a scalar
					// variable, we instead need to replace the load operations which
					// are executed on that variable.
					// Example:
					//  %0 = modelica.alloca : modelica.ptr<int>
					//  equation:
					//    %1 = modelica.load %0 : int
					//    modelica.equation_sides (%1, ...)
					// needs to become
					//  %0 = modelica.alloca : modelica.ptr<int>
					//  equation:
					//    modelica.equation_sides (%newValue, ...)

					loadOp->replaceAllUsesWith(clonedTerminator.rhs()[0].getDefiningOp());
					loadOp->erase();
				}
				else
				{
					use.set(clonedTerminator.rhs()[0]);
				}
			}

			{
				// Prune the remaining useless operations
				mlir::Block::reverse_iterator it(original.getOp().body()->getTerminator());
				auto end = original.getOp().body()->rend();

				while (it != end)
				{
					if (it->getNumResults() != 0 && it->getUses().empty())
					{
						auto curr = it;
						++it;
						curr->erase();
					}
					else
					{
						++it;
					}
				}
			}

			composed.getOp()->erase();
		}

		original.update();
	}

	mlir::LogicalResult linearySolve(mlir::OpBuilder& builder, llvm::SmallVectorImpl<Equation>& equations)
	{
		for (auto eq = equations.rbegin(); eq != equations.rend(); eq++)
			for (auto eq2 = eq + 1; eq2 != equations.rend(); eq2++)
				replaceUses(builder, *eq, *eq2);

		//for (auto& eq : equations)
		//	eq = eq.groupLeftHand();

		return mlir::success();
	}
}

