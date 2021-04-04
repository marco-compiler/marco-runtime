#include <mlir/Dialect/StandardOps/IR/Ops.h>
#include <mlir/IR/BlockAndValueMapping.h>
#include <mlir/Transforms/DialectConversion.h>
#include <modelica/mlirlowerer/ModelicaBuilder.h>
#include <modelica/mlirlowerer/ModelicaDialect.h>
#include <modelica/mlirlowerer/passes/SolveModel.h>
#include <modelica/mlirlowerer/passes/matching/Matching.h>
#include <modelica/mlirlowerer/passes/matching/SCCCollapsing.h>
#include <modelica/mlirlowerer/passes/matching/Schedule.h>
#include <modelica/mlirlowerer/passes/matching/SCCDependencyGraph.h>
#include <modelica/mlirlowerer/passes/matching/SVarDependencyGraph.h>
#include <modelica/mlirlowerer/passes/matching/VVarDependencyGraph.h>
#include <modelica/mlirlowerer/passes/model/Equation.h>
#include <modelica/mlirlowerer/passes/model/Expression.h>
#include <modelica/mlirlowerer/passes/model/Model.h>
#include <modelica/mlirlowerer/passes/model/ModelBuilder.h>
#include <modelica/mlirlowerer/passes/model/SolveDer.h>
#include <modelica/mlirlowerer/passes/model/Variable.h>
#include <modelica/utils/Interval.hpp>
#include <variant>

using namespace modelica;
using namespace codegen;
using namespace model;

struct EquationOpScalarizePattern : public mlir::OpRewritePattern<EquationOp>
{
	using mlir::OpRewritePattern<EquationOp>::OpRewritePattern;

	mlir::LogicalResult matchAndRewrite(EquationOp op, mlir::PatternRewriter& rewriter) const override
	{
		mlir::Location location = op->getLoc();
		auto sides = mlir::cast<EquationSidesOp>(op.body().front().getTerminator());
		assert(sides.lhs().size() == 1 && sides.rhs().size() == 1);

		auto lhs = sides.lhs()[0];
		auto rhs = sides.rhs()[0];

		auto lhsPointerType = lhs.getType().cast<PointerType>();
		auto rhsPointerType = rhs.getType().cast<PointerType>();

		llvm::SmallVector<mlir::Value, 3> inductions;

		for (const auto& [left, right] : llvm::zip(lhsPointerType.getShape(), rhsPointerType.getShape()))
		{
			assert(left != -1 || right != -1);

			if (left != -1 && right != -1)
				assert(left == right);

			long size = std::max(left, right);
			mlir::Value induction = rewriter.create<InductionOp>(location, 1, size);
			inductions.push_back(induction);
		}

		auto forEquation = rewriter.create<ForEquationOp>(location, inductions);
		rewriter.mergeBlocks(&op.body().front(), &forEquation.body().front());
		rewriter.setInsertionPoint(sides);

		llvm::SmallVector<mlir::Value, 1> newLhs;
		llvm::SmallVector<mlir::Value, 1> newRhs;

		for (auto [lhs, rhs] : llvm::zip(sides.lhs(), sides.rhs()))
		{
			auto leftSubscription = rewriter.create<SubscriptionOp>(location, lhs, inductions);
			newLhs.push_back(rewriter.create<LoadOp>(location, leftSubscription));

			auto rightSubscription = rewriter.create<SubscriptionOp>(location, rhs, inductions);
			newRhs.push_back(rewriter.create<LoadOp>(location, rightSubscription));
		}

		rewriter.setInsertionPointAfter(sides);
		rewriter.replaceOpWithNewOp<EquationSidesOp>(sides, newLhs, newRhs);
		rewriter.eraseOp(op);
		return mlir::success();
	}
};

struct ForEquationOpScalarizePattern : public mlir::OpRewritePattern<ForEquationOp>
{
	using mlir::OpRewritePattern<ForEquationOp>::OpRewritePattern;

	mlir::LogicalResult matchAndRewrite(ForEquationOp op, mlir::PatternRewriter& rewriter) const override
	{
		mlir::Location location = op->getLoc();
		// TODO
		return mlir::failure();
	}
};

struct SimulationOpPattern : public mlir::OpRewritePattern<SimulationOp>
{
	using mlir::OpRewritePattern<SimulationOp>::OpRewritePattern;

	mlir::LogicalResult matchAndRewrite(SimulationOp op, mlir::PatternRewriter& rewriter) const override
	{
		mlir::Location location = op->getLoc();

		// Initiate the time variable
		mlir::Value startTime = rewriter.create<ConstantOp>(location, op.startTime());
		rewriter.create<StoreOp>(location, startTime, op.time());
		mlir::Value timeStep = rewriter.create<ConstantOp>(location, op.timeStep());

		// Create the loop
		auto loop = rewriter.create<ForOp>(location);

		{
			// Condition
			rewriter.setInsertionPointToStart(&loop.condition().front());
			mlir::Value currentTime = rewriter.create<LoadOp>(location, op.time());
			mlir::Value endTime = rewriter.create<ConstantOp>(location, op.endTime());
			mlir::Value condition = rewriter.create<LtOp>(location, BooleanType::get(op.getContext()), currentTime, endTime);
			rewriter.create<ConditionOp>(location, condition);
		}

		{
			// Body
			assert(op.body().getBlocks().size() == 1);
			rewriter.mergeBlocks(&op.body().front(), &loop.body().front());
		}

		{
			// Step
			rewriter.setInsertionPointToStart(&loop.step().front());
			mlir::Value currentTime = rewriter.create<LoadOp>(location, op.time());
			mlir::Value increasedTime = rewriter.create<AddOp>(location, currentTime.getType(), currentTime, timeStep);
			rewriter.create<StoreOp>(location, increasedTime, op.time());
			rewriter.create<YieldOp>(location);
		}

		rewriter.eraseOp(op);
		return mlir::success();
	}
};

struct EquationOpPattern : public mlir::OpRewritePattern<EquationOp>
{
	using mlir::OpRewritePattern<EquationOp>::OpRewritePattern;

	mlir::LogicalResult matchAndRewrite(EquationOp op, mlir::PatternRewriter& rewriter) const override
	{
		// Create the assignment
		auto sides = mlir::cast<EquationSidesOp>(op.body().front().getTerminator());
		rewriter.setInsertionPoint(sides);

		for (auto [lhs, rhs] : llvm::zip(sides.lhs(), sides.rhs()))
		{
			if (auto loadOp = mlir::dyn_cast<LoadOp>(lhs.getDefiningOp()))
			{
				assert(loadOp.indexes().empty());
				rewriter.create<AssignmentOp>(sides.getLoc(), rhs, loadOp.memory());
			}
			else
			{
				rewriter.create<AssignmentOp>(sides->getLoc(), rhs, lhs);
			}
		}

		rewriter.eraseOp(sides);

		// Inline the equation body
		rewriter.setInsertionPoint(op);
		rewriter.mergeBlockBefore(&op.body().front(), op);

		rewriter.eraseOp(op);
		return mlir::success();
	}
};

struct ForEquationOpPattern : public mlir::OpRewritePattern<ForEquationOp>
{
	using mlir::OpRewritePattern<ForEquationOp>::OpRewritePattern;

	mlir::LogicalResult matchAndRewrite(ForEquationOp op, mlir::PatternRewriter& rewriter) const override
	{

		// Create the assignment
		auto sides = mlir::cast<EquationSidesOp>(op.body().front().getTerminator());
		rewriter.setInsertionPoint(sides);

		for (auto [lhs, rhs] : llvm::zip(sides.lhs(), sides.rhs()))
		{
			if (auto loadOp = mlir::dyn_cast<LoadOp>(lhs.getDefiningOp()))
			{
				assert(loadOp.indexes().empty());
				rewriter.create<AssignmentOp>(sides.getLoc(), rhs, loadOp.memory());
			}
			else
			{
				rewriter.create<AssignmentOp>(sides->getLoc(), rhs, lhs);
			}
		}

		rewriter.eraseOp(sides);

		// Create the loop
		rewriter.setInsertionPoint(op);
		llvm::SmallVector<mlir::Block*, 3> bodies;
		llvm::SmallVector<mlir::Value, 3> inductions;
		mlir::Block* innermostBody = rewriter.getInsertionBlock();

		for (auto induction : op.inductions())
		{
			auto inductionOp = mlir::cast<InductionOp>(induction.getDefiningOp());

			// Init value
			mlir::Value start = rewriter.create<ConstantOp>(induction.getLoc(), rewriter.getIndexAttr(inductionOp.start()));
			start = rewriter.create<SubOp>(induction.getLoc(), start.getType(), start,
																		 rewriter.create<ConstantOp>(induction.getLoc(), rewriter.getIndexAttr(1)));

			auto loop = rewriter.create<ForOp>(induction.getLoc(), start);
			bodies.push_back(&loop.body().front());
			inductions.push_back(loop.body().getArgument(0));

			{
				// Condition
				rewriter.setInsertionPointToStart(&loop.condition().front());
				mlir::Value end = rewriter.create<ConstantOp>(induction.getLoc(), rewriter.getIndexAttr(inductionOp.end()));
				mlir::Value condition = rewriter.create<LtOp>(induction.getLoc(), BooleanType::get(op->getContext()), loop.condition().getArgument(0), end);
				rewriter.create<ConditionOp>(induction.getLoc(), condition, loop.condition().getArgument(0));
			}

			{
				// Step
				rewriter.setInsertionPointToStart(&loop.step().front());
				mlir::Value newInductionValue = rewriter.create<AddOp>(induction.getLoc(),
																															 loop.step().getArgument(0).getType(),
																															 loop.step().getArgument(0),
																															 rewriter.create<ConstantOp>(induction.getLoc(), rewriter.getIndexAttr(1)));
				rewriter.create<YieldOp>(induction.getLoc(), newInductionValue);
			}

			// The next loop will be built inside the current body
			innermostBody = &loop.body().front();
			rewriter.setInsertionPointToStart(innermostBody);
		}

		rewriter.mergeBlocks(&op.body().front(), innermostBody);

		// Add the terminator to each body block
		for (auto [block, induction] : llvm::zip(bodies, inductions))
		{
			rewriter.setInsertionPointToEnd(block);
			rewriter.create<YieldOp>(induction.getLoc(), induction);
		}

		// Remove the original inductions and use the loops ones
		for (auto [original, newValue] : llvm::zip(op.inductions(), inductions))
			rewriter.replaceOp(original.getDefiningOp(), newValue);

		rewriter.eraseOp(op);
		return mlir::success();
	}
};

/**
 * Model solver pass.
 * Its objective is to convert a descriptive (and thus not sequential) model
 * into an algorithmic one.
 */
class SolveModelPass: public mlir::PassWrapper<SolveModelPass, mlir::OperationPass<mlir::ModuleOp>>
{
	public:
	void getDependentDialects(mlir::DialectRegistry &registry) const override
	{
		registry.insert<ModelicaDialect>();
	}

	void runOnOperation() override
	{
		auto module = getOperation();

		// Scalarize the vector assignments
		if (failed(scalarizeArrayEquations()))
			return signalPassFailure();

		module->walk([&](SimulationOp simulation) {
			// Create the model
			Model model(simulation, {}, {});
			ModelBuilder builder(model);

			simulation.walk([&](EquationOp equation) {
				builder.lower(equation);
			});

			simulation.walk([&](ForEquationOp forEquation) {
				builder.lower(forEquation);
			});

			// Remove the derivative operations and allocate the appropriate buffers
			DerSolver solver(simulation, model);
			solver.solve();

			// Match
			if (failed(match(model, 1000)))
				return signalPassFailure();

			// Solve SCC
			if (failed(solveSCC(model, 1000)))
				return signalPassFailure();

			// Schedule
			if (failed(schedule(model)))
				return signalPassFailure();

			// Explicitate the equations so that the updated variable is the only
			// one on the left side of the equation.
			if (failed(explicitateEquations(model)))
				return signalPassFailure();

			// Print the variables
			if (failed(printVariables(model)))
				return signalPassFailure();

			// Calculate the values that the state variables will have in the next
			// iteration.
			if (failed(updateStates(model)))
				return signalPassFailure();
		});

		// The model has been solved and we can now proceed to inline the
		// equations body and create the main simulation loop.

		mlir::ConversionTarget target(getContext());
		target.addLegalDialect<ModelicaDialect>();
		target.addIllegalOp<SimulationOp, EquationOp, ForEquationOp>();

		mlir::OwningRewritePatternList patterns;
		patterns.insert<SimulationOpPattern, EquationOpPattern, ForEquationOpPattern>(&getContext());

		if (failed(applyPartialConversion(module, target, std::move(patterns))))
			return signalPassFailure();
	}

	/**
	 * If and equation consists in an assignment between two arrays, then
	 * convert it into a for equation, in order to scalarize the assignments.
	 */
	mlir::LogicalResult scalarizeArrayEquations()
	{
		mlir::ConversionTarget target(getContext());
		target.addLegalDialect<ModelicaDialect>();

		target.addDynamicallyLegalOp<EquationOp>([](EquationOp op) {
			auto sides = mlir::cast<EquationSidesOp>(op.body().front().getTerminator());
			auto pairs = llvm::zip(sides.lhs(), sides.rhs());

			return std::all_of(pairs.begin(), pairs.end(), [](const auto& pair) {
				mlir::Type lhs = std::get<0>(pair).getType();
				mlir::Type rhs = std::get<1>(pair).getType();

				return !lhs.isa<PointerType>() && !rhs.isa<PointerType>();
			});
		});

		target.addDynamicallyLegalOp<ForEquationOp>([](ForEquationOp op) {
			auto sides = mlir::cast<EquationSidesOp>(op.body().front().getTerminator());
			auto pairs = llvm::zip(sides.lhs(), sides.rhs());

			return std::all_of(pairs.begin(), pairs.end(), [](const auto& pair) {
				mlir::Type lhs = std::get<0>(pair).getType();
				mlir::Type rhs = std::get<1>(pair).getType();

				return !lhs.isa<PointerType>() && !rhs.isa<PointerType>();
			});
		});

		mlir::OwningRewritePatternList patterns;
		patterns.insert<EquationOpScalarizePattern, ForEquationOpScalarizePattern>(&getContext());

		if (auto res = applyPartialConversion(getOperation(), target, std::move(patterns)); failed(res))
			return res;

		return mlir::success();
	}

	mlir::LogicalResult explicitateEquations(Model& model)
	{
		for (auto& equation : model.getEquations())
			if (auto res = equation->explicitate(); failed(res))
				return res;

		return mlir::success();
	}

	mlir::LogicalResult printVariables(Model& model)
	{
		// TODO
		return mlir::success();
	}

	mlir::LogicalResult updateStates(Model& model)
	{
		mlir::OpBuilder builder(model.getOp());
		mlir::Location location = model.getOp()->getLoc();
		mlir::Value timeStep = builder.create<ConstantOp>(location, model.getOp().timeStep());

		builder.setInsertionPoint(model.getOp().body().back().getTerminator());

		for (auto& variable : model.getVariables())
		{
			if (!variable->isState())
				continue;

			mlir::Value var = variable->getReference();
			mlir::Value der = variable->getDer();

			if (auto pointerType = var.getType().cast<PointerType>(); pointerType.getRank() == 0)
				var = builder.create<LoadOp>(location, variable->getReference());

			if (auto pointerType = der.getType().cast<PointerType>(); pointerType.getRank() == 0)
				der = builder.create<LoadOp>(location, variable->getReference());

			mlir::Value newValue = builder.create<MulOp>(location, der.getType(), der, timeStep);
			newValue = builder.create<AddOp>(location, var.getType(), newValue, var);
			builder.create<AssignmentOp>(location, newValue, variable->getReference());
		}

		return mlir::success();
	}
};

std::unique_ptr<mlir::Pass> modelica::codegen::createSolveModelPass()
{
	return std::make_unique<SolveModelPass>();
}
