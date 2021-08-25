#include "llvm/Support/CommandLine.h"
#include "llvm/Support/Error.h"
#include "llvm/Support/FileSystem.h"
#include "llvm/Support/MemoryBuffer.h"
#include "marco/frontend/Parser.h"
#include "marco/frontend/Passes.h"
#include "marco/frontend/SymbolTable.hpp"
#include "marco/frontend/passes/TypeCheckingPass.h"
#include "marco/lowerer/Lowerer.hpp"
#include "marco/matching/Matching.hpp"
#include "marco/matching/SccCollapsing.hpp"
#include "marco/matching/Schedule.hpp"
#include "marco/model/AssignModel.hpp"
#include "marco/model/ModVariable.hpp"
#include "marco/model/ScheduledModel.hpp"
#include "marco/omcToModel/OmcToModelPass.hpp"
#include "marco/passes/CleverDAE.hpp"
#include "marco/passes/ConstantFold.hpp"
#include "marco/passes/ForwardEuler.hpp"
#include "marco/passes/SolveModel.hpp"

using namespace marco;
using namespace llvm;
using namespace std;
using namespace cl;

cl::OptionCategory omcCCat("OmcC options");
cl::opt<string> InputFileName(
		cl::Positional, cl::desc("<input-file>"), cl::init("-"), cl::cat(omcCCat));

opt<bool> externalLinkage(
		"publicSymbols",
		cl::desc("globals symbols are set as external linkage"),
		cl::init(false),
		cl::cat(omcCCat));

opt<int> simulationTime(
		"simTime",
		cl::desc("how many ticks the simulation must perform"),
		cl::init(10),
		cl::cat(omcCCat));

opt<double> timeStep(
		"timeStep",
		cl::desc("how long in seconds a ticks in simulation lasts"),
		cl::init(0.1F),
		cl::cat(omcCCat));

opt<bool> dumpModel(
		"d", cl::desc("dump model"), cl::init(false), cl::cat(omcCCat));

opt<bool> dumpFolded(
		"df", cl::desc("dump constant folded"), cl::init(false), cl::cat(omcCCat));

opt<bool> dumpSolvedModel(
		"db",
		cl::desc("dump before lowering to ir model"),
		cl::init(false),
		cl::cat(omcCCat));

opt<bool> dumpSolvedDerModel(
		"dd",
		cl::desc("dump after having removed the derivatives"),
		cl::init(false),
		cl::cat(omcCCat));

opt<bool> dumpCollapsed(
		"dcol",
		cl::desc("dump after having collapsed scc"),
		cl::init(false),
		cl::cat(omcCCat));

opt<bool> dumpMatched(
		"dm", cl::desc("dump matched model"), cl::init(false), cl::cat(omcCCat));

opt<bool> dumpScheduled(
		"dsched",
		cl::desc("dump scheduled model"),
		cl::init(false),
		cl::cat(omcCCat));

opt<bool> dumpTypeChecked(
		"dc",
		cl::desc("dump type checked model"),
		cl::init(false),
		cl::cat(omcCCat));

opt<bool> dumpLowered(
		"l",
		cl::desc("dump lowered model and exit"),
		cl::init(false),
		cl::cat(omcCCat));

opt<string> outputFile(
		"o", cl::desc("<output-file>"), cl::init("-"), cl::cat(omcCCat));

opt<string> entryPointName(
		"entryPoint",
		cl::desc("entry point function name generated by the compiler"),
		cl::init("main"),
		cl::cat(omcCCat));

enum Solver
{
	ForwardEuler,
	CleverDAE
};

opt<Solver> solverName(
		cl::desc("Solvers:"),
		cl::values(
				clEnumValN(ForwardEuler, "forward-euler", "Forward Euler (default)"),
				clEnumValN(CleverDAE, "clever-dae", "Clever DAE")),
		cl::init(ForwardEuler),
		cl::cat(omcCCat));

Expected<AssignModel> selectSolver(ScheduledModel scheduled)
{
	if (solverName == ForwardEuler)
		return addApproximation(scheduled, timeStep);
	if (solverName == CleverDAE)
		return addBltBlocks(scheduled);
	return createStringError(
			errc::executable_format_error, "Could not find the chosen solver");
}

SmallVector<Assigment, 2> toAssign(SmallVector<ModEquation, 2>&& equs)
{
	SmallVector<Assigment, 2> assign;

	for (ModEquation& eq : equs)
	{
		assert(eq.getLeft().isReference() || eq.getLeft().isReferenceAccess());
		assign.emplace_back(
				eq.getTemplate(), move(eq.getInductions()), eq.isForward());
	}

	return assign;
}

ExitOnError exitOnErr;
int main(int argc, char* argv[])
{
	cl::ParseCommandLineOptions(argc, argv);
	auto errorOrBuffer = MemoryBuffer::getFileOrSTDIN(InputFileName);
	error_code error;
	raw_fd_ostream OS(outputFile, error, sys::fs::OF_None);
	if (error)
	{
		errs() << error.message();
	}

	auto buffer = exitOnErr(errorOrToExpected(move(errorOrBuffer)));
	frontend::Parser parser(buffer->getBufferStart());
	auto ast = exitOnErr(parser.classDefinition());

	frontend::TypeChecker checker;
	exitOnErr(checker.run(ast));
	if (dumpTypeChecked)
	{
		ast->dump(OS);
		return 0;
	}

	frontend::ConstantFolder folder;
	exitOnErr(folder.run(ast));
	if (dumpFolded)
	{
		ast->dump(OS);
		return 0;
	}

	Model model;
	OmcToModelPass pass(model);
	exitOnErr(pass.lower<frontend::Class>(*ast, frontend::SymbolTable()));
	if (dumpModel)
	{
		model.dump(OS);
		return 0;
	}

	auto foldedModel = exitOnErr(constantFold(move(model)));
	exitOnErr(solveDer(foldedModel));
	if (dumpSolvedDerModel)
	{
		foldedModel.dump(OS);
		return 0;
	}
	auto matchedModel = exitOnErr(match(move(foldedModel), 1000));
	if (dumpMatched)
	{
		matchedModel.dump(OS);
		return 0;
	}

	auto collapsed = exitOnErr(solveScc(move(matchedModel), 1000));
	if (dumpCollapsed)
	{
		collapsed.dump(OS);
		return 0;
	}

	auto scheduled = schedule(move(collapsed));
	if (dumpScheduled)
	{
		scheduled.dump(OS);
		return 0;
	}

	auto assModel = exitOnErr(selectSolver(scheduled));
	if (dumpSolvedModel)
	{
		assModel.dump(OS);
		return 0;
	}

	LLVMContext context;
	Lowerer sim(
			context,
			move(assModel.getVars()),
			move(assModel.getUpdates()),
			"Modelica Model",
			entryPointName,
			simulationTime,
			true);

	if (externalLinkage)
		sim.setVarsLinkage(GlobalValue::LinkageTypes::ExternalLinkage);

	if (dumpLowered)
	{
		sim.dump(OS);
		return 0;
	}
	exitOnErr(sim.lower());

	sim.verify();
	sim.dumpBC(OS);

	return 0;
}
