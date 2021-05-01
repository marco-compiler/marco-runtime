#include <gtest/gtest.h>
#include <mlir/Dialect/StandardOps/IR/Ops.h>
#include <mlir/IR/BuiltinOps.h>
#include <modelica/frontend/AST.h>
#include <modelica/mlirlowerer/CodeGen.h>
#include <modelica/mlirlowerer/Runner.h>
#include <modelica/utils/CRunnerUtils.h>
#include <modelica/utils/SourcePosition.h>

using namespace modelica;
using namespace frontend;
using namespace codegen;
using namespace std;

TEST(BuiltInOps, ndims)	 // NOLINT
{
	/**
	 * function main
	 *   input Integer[3] x;
	 *   output Integer y;
	 *
	 *   algorithm
	 *     y := ndims(x);
	 * end main
	 */

	SourcePosition location = SourcePosition::unknown();

	Member xMember(location, "x", makeType<int>(3), TypePrefix(ParameterQualifier::none, IOQualifier::input));
	Member yMember(location, "y", makeType<int>(), TypePrefix(ParameterQualifier::none, IOQualifier::output));

	Statement assignment = AssignmentStatement(
			location,
			Expression::reference(location, makeType<int>(), "y"),
			Expression::call(location, makeType<int>(),
			    Expression::reference(location, makeType<int>(), "ndims"),
											 Expression::reference(location, makeType<int>(3), "x")));

	ClassContainer cls(Function(location, "main", true,
															{ xMember, yMember },
															Algorithm(location, assignment)));

	mlir::MLIRContext context;

	ModelicaOptions modelicaOptions;
	modelicaOptions.x64 = false;
	MLIRLowerer lowerer(context, modelicaOptions);

	auto module = lowerer.lower(cls);

	ModelicaLoweringOptions loweringOptions;
	loweringOptions.llvmOptions.emitCWrappers = true;
	ASSERT_TRUE(module && !failed(lowerer.convertToLLVMDialect(*module, loweringOptions)));

	array<int, 3> x = { 10, 23, -57 };
	ArrayDescriptor<int, 1> xPtr(x.data(), { 3 });

	int y = 0;

	jit::Runner runner(*module);
	ASSERT_TRUE(mlir::succeeded(runner.run("main", xPtr, jit::Runner::result(y))));
	EXPECT_EQ(y, xPtr.getRank());
}

TEST(BuiltInOps, sizeSpecificArrayDimension)	 // NOLINT
{
	/**
	 * function main
	 *   input Integer[3, 2] x;
	 *   output Integer y;
	 *
	 *   algorithm
	 *     y := size(x, 2);
	 * end main
	 */

	SourcePosition location = SourcePosition::unknown();

	Member xMember(location, "x", makeType<int>(3, 2), TypePrefix(ParameterQualifier::none, IOQualifier::input));
	Member yMember(location, "y", makeType<int>(), TypePrefix(ParameterQualifier::none, IOQualifier::output));

	Statement assignment = AssignmentStatement(
			location,
			Expression::reference(location, makeType<int>(), "y"),
			Expression::call(location, makeType<int>(),
											 Expression::reference(location, makeType<int>(), "size"),
											 Expression::reference(location, makeType<int>(3, 2), "x"),
											 Expression::constant(location, makeType<int>(), 2)));

	ClassContainer cls(Function(location, "main", true,
															{ xMember, yMember },
															Algorithm(location, assignment)));

	mlir::MLIRContext context;

	ModelicaOptions modelicaOptions;
	modelicaOptions.x64 = false;
	MLIRLowerer lowerer(context, modelicaOptions);

	auto module = lowerer.lower(cls);

	ModelicaLoweringOptions loweringOptions;
	loweringOptions.llvmOptions.emitCWrappers = true;
	ASSERT_TRUE(module && !failed(lowerer.convertToLLVMDialect(*module, loweringOptions)));

	array<int, 6> x = { 1, 2, 3, 4, 5, 6 };
	ArrayDescriptor<int, 2> xPtr(x.data(), { 3, 2 });

	int y = 0;

	jit::Runner runner(*module);
	ASSERT_TRUE(mlir::succeeded(runner.run("main", xPtr, jit::Runner::result(y))));
	EXPECT_EQ(y, xPtr.getSize(1));
}

TEST(BuiltInOps, sizeAllArrayDimensions)	 // NOLINT
{
	/**
	 * function main
	 *   input Integer[4, 3] x;
	 *   output Integer[2] y;
	 *
	 *   algorithm
	 *     y := size(x);
	 * end main
	 */

	SourcePosition location = SourcePosition::unknown();

	Member xMember(location, "x", makeType<int>(4, 3), TypePrefix(ParameterQualifier::none, IOQualifier::input));
	Member yMember(location, "y", makeType<int>(2), TypePrefix(ParameterQualifier::none, IOQualifier::output));

	Statement assignment = AssignmentStatement(
			location,
			Expression::reference(location, makeType<int>(2), "y"),
			Expression::call(location, makeType<int>(2),
											 Expression::reference(location, makeType<int>(2), "size"),
											 Expression::reference(location, makeType<int>(4, 3), "x")));

	ClassContainer cls(Function(location, "main", true,
															{ xMember, yMember },
															Algorithm(location, assignment)));

	mlir::MLIRContext context;

	ModelicaOptions modelicaOptions;
	modelicaOptions.x64 = false;
	MLIRLowerer lowerer(context, modelicaOptions);

	auto module = lowerer.lower(cls);

	ModelicaLoweringOptions loweringOptions;
	loweringOptions.llvmOptions.emitCWrappers = true;
	ASSERT_TRUE(module && !failed(lowerer.convertToLLVMDialect(*module, loweringOptions)));

	array<int, 12> x = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
	array<int, 2> y = { 0, 0 };

	ArrayDescriptor<int, 2> xPtr(x.data(), { 4, 3 });
	ArrayDescriptor<int, 1> yPtr(y.data(), { 1 });

	jit::Runner runner(*module);
	ASSERT_TRUE(mlir::succeeded(runner.run("main", xPtr, jit::Runner::result(yPtr))));

	EXPECT_EQ(yPtr[0], xPtr.getSize(0));
	EXPECT_EQ(yPtr[1], xPtr.getSize(1));
}

TEST(BuiltInOps, identityMatrix)	 // NOLINT
{
	/**
	 * function main
	 *   input Integer x;
	 *   output Integer[:,:] y;
	 *
	 *   algorithm
	 *     y := identity(x);
	 * end main
	 */

	SourcePosition location = SourcePosition::unknown();

	Member xMember(location, "x", makeType<int>(), TypePrefix(ParameterQualifier::none, IOQualifier::input));
	Member yMember(location, "y", makeType<int>(-1, -1), TypePrefix(ParameterQualifier::none, IOQualifier::output));

	Statement assignment = AssignmentStatement(
			location,
			Expression::reference(location, makeType<int>(-1, -1), "y"),
			Expression::call(location, makeType<int>(-1, -1),
											 Expression::reference(location, makeType<int>(-1, -1), "identity"),
											 Expression::reference(location, makeType<int>(), "x")));

	ClassContainer cls(Function(location, "main", true,
															{ xMember, yMember },
															Algorithm(location, assignment)));

	mlir::MLIRContext context;

	ModelicaOptions modelicaOptions;
	modelicaOptions.x64 = false;
	MLIRLowerer lowerer(context, modelicaOptions);

	auto module = lowerer.lower(cls);

	ModelicaLoweringOptions loweringOptions;
	loweringOptions.llvmOptions.emitCWrappers = true;
	ASSERT_TRUE(module && !failed(lowerer.convertToLLVMDialect(*module, loweringOptions)));

	int x = 3;
	array<int, 3> y = { 2, 2, 2 };

	ArrayDescriptor<int, 1> yPtr(y.data(), { 3 });

	jit::Runner runner(*module);
	ASSERT_TRUE(mlir::succeeded(runner.run("main", x, jit::Runner::result(yPtr))));

	EXPECT_EQ(yPtr.getRank(), 2);

	for (long i = 0; i < 3; ++i)
		for (long j = 0; j < 3; ++j)
			EXPECT_EQ(yPtr.get(i, j), i == j ? 1 : 0);
}

TEST(BuiltInOps, sumOfIntegerStaticArrayValues)	 // NOLINT
{
	/**
	 * function main
	 *   input Integer[3] x;
	 *   output Integer y;
	 *
	 *   algorithm
	 *     y := sum(x);
	 * end main
	 */

	// TODO

	/*
	SourcePosition location = SourcePosition::unknown();

	Member xMember(location, "x", makeType<int>(3), TypePrefix(ParameterQualifier::none, IOQualifier::input));
	Member yMember(location, "y", makeType<int>(), TypePrefix(ParameterQualifier::none, IOQualifier::output));

	Statement assignment = AssignmentStatement(
			location,
			Expression::reference(location, makeType<int>(), "y"),
			Expression::call(location, makeType<int>(),
			    Expression::reference(location, makeType<int>(), "sum"),
											 Expression::reference(location, makeType<int>(3), "x")));

	ClassContainer cls(Function(location, "main", true,
															{ xMember, yMember },
															Algorithm(location, assignment)));

	mlir::MLIRContext context;

	ModelicaOptions modelicaOptions;
	modelicaOptions.x64 = false;
	MLIRLowerer lowerer(context, modelicaOptions);
	
	auto module = lowerer.lower(cls);
	
	ModelicaLoweringOptions loweringOptions;
	loweringOptions.llvmOptions.emitCWrappers = true;
	ASSERT_TRUE(module && !failed(lowerer.convertToLLVMDialect(*module, loweringOptions)));

	jit::Runner runner(module);

	array<int, 3> x = { 1, 2, 3 };
	ArrayDescriptor<int, 1> xPtr(x.data(), { 3 });

	int y = 0;

	if (failed(runner.run("main", xPtr, jit::Runner::result(y))))
		FAIL();

	EXPECT_EQ(y, x[0] + x[1] + x[2]);
	 */
}
