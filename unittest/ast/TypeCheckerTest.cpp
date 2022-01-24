#include "gtest/gtest.h"
#include "llvm/ADT/StringMap.h"
#include "llvm/Support/Error.h"
#include "marco/ast/AST.h"
#include "marco/ast/Parser.h"
#include "marco/ast/passes/TypeCheckingPass.h"
#include "marco/utils/ErrorTest.h"

using namespace llvm;
using namespace marco;
using namespace marco::ast;
using namespace std;

TEST(TypeCheckTest, sumOfIntShouldProduceInt)	 // NOLINT
{
	auto exp = Expression::operation(
			SourceRange::unknown(),
			Type::unknown(),
			OperationKind::add,
			llvm::ArrayRef({
					Expression::constant(SourceRange::unknown(), makeType<int>(), 0),
					Expression::constant(SourceRange::unknown(), makeType<int>(), 4)
			}));

	TypeChecker checker;

	if (checker.run<Expression>(*exp))
		FAIL();

	EXPECT_EQ(exp->getType(), makeType<int>());
}

TEST(TypeCheckTest, andOfBoolShouldProduceBool)	 // NOLINT
{
	auto exp = Expression::operation(
			SourceRange::unknown(),
			Type::unknown(),
			OperationKind::add,
			llvm::ArrayRef({
					Expression::constant(SourceRange::unknown(), makeType<bool>(), true),
					Expression::constant(SourceRange::unknown(), makeType<bool>(), false)
			}));

	TypeChecker checker;

	if (checker.run<Expression>(*exp))
		FAIL();

	EXPECT_EQ(exp->getType(), makeType<bool>());
}

/*
TEST(TypeCheckerTest, tupleExpressionType)	// NOLINT
{
	SmallVector<Expression, 3> expressions;
	expressions.emplace_back(Expression::reference(SourceRange::unknown(), makeType<int>(), "x"));
	expressions.emplace_back(Expression::reference(SourceRange::unknown(), makeType<float>(), "y"));

	Expression exp = Expression::tuple(
			SourceRange::unknown(),
			Type::unknown(),
			expressions);

	Member x(SourceRange::unknown(), "x", makeType<int>(), TypePrefix::none());
	Member y(SourceRange::unknown(), "y", makeType<float>(), TypePrefix::none());

	TypeChecker typeChecker;
	auto& symbolTable = typeChecker.getSymbolTable();
	TypeChecker::SymbolTableScope scope(symbolTable);

	symbolTable.insert(x.getName(), Symbol(x));
	symbolTable.insert(y.getName(), Symbol(y));

	if (typeChecker.check<Tuple>(exp))
		FAIL();

	Type expected({ makeType<int>(), makeType<float>() });
	ASSERT_EQ(expected, exp.getType());
}
*/