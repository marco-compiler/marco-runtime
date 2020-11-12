#include <gtest/gtest.h>
#include <llvm/Support/Error.h>
#include <modelica/frontend/Expression.hpp>
#include <modelica/frontend/Parser.hpp>
#include <modelica/frontend/Type.hpp>
#include <modelica/frontend/TypeChecker.hpp>
#include <modelica/utils/ErrorTest.hpp>

using namespace std;
using namespace modelica;
using namespace llvm;

TEST(TypeCheckTest, sumOfIntShouldProduceInt)
{
	Expression exp = Expression::op<OperationKind::add>(
			Type::unknown(),
			Expression(makeType<int>(), 0),
			Expression(makeType<int>(), 4));

	TypeChecker checker;

	if (checker.checkType<Expression>(exp, {}))
		FAIL();

	EXPECT_EQ(exp.getType(), makeType<int>());
}

TEST(TypeCheckTest, andOfBoolShouldProduceBool)
{
	Expression exp = Expression::op<OperationKind::add>(
			Type::unknown(),
			Expression(makeType<bool>(), true),
			Expression(makeType<bool>(), false));

	TypeChecker checker;

	if (checker.checkType<Expression>(exp, {}))
		FAIL();

	EXPECT_EQ(exp.getType(), makeType<bool>());
}
