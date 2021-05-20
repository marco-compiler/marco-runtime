#include <gtest/gtest.h>
#include <modelica/frontend/Parser.h>

using namespace modelica;
using namespace frontend;

TEST(Parser, algorithmWithNoStatements)	 // NOLINT
{
	Parser parser("algorithm");

	auto ast = parser.algorithmSection();
	ASSERT_FALSE(!ast);

	EXPECT_TRUE((*ast)->getBody().empty());
}

TEST(Parser, algorithmStatementsCount)	// NOLINT
{
	Parser parser("algorithm"
								"	x := 1;"
								"	y := 2;"
								"	z := 3;");

	auto ast = parser.algorithmSection();
	ASSERT_FALSE(!ast);

	EXPECT_EQ((*ast)->getBody().size(), 3);
}
