#include "gtest/gtest.h"

#include "modelica/model/ModMatchers.hpp"
#include "modelica/model/ModParser.hpp"

using namespace modelica;
using namespace llvm;
using namespace std;

TEST(ModMatchersTest, scalarRefVariableCollectorTest)	 // NOLINT
{
	ModParser parser(
			"intVector = INT[3,1] (+ INT[3,1] {3, 3, 3}, INT[3,1] intVector)");

	auto assigment = parser.updateStatement();
	if (!assigment)
		FAIL();

	ReferenceMatcher visitor;

	visit(assigment->getLeft(), visitor);
	visit(assigment->getRight(), visitor);

	EXPECT_EQ(2, visitor.size());
	EXPECT_TRUE(visitor.at(0).isReference());
	EXPECT_TRUE(visitor.at(1).isReference());
	EXPECT_EQ(visitor.at(1).getReference(), "intVector");
	EXPECT_EQ(visitor.at(0).getReference(), "intVector");
}

TEST(ModMatchersTest, arrayRefVariableCollectorTest)	// NOLINT
{
	ModParser parser("negation = INT[1](! INT[1](at INT[2] (at INT[2,2] "
									 "int2Vector, INT[1]{0}), INT[1]{0}))");

	auto assigment = parser.updateStatement();
	if (!assigment)
		FAIL();

	ReferenceMatcher visitor;

	visit(assigment->getLeft(), visitor);
	visit(assigment->getRight(), visitor);

	EXPECT_EQ(2, visitor.size());
	EXPECT_TRUE(visitor.at(0).isReference());
	EXPECT_TRUE(visitor.at(1).isOperation());
	EXPECT_TRUE(visitor.at(1).getLeftHand().isOperation());
	EXPECT_TRUE(visitor.at(1).getLeftHand().getLeftHand().isReference());
}