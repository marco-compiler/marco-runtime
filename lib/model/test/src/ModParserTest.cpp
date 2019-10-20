#include "gtest/gtest.h"

#include "llvm/Support/Error.h"
#include "modelica/model/ModParser.hpp"

using namespace modelica;
using namespace std;

TEST(ModParserTest, contIntVectorShouldParse)
{
	auto parser = ModParser("{1, 2, 3}");

	auto vec = parser.intVector();
	if (!vec)
		FAIL();

	auto constVector = *vec;

	EXPECT_EQ(constVector.size(), 3);
	EXPECT_EQ(constVector.get(0), 1);
	EXPECT_EQ(constVector.get(1), 2);
	EXPECT_EQ(constVector.get(2), 3);
}

TEST(ModParserTest, contFloatVectorShouldParse)
{
	auto parser = ModParser("{1.4, 2.1, 3.9}");

	auto vec = parser.floatVector();
	if (!vec)
		FAIL();

	auto constVector = *vec;

	EXPECT_EQ(constVector.size(), 3);
	EXPECT_NEAR(constVector.get(0), 1.4f, 0.1f);
	EXPECT_NEAR(constVector.get(1), 2.1f, 0.1f);
	EXPECT_NEAR(constVector.get(2), 3.9f, 0.1f);
}

TEST(ModParserTest, contBoolVectorShouldParse)
{
	auto parser = ModParser("{1, 2, 0}");

	auto vec = parser.boolVector();
	if (!vec)
		FAIL();

	auto constVector = *vec;

	EXPECT_EQ(constVector.size(), 3);
	EXPECT_EQ(constVector.get(0), true);
	EXPECT_EQ(constVector.get(1), true);
	EXPECT_EQ(constVector.get(2), false);
}

TEST(ModParserTest, constExp)
{
	auto parser = ModParser("INT[1]{4, 1, 9}");

	auto vec = parser.expression();
	if (!vec)
		FAIL();

	auto exp = *vec;
	EXPECT_TRUE(exp.isConstant<int>());

	auto& constant = exp.getConstant<int>();

	EXPECT_EQ(constant.size(), 3);
	EXPECT_EQ(constant.get(0), 4);
	EXPECT_EQ(constant.get(1), 1);
	EXPECT_EQ(constant.get(2), 9);
}

TEST(ModParserTest, simCall)
{
	auto parser = ModParser("call fun INT[1](INT[1]{1}, INT[1]{2}, INT[1]{3})");

	auto vec = parser.call();
	if (!vec)
		FAIL();

	auto call = *vec;
	EXPECT_EQ("fun", call.getName());
	EXPECT_EQ(call.getType(), ModType(BultinModTypes::INT));

	EXPECT_EQ(call.argsSize(), 3);
	EXPECT_TRUE(call.at(0).isConstant<int>());
	EXPECT_TRUE(call.at(1).isConstant<int>());
	EXPECT_TRUE(call.at(2).isConstant<int>());
}

TEST(ModParserTest, simCallExp)
{
	auto parser =
			ModParser("FLOAT[1] call fun INT[1](INT[1]{1}, INT[1]{2}, INT[1]{3})");

	auto vec = parser.expression();
	if (!vec)
		FAIL();

	auto exp = *vec;
	EXPECT_TRUE(exp.isCall());
	auto& call = exp.getCall();
	EXPECT_EQ("fun", call.getName());
	EXPECT_EQ(call.getType(), ModType(BultinModTypes::INT));

	EXPECT_EQ(call.argsSize(), 3);
	EXPECT_TRUE(call.at(0).isConstant<int>());
	EXPECT_TRUE(call.at(1).isConstant<int>());
	EXPECT_TRUE(call.at(2).isConstant<int>());
}

TEST(ModParserTest, simRefExp)
{
	auto parser = ModParser("FLOAT[1] ref");

	auto vec = parser.expression();
	if (!vec)
		FAIL();

	auto exp = *vec;
	EXPECT_TRUE(exp.isReference());
	EXPECT_EQ("ref", exp.getReference());
}

TEST(ModParserTest, simOperation)
{
	auto parser = ModParser("FLOAT[1] (+ INT[1]{1}, INT[1]{2})");

	auto vec = parser.expression();
	if (!vec)
		FAIL();

	auto exp = *vec;
	EXPECT_TRUE(exp.isOperation());
	EXPECT_EQ(ModExpKind::add, exp.getKind());
}

TEST(ModParserTest, statement)
{
	auto parser = ModParser("id = FLOAT[1] (+ INT[1]{1}, INT[1]{2})");

	auto vec = parser.statement();
	if (!vec)
		FAIL();

	auto [name, exp] = *vec;
	EXPECT_EQ("id", name);
	EXPECT_TRUE(exp.isOperation());
	EXPECT_EQ(ModExpKind::add, exp.getKind());
}

TEST(ModParserTest, forUpdateStatement)
{
	auto parser =
			ModParser("for [1,3][1,4]id = FLOAT[1] (+ INT[1]{1}, INT[1]{2})");

	auto vec = parser.updateStatement();
	if (!vec)
		FAIL();

	EXPECT_EQ("id", vec->getVarName().getReference());
	EXPECT_TRUE(vec->getExpression().isOperation());
	EXPECT_EQ(ModExpKind::add, vec->getExpression().getKind());
	EXPECT_EQ(vec->getInductionVar(0).begin(), 1);
	EXPECT_EQ(vec->getInductionVar(0).end(), 3);
	EXPECT_EQ(vec->getInductionVar(1).begin(), 1);
	EXPECT_EQ(vec->getInductionVar(1).end(), 4);
}

TEST(ModParserTest, sectionStatement)
{
	auto parser = ModParser("init id = FLOAT[1] (+ INT[1]{1}, INT[1]{2})");

	auto vec = parser.initSection();
	if (!vec)
		FAIL();

	EXPECT_TRUE(vec->find("id") != vec->end());
}

TEST(ModParserTest, updateSection)
{
	auto parser = ModParser("update id = FLOAT[1] (+ INT[1]{1}, INT[1]{2})");

	auto vec = parser.updateSection();
	if (!vec)
		FAIL();

	EXPECT_TRUE(vec.get()[0].getVarName().getReference() == "id");
}

TEST(ModParserTest, simulation)
{
	auto parser = ModParser("init id = FLOAT[1] (+ INT[1]{1}, INT[1]{2}) update "
													"id = FLOAT[1] (+ INT[1]{1}, INT[1]{2})");

	auto vec = parser.simulation();
	if (!vec)
		FAIL();

	auto [init, update] = move(*vec);

	EXPECT_TRUE(init.find("id") != init.end());
	EXPECT_TRUE(update[0].getVarName().getReference() == "id");
	EXPECT_TRUE(init.find("id")->second == update[0].getExpression());
}
