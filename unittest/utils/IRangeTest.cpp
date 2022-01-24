#include "gtest/gtest.h"
#include "marco/utils/IRange.h"
#include <algorithm>

using namespace marco;
using namespace llvm;

TEST(IRangeTest, irangeCanBeUsedInRangeFor)
{
	int a = 0;
	for (auto num : irange(0, 5))
		EXPECT_EQ(num, a++);
	EXPECT_EQ(a, 5);
}

TEST(IRangeTest, irangeCanBeUnsigned)
{
	unsigned a = 0;
	for (auto num : irange<unsigned>(0, 5))
		EXPECT_EQ(num, a++);

	EXPECT_EQ(a, 5);
}

TEST(IRangeTest, irangeWithSingleArgument)
{
	int a = 0;
	for (auto num : irange(5))
		EXPECT_EQ(num, a++);
	EXPECT_EQ(a, 5);
}
