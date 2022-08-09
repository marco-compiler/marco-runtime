#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "marco/Dialect/Modelica/ModelicaDialect.h"

using namespace ::mlir::modelica;

TEST(ArrayType, staticShape)
{
  mlir::MLIRContext context;
  context.loadDialect<ModelicaDialect>();

  mlir::Type elementType = IntegerType::get(&context);
  auto arrayType = ArrayType::get({ 3, 5 }, elementType);

  EXPECT_EQ(arrayType.getRank(), 2);
  EXPECT_EQ(arrayType.getNumElements(), 15);
  EXPECT_FALSE(arrayType.isDynamicDim(0));
  EXPECT_FALSE(arrayType.isDynamicDim(1));
  EXPECT_TRUE(arrayType.hasStaticShape());
  EXPECT_EQ(arrayType.getNumDynamicDims(), 0);
  EXPECT_EQ(arrayType.getDimSize(0), 3);
  EXPECT_EQ(arrayType.getDimSize(1), 5);
}

TEST(ArrayType, dynamicShape)
{
  mlir::MLIRContext context;
  context.loadDialect<ModelicaDialect>();

  mlir::Type elementType = IntegerType::get(&context);
  auto arrayType = ArrayType::get({ 3, ArrayType::kDynamicSize }, elementType);

  EXPECT_EQ(arrayType.getRank(), 2);
  EXPECT_FALSE(arrayType.isDynamicDim(0));
  EXPECT_TRUE(arrayType.isDynamicDim(1));
  EXPECT_FALSE(arrayType.hasStaticShape());
  EXPECT_EQ(arrayType.getNumDynamicDims(), 1);
  EXPECT_EQ(arrayType.getDimSize(0), 3);
  EXPECT_EQ(arrayType.getDimSize(1), ArrayType::kDynamicSize);
  EXPECT_EQ(arrayType.getDynamicDimIndex(1), 0);
}
