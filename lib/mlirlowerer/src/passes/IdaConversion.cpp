#include <mlir/Conversion/SCFToStandard/SCFToStandard.h>
#include <mlir/Dialect/Math/IR/Math.h>
#include <mlir/Dialect/MemRef/IR/MemRef.h>
#include <mlir/Dialect/SCF/SCF.h>
#include <mlir/Dialect/StandardOps/IR/Ops.h>
#include <mlir/Transforms/DialectConversion.h>
#include <marco/mlirlowerer/dialects/ida/IdaDialect.h>
#include <marco/mlirlowerer/dialects/modelica/ModelicaDialect.h>
#include <marco/mlirlowerer/passes/IdaConversion.h>
#include <marco/mlirlowerer/passes/TypeConverter.h>
#include <numeric>

using namespace marco::codegen;
using namespace ida;

static mlir::FuncOp getOrDeclareFunction(mlir::OpBuilder& builder, mlir::ModuleOp module, llvm::StringRef name, mlir::TypeRange results, mlir::TypeRange args)
{
	if (auto foo = module.lookupSymbol<mlir::FuncOp>(name))
		return foo;

	mlir::PatternRewriter::InsertionGuard insertGuard(builder);
	builder.setInsertionPointToStart(module.getBody());
	auto foo = builder.create<mlir::FuncOp>(module.getLoc(), name, builder.getFunctionType(args, results));
	foo.setPrivate();
	return foo;
}

static mlir::FuncOp getOrDeclareFunction(mlir::OpBuilder& builder, mlir::ModuleOp module, llvm::StringRef name, mlir::TypeRange results, mlir::ValueRange args)
{
	return getOrDeclareFunction(builder, module, name, results, args.getTypes());
}

/**
 * Generic conversion pattern that provides some utility functions.
 *
 * @tparam FromOp type of the operation to be converted
 */
template<typename FromOp>
class IdaOpConversion : public mlir::OpConversionPattern<FromOp>
{
	public:
	IdaOpConversion(mlir::MLIRContext* ctx, TypeConverter& typeConverter)
			: mlir::OpConversionPattern<FromOp>(typeConverter, ctx, 1)
	{
	}

	[[nodiscard]] marco::codegen::TypeConverter& typeConverter() const
	{
		return *static_cast<marco::codegen::TypeConverter *>(this->getTypeConverter());
	}

	[[nodiscard]] mlir::Type convertType(mlir::Type type) const
	{
		return typeConverter().convertType(type);
	}
};

struct ConstantValueOpLowering : public IdaOpConversion<ConstantValueOp>
{
	using IdaOpConversion<ConstantValueOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(ConstantValueOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		llvm::Optional<mlir::Attribute> attribute = convertAttribute(rewriter, op.resultType(), op.value());

		if (!attribute)
			return rewriter.notifyMatchFailure(op, "Unknown attribute type");

		rewriter.replaceOpWithNewOp<mlir::ConstantOp>(op, *attribute);
		return mlir::success();
	}

	private:
	llvm::Optional<mlir::Attribute> convertAttribute(mlir::OpBuilder& builder, mlir::Type resultType, mlir::Attribute attribute) const
	{
		resultType = getTypeConverter()->convertType(resultType);

		if (auto booleanAttribute = attribute.dyn_cast<BooleanAttribute>())
			return builder.getBoolAttr(booleanAttribute.getValue());

		if (auto integerAttribute = attribute.dyn_cast<IntegerAttribute>())
			return builder.getIntegerAttr(resultType, integerAttribute.getValue());

		if (auto realAttribute = attribute.dyn_cast<RealAttribute>())
			return builder.getFloatAttr(resultType, realAttribute.getValue());

		if (auto realAttribute = attribute.dyn_cast<modelica::RealAttribute>())
			return builder.getFloatAttr(resultType, realAttribute.getValue());

		return llvm::None;
	}
};

struct AllocUserDataOpLowering : public IdaOpConversion<AllocUserDataOp>
{
	using IdaOpConversion<AllocUserDataOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AllocUserDataOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		OpaquePointerType result = OpaquePointerType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"allocIdaUserData",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct FreeUserDataOpLowering : public IdaOpConversion<FreeUserDataOp>
{
	using IdaOpConversion<FreeUserDataOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(FreeUserDataOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		BooleanType result = BooleanType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"freeIdaUserData",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct SetInitialValueOpLowering : public IdaOpConversion<SetInitialValueOp>
{
	using IdaOpConversion<SetInitialValueOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(SetInitialValueOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"setInitialValue",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct InitOpLowering : public IdaOpConversion<InitOp>
{
	using IdaOpConversion<InitOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(InitOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		BooleanType result = BooleanType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"idaInit",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct StepOpLowering : public IdaOpConversion<StepOp>
{
	using IdaOpConversion<StepOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(StepOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		BooleanType result = BooleanType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"idaStep",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct AddTimeOpLowering : public IdaOpConversion<AddTimeOp>
{
	using IdaOpConversion<AddTimeOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddTimeOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addTime",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct AddToleranceOpLowering : public IdaOpConversion<AddToleranceOp>
{
	using IdaOpConversion<AddToleranceOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddToleranceOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addTolerance",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct AddRowLengthOpLowering : public IdaOpConversion<AddRowLengthOp>
{
	using IdaOpConversion<AddRowLengthOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddRowLengthOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addRowLength",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct AddDimensionOpLowering : public IdaOpConversion<AddDimensionOp>
{
	using IdaOpConversion<AddDimensionOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddDimensionOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addDimension",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct AddResidualOpLowering : public IdaOpConversion<AddResidualOp>
{
	using IdaOpConversion<AddResidualOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddResidualOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addResidual",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct AddJacobianOpLowering : public IdaOpConversion<AddJacobianOp>
{
	using IdaOpConversion<AddJacobianOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddJacobianOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addJacobian",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct GetTimeOpLowering : public IdaOpConversion<GetTimeOp>
{
	using IdaOpConversion<GetTimeOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(GetTimeOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		modelica::RealType result = modelica::RealType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"getIdaTime",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct GetVariableOpLowering : public IdaOpConversion<GetVariableOp>
{
	using IdaOpConversion<GetVariableOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(GetVariableOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		modelica::RealType result = modelica::RealType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"getIdaVariable",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct AddNewLambdaAccessOpLowering : public IdaOpConversion<AddNewLambdaAccessOp>
{
	using IdaOpConversion<AddNewLambdaAccessOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddNewLambdaAccessOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		IntegerType result = IntegerType::get(op.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addNewLambdaAccess",
				result,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), result, op.args());
		return mlir::success();
	}
};

struct AddLambdaAccessOpLowering : public IdaOpConversion<AddLambdaAccessOp>
{
	using IdaOpConversion<AddLambdaAccessOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddLambdaAccessOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addLambdaAccess",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

struct AddLambdaDimensionOpLowering : public IdaOpConversion<AddLambdaDimensionOp>
{
	using IdaOpConversion<AddLambdaDimensionOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(AddLambdaDimensionOp op, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				"addLambdaDimension",
				llvm::None,
				op.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(op, callee.getName(), llvm::None, op.args());
		return mlir::success();
	}
};

template<typename FromOp, typename FromOpLowering>
struct LambdaLikeLowering : public IdaOpConversion<FromOp>
{
	using IdaOpConversion<FromOp>::IdaOpConversion;

	mlir::LogicalResult matchAndRewrite(FromOp fromOp, llvm::ArrayRef<mlir::Value> operands, mlir::ConversionPatternRewriter& rewriter) const override
	{
		mlir::Operation* op = static_cast<mlir::Operation*>(fromOp);
		IntegerType result = IntegerType::get(fromOp.getContext());

		mlir::FuncOp callee = getOrDeclareFunction(
				rewriter,
				op->getParentOfType<mlir::ModuleOp>(),
				FromOpLowering::operationName,
				result,
				fromOp.args());

		rewriter.replaceOpWithNewOp<mlir::CallOp>(fromOp, callee.getName(), result, fromOp.args());
		return mlir::success();
	}
};

struct LambdaConstantOpLowering : public LambdaLikeLowering<LambdaConstantOp, LambdaConstantOpLowering>
{
	using LambdaLikeLowering<LambdaConstantOp, LambdaConstantOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaConstant";
};

struct LambdaTimeOpLowering : public LambdaLikeLowering<LambdaTimeOp, LambdaTimeOpLowering>
{
	using LambdaLikeLowering<LambdaTimeOp, LambdaTimeOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaTime";
};

struct LambdaScalarVariableOpLowering : public LambdaLikeLowering<LambdaScalarVariableOp, LambdaScalarVariableOpLowering>
{
	using LambdaLikeLowering<LambdaScalarVariableOp, LambdaScalarVariableOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaScalarVariable";
};

struct LambdaScalarDerivativeOpLowering : public LambdaLikeLowering<LambdaScalarDerivativeOp, LambdaScalarDerivativeOpLowering>
{
	using LambdaLikeLowering<LambdaScalarDerivativeOp, LambdaScalarDerivativeOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaScalarDerivative";
};

struct LambdaVectorVariableOpLowering : public LambdaLikeLowering<LambdaVectorVariableOp, LambdaVectorVariableOpLowering>
{
	using LambdaLikeLowering<LambdaVectorVariableOp, LambdaVectorVariableOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaVectorVariable";
};

struct LambdaVectorDerivativeOpLowering : public LambdaLikeLowering<LambdaVectorDerivativeOp, LambdaVectorDerivativeOpLowering>
{
	using LambdaLikeLowering<LambdaVectorDerivativeOp, LambdaVectorDerivativeOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaVectorDerivative";
};

struct LambdaAddOpLowering : public LambdaLikeLowering<LambdaAddOp, LambdaAddOpLowering>
{
	using LambdaLikeLowering<LambdaAddOp, LambdaAddOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaAdd";
};

struct LambdaSubOpLowering : public LambdaLikeLowering<LambdaSubOp, LambdaSubOpLowering>
{
	using LambdaLikeLowering<LambdaSubOp, LambdaSubOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaSub";
};

struct LambdaMulOpLowering : public LambdaLikeLowering<LambdaMulOp, LambdaMulOpLowering>
{
	using LambdaLikeLowering<LambdaMulOp, LambdaMulOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaMul";
};

struct LambdaDivOpLowering : public LambdaLikeLowering<LambdaDivOp, LambdaDivOpLowering>
{
	using LambdaLikeLowering<LambdaDivOp, LambdaDivOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaDiv";
};

struct LambdaPowOpLowering : public LambdaLikeLowering<LambdaPowOp, LambdaPowOpLowering>
{
	using LambdaLikeLowering<LambdaPowOp, LambdaPowOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaPow";
};

struct LambdaNegateOpLowering : public LambdaLikeLowering<LambdaNegateOp, LambdaNegateOpLowering>
{
	using LambdaLikeLowering<LambdaNegateOp, LambdaNegateOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaNegate";
};

struct LambdaAbsOpLowering : public LambdaLikeLowering<LambdaAbsOp, LambdaAbsOpLowering>
{
	using LambdaLikeLowering<LambdaAbsOp, LambdaAbsOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaAbs";
};

struct LambdaSignOpLowering : public LambdaLikeLowering<LambdaSignOp, LambdaSignOpLowering>
{
	using LambdaLikeLowering<LambdaSignOp, LambdaSignOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaSign";
};

struct LambdaSqrtOpLowering : public LambdaLikeLowering<LambdaSqrtOp, LambdaSqrtOpLowering>
{
	using LambdaLikeLowering<LambdaSqrtOp, LambdaSqrtOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaSqrt";
};

struct LambdaExpOpLowering : public LambdaLikeLowering<LambdaExpOp, LambdaExpOpLowering>
{
	using LambdaLikeLowering<LambdaExpOp, LambdaExpOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaExp";
};

struct LambdaLogOpLowering : public LambdaLikeLowering<LambdaLogOp, LambdaLogOpLowering>
{
	using LambdaLikeLowering<LambdaLogOp, LambdaLogOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaLog";
};

struct LambdaLog10OpLowering : public LambdaLikeLowering<LambdaLog10Op, LambdaLog10OpLowering>
{
	using LambdaLikeLowering<LambdaLog10Op, LambdaLog10OpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaLog10";
};

struct LambdaSinOpLowering : public LambdaLikeLowering<LambdaSinOp, LambdaSinOpLowering>
{
	using LambdaLikeLowering<LambdaSinOp, LambdaSinOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaSin";
};

struct LambdaCosOpLowering : public LambdaLikeLowering<LambdaCosOp, LambdaCosOpLowering>
{
	using LambdaLikeLowering<LambdaCosOp, LambdaCosOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaCos";
};

struct LambdaTanOpLowering : public LambdaLikeLowering<LambdaTanOp, LambdaTanOpLowering>
{
	using LambdaLikeLowering<LambdaTanOp, LambdaTanOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaTan";
};

struct LambdaAsinOpLowering : public LambdaLikeLowering<LambdaAsinOp, LambdaAsinOpLowering>
{
	using LambdaLikeLowering<LambdaAsinOp, LambdaAsinOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaAsin";
};

struct LambdaAcosOpLowering : public LambdaLikeLowering<LambdaAcosOp, LambdaAcosOpLowering>
{
	using LambdaLikeLowering<LambdaAcosOp, LambdaAcosOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaAcos";
};

struct LambdaAtanOpLowering : public LambdaLikeLowering<LambdaAtanOp, LambdaAtanOpLowering>
{
	using LambdaLikeLowering<LambdaAtanOp, LambdaAtanOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaAtan";
};

struct LambdaSinhOpLowering : public LambdaLikeLowering<LambdaSinhOp, LambdaSinhOpLowering>
{
	using LambdaLikeLowering<LambdaSinhOp, LambdaSinhOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaSinh";
};

struct LambdaCoshOpLowering : public LambdaLikeLowering<LambdaCoshOp, LambdaCoshOpLowering>
{
	using LambdaLikeLowering<LambdaCoshOp, LambdaCoshOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaCosh";
};

struct LambdaTanhOpLowering : public LambdaLikeLowering<LambdaTanhOp, LambdaTanhOpLowering>
{
	using LambdaLikeLowering<LambdaTanhOp, LambdaTanhOpLowering>::LambdaLikeLowering;
	static constexpr llvm::StringRef operationName = "lambdaTanh";
};

static void populateIdaConversionPatterns(
		mlir::OwningRewritePatternList& patterns,
		mlir::MLIRContext* context,
		marco::codegen::TypeConverter& typeConverter)
{
	// Allocation, initialization, usage and deletion.
	patterns.insert<
			ConstantValueOpLowering,
			AllocUserDataOpLowering,
			FreeUserDataOpLowering,
			SetInitialValueOpLowering,
			InitOpLowering,
			StepOpLowering>(context, typeConverter);

	// Setters.
	patterns.insert<
			AddTimeOpLowering,
			AddToleranceOpLowering,
			AddRowLengthOpLowering,
			AddDimensionOpLowering,
			AddResidualOpLowering,
			AddJacobianOpLowering>(context, typeConverter);
	
	// Getters.
	patterns.insert<
			GetTimeOpLowering,
			GetVariableOpLowering>(context, typeConverter);

	// Lambda helpers.
	patterns.insert<
			AddNewLambdaAccessOpLowering,
			AddLambdaAccessOpLowering,
			AddLambdaDimensionOpLowering>(context, typeConverter);

	// Lambda constructions.
	patterns.insert<
			LambdaConstantOpLowering,
			LambdaTimeOpLowering,
			LambdaScalarVariableOpLowering,
			LambdaScalarDerivativeOpLowering,
			LambdaVectorVariableOpLowering,
			LambdaVectorDerivativeOpLowering>(context, typeConverter);

	patterns.insert<
			LambdaAddOpLowering,
			LambdaSubOpLowering,
			LambdaMulOpLowering,
			LambdaDivOpLowering,
			LambdaPowOpLowering,
			LambdaNegateOpLowering,
			LambdaAbsOpLowering,
			LambdaSignOpLowering,
			LambdaSqrtOpLowering,
			LambdaExpOpLowering,
			LambdaLogOpLowering,
			LambdaLog10OpLowering>(context, typeConverter);

	patterns.insert<
			LambdaSinOpLowering,
			LambdaCosOpLowering,
			LambdaTanOpLowering,
			LambdaAsinOpLowering,
			LambdaAcosOpLowering,
			LambdaAtanOpLowering,
			LambdaSinhOpLowering,
			LambdaCoshOpLowering,
			LambdaTanhOpLowering>(context, typeConverter);
}

class IdaConversionPass : public mlir::PassWrapper<IdaConversionPass, mlir::OperationPass<mlir::ModuleOp>>
{
	public:
	explicit IdaConversionPass(unsigned int bitWidth)
			: bitWidth(bitWidth)
	{
	}

	void getDependentDialects(mlir::DialectRegistry &registry) const override
	{
		registry.insert<IdaDialect>();
	}

	void runOnOperation() override
	{
		if (mlir::failed(convertOperations()))
		{
			mlir::emitError(getOperation().getLoc(), "Error in converting the Ida operations\n");
			return signalPassFailure();
		}
	}

	private:
	mlir::LogicalResult convertOperations()
	{
		auto module = getOperation();
		mlir::ConversionTarget target(getContext());

		// Allocation, initialization, usage and deletion.
		target.addIllegalOp<
				ConstantValueOp,
				AllocUserDataOp,
				FreeUserDataOp,
				SetInitialValueOp,
				InitOp,
				StepOp>();

		// Setters.
		target.addIllegalOp<
				AddTimeOp,
				AddToleranceOp,
				AddRowLengthOp,
				AddDimensionOp,
				AddResidualOp,
				AddJacobianOp>();

		// Getters.
		target.addIllegalOp<GetTimeOp, GetVariableOp>();

		// Lambda helpers.
		target.addIllegalOp<
				AddNewLambdaAccessOp,
				AddLambdaAccessOp,
				AddLambdaDimensionOp>();

		// Lambda constructions.
		target.addIllegalOp<
				LambdaConstantOp,
				LambdaTimeOp,
				LambdaScalarVariableOp,
				LambdaScalarDerivativeOp,
				LambdaVectorVariableOp,
				LambdaVectorDerivativeOp>();

		target.addIllegalOp<
				LambdaAddOp,
				LambdaSubOp,
				LambdaMulOp,
				LambdaDivOp,
				LambdaPowOp,
				LambdaNegateOp,
				LambdaAbsOp,
				LambdaSignOp,
				LambdaSqrtOp,
				LambdaExpOp,
				LambdaLogOp,
				LambdaLog10Op>();

		target.addIllegalOp<
				LambdaSinOp,
				LambdaCosOp,
				LambdaTanOp,
				LambdaAsinOp,
				LambdaAcosOp,
				LambdaAtanOp,
				LambdaSinhOp,
				LambdaCoshOp,
				LambdaTanhOp>();

		target.markUnknownOpDynamicallyLegal([](mlir::Operation* op) { return true; });

		mlir::LowerToLLVMOptions llvmLoweringOptions(&getContext());
		TypeConverter typeConverter(&getContext(), llvmLoweringOptions, bitWidth);

		mlir::OwningRewritePatternList patterns(&getContext());
		populateIdaConversionPatterns(patterns, &getContext(), typeConverter);

		if (auto status = applyPartialConversion(module, target, std::move(patterns)); failed(status))
			return status;

		return mlir::success();
	}

	unsigned int bitWidth;
};

std::unique_ptr<mlir::Pass> marco::codegen::createIdaConversionPass(unsigned int bitWidth)
{
	return std::make_unique<IdaConversionPass>(bitWidth);
}
