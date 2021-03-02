#pragma once

#include <mlir/IR/OpDefinition.h>
#include <modelica/mlirlowerer/ops/OpTrait.h>
#include <modelica/mlirlowerer/Type.h>

#include "ops/MathOps.h"

namespace modelica
{
	/**
	 * Generic operation adaptor.
	 *
	 * The purpose of an adaptor is to allow to access specific values in both
	 * the operation and the conversion pattern without relying on hard-coded
	 * constants in both places.
	 *
	 * @tparam OpType operation class
	 */
	template<typename OpType>
	class OpAdaptor
	{
		public:
		OpAdaptor(mlir::ValueRange values, mlir::DictionaryAttr attrs = nullptr)
				: values(values), attrs(attrs)
		{
		}

		OpAdaptor(OpType& op)
				: values(op->getOperands()), attrs(op->getAttrDictionary())
		{
		}

		protected:
		[[nodiscard]] mlir::ValueRange getValues() const
		{
			return values;
		}

		[[nodiscard]] mlir::DictionaryAttr getAttrs() const
		{
			return attrs;
		}

		private:
		mlir::ValueRange values;
		mlir::DictionaryAttr attrs;
	};

	//===----------------------------------------------------------------------===//
	// Modelica::AssignmentOp
	//===----------------------------------------------------------------------===//

	class AssignmentOp;

	class AssignmentOpAdaptor : public OpAdaptor<AssignmentOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value source();
		mlir::Value destination();
	};

	class AssignmentOp : public mlir::Op<AssignmentOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::ZeroResult, mlir::OpTrait::VariadicOperands> {
		public:
		using Op::Op;
		using Adaptor = AssignmentOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder &odsBuilder, mlir::OperationState &odsState, mlir::Value source, mlir::Value destination);
		void print(mlir::OpAsmPrinter &p);

		mlir::Value source();
		mlir::Value destination();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::AllocaOp
	//===----------------------------------------------------------------------===//

	class AllocaOp;

	class AllocaOpAdaptor : public OpAdaptor<AllocaOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::ValueRange dimensions();
	};

	class AllocaOp : public mlir::Op<AllocaOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::VariadicOperands, mlir::OpTrait::OneResult>
	{
		public:
		using Op::Op;
		using Adaptor = AllocaOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type elementType, PointerType::Shape shape = {}, mlir::ValueRange dimensions = {});
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::AllocOp
	//===----------------------------------------------------------------------===//

	class AllocOp;

	class AllocOpAdaptor : public OpAdaptor<AllocOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::ValueRange dimensions();
	};

	class AllocOp : public mlir::Op<AllocOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::VariadicOperands, mlir::OpTrait::OneResult>
	{
		public:
		using Op::Op;
		using Adaptor = AllocOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type elementType, PointerType::Shape shape = {}, mlir::ValueRange dimensions = {});
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::FreeOp
	//===----------------------------------------------------------------------===//

	class FreeOp;

	class FreeOpAdaptor : public OpAdaptor<FreeOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value memory();
	};

	class FreeOp : public mlir::Op<FreeOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::OneOperand, mlir::OpTrait::ZeroResult>
	{
		public:
		using Op::Op;
		using Adaptor = FreeOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value memory);
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value memory();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::DimOp
	//===----------------------------------------------------------------------===//

	class DimOp;

	class DimOpAdaptor : public OpAdaptor<DimOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value memory();
		mlir::Value dimension();
	};

	class DimOp : public mlir::Op<DimOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult>
	{
		public:
		using Op::Op;
		using Adaptor = DimOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value memory, mlir::Value dimension);
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
		mlir::Value memory();
		mlir::Value dimension();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::SubscriptionOp
	//===----------------------------------------------------------------------===//

	class SubscriptionOp;

	class SubscriptionOpAdaptor : public OpAdaptor<SubscriptionOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value source();
		mlir::ValueRange indexes();
	};

	class SubscriptionOp : public mlir::Op<SubscriptionOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::AtLeastNOperands<2>::Impl, mlir::OpTrait::OneResult>
	{
		public:
		using Op::Op;
		using Adaptor = SubscriptionOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value source, mlir::ValueRange indexes);
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
		mlir::Value source();
		mlir::ValueRange indexes();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::LoadOp
	//===----------------------------------------------------------------------===//

	class LoadOp;

	class LoadOpAdaptor : public OpAdaptor<LoadOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value memory();
		mlir::ValueRange indexes();
	};

	class LoadOp : public mlir::Op<LoadOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::AtLeastNOperands<1>::Impl, mlir::OpTrait::OneResult>
	{
		public:
		using Op::Op;
		using Adaptor = LoadOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value memory, mlir::ValueRange indexes = {});
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
		mlir::Value memory();
		mlir::ValueRange indexes();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::StoreOp
	//===----------------------------------------------------------------------===//

	class StoreOp;

	class StoreOpAdaptor : public OpAdaptor<StoreOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value value();
		mlir::Value memory();
		mlir::ValueRange indexes();
	};

	class StoreOp :public mlir::Op<StoreOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::AtLeastNOperands<2>::Impl, mlir::OpTrait::ZeroResult>
	{
		public:
		using Op::Op;
		using Adaptor = StoreOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value value, mlir::Value memory, mlir::ValueRange indexes = {});
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
		mlir::Value value();
		mlir::Value memory();
		mlir::ValueRange indexes();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::ArrayCopyOp
	//===----------------------------------------------------------------------===//

	class ArrayCopyOp;

	class ArrayCopyOpAdaptor : public OpAdaptor<ArrayCopyOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value source();
	};

	class ArrayCopyOp :public mlir::Op<ArrayCopyOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::OneOperand, mlir::OpTrait::OneResult>
	{
		public:
		using Op::Op;
		using Adaptor = ArrayCopyOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value source, bool heap);
		void print(mlir::OpAsmPrinter& printer);

		PointerType getPointerType();
		mlir::Value source();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::IfOp
	//===----------------------------------------------------------------------===//

	class IfOp;

	class IfOpAdaptor : public OpAdaptor<IfOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value condition();
	};

	class IfOp : public mlir::Op<IfOp, mlir::OpTrait::NRegions<2>::Impl, mlir::OpTrait::ZeroResult, mlir::OpTrait::ZeroSuccessor, mlir::OpTrait::OneOperand> {
		public:
		using Op::Op;
		using Adaptor = IfOpAdaptor;

		static ::llvm::StringRef getOperationName();
		static void build(::mlir::OpBuilder &odsBuilder, ::mlir::OperationState &odsState, mlir::Value cond, bool withElseRegion = false);
		void print(::mlir::OpAsmPrinter &p);

		mlir::Value condition();
		mlir::Region& thenRegion();
		mlir::Region& elseRegion();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::ForOp
	//===----------------------------------------------------------------------===//

	class ForOp;

	class ForOpAdaptor : public OpAdaptor<ForOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value breakCondition();
		mlir::Value returnCondition();
		mlir::ValueRange args();
	};

	class ForOp : public mlir::Op<ForOp, mlir::OpTrait::NRegions<3>::Impl, mlir::OpTrait::AtLeastNOperands<2>::Impl, mlir::OpTrait::ZeroResult>
	{
		public:
		using Op::Op;
		using Adaptor = ForOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value breakCondition, mlir::Value returnCondition);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value breakCondition, mlir::Value returnCondition, mlir::ValueRange args);
		void print(mlir::OpAsmPrinter& printer);

		mlir::Region& condition();
		mlir::Region& step();
		mlir::Region& body();

		mlir::Value breakCondition();
		mlir::Value returnCondition();
		mlir::ValueRange args();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::WhileOp
	//===----------------------------------------------------------------------===//

	class WhileOp;

	class WhileOpAdaptor : public OpAdaptor<WhileOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value breakCondition();
		mlir::Value returnCondition();
	};

	class WhileOp : public mlir::Op<WhileOp, mlir::OpTrait::NRegions<3>::Impl, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::ZeroResult>
	{
		public:
		using Op::Op;
		using Adaptor = WhileOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value breakCondition, mlir::Value returnCondition);
		void print(mlir::OpAsmPrinter& printer);

		mlir::Region& condition();
		mlir::Region& body();

		mlir::Value breakCondition();
		mlir::Value returnCondition();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::ConditionOp
	//===----------------------------------------------------------------------===//

	class ConditionOp;

	class ConditionOpAdaptor : public OpAdaptor<ConditionOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value condition();
		mlir::ValueRange args();
	};

	class ConditionOp : public mlir::Op<ConditionOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::VariadicOperands, mlir::OpTrait::ZeroResult, mlir::OpTrait::ZeroSuccessor, mlir::OpTrait::HasParent<ForOp, WhileOp>::Impl, mlir::OpTrait::IsTerminator> {
		public:
		using Op::Op;
		using Adaptor = ConditionOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder &odsBuilder, ::mlir::OperationState &odsState, mlir::Value condition, mlir::ValueRange args = {});
		void print(::mlir::OpAsmPrinter &p);

		mlir::Value condition();
		mlir::ValueRange args();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::YieldOp
	//===----------------------------------------------------------------------===//

	class YieldOp;

	class YieldOpAdaptor : public OpAdaptor<YieldOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::ValueRange args();
	};

	class YieldOp : public mlir::Op<YieldOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::VariadicOperands, mlir::OpTrait::ZeroResult, mlir::OpTrait::HasParent<IfOp, ForOp, WhileOp>::Impl, mlir::OpTrait::IsTerminator>
	{
		public:
		using Op::Op;
		using Adaptor = YieldOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::ValueRange args = {});
		void print(mlir::OpAsmPrinter& printer);

		mlir::ValueRange args();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::CastOp
	//===----------------------------------------------------------------------===//

	class CastOp;

	class CastOpAdaptor : public OpAdaptor<CastOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value value();
	};

	class CastOp : public mlir::Op<CastOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::OneOperand, mlir::OpTrait::OneResult> {
		public:
		using Op::Op;
		using Adaptor = CastOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& Builder, mlir::OperationState& state, mlir::Value value, mlir::Type resultType);
		void print(mlir::OpAsmPrinter &p);

		mlir::Value value();
		mlir::Type resultType();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::CastCommonOp
	//===----------------------------------------------------------------------===//

	class CastCommonOp;

	class CastCommonOpAdaptor : public OpAdaptor<CastCommonOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::ValueRange operands();
	};

	class CastCommonOp : public mlir::Op<CastCommonOp, mlir::OpTrait::ZeroRegion, mlir::OpTrait::VariadicOperands, mlir::OpTrait::VariadicResults> {
		public:
		using Op::Op;
		using Adaptor = CastCommonOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::ValueRange values);
		void print(mlir::OpAsmPrinter &p);

		mlir::Type resultType();
		mlir::ValueRange operands();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::NegateOp
	//===----------------------------------------------------------------------===//

	class NegateOp;

	class NegateOpAdaptor : public OpAdaptor<NegateOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value operand();
	};

	class NegateOp : public mlir::Op<NegateOp, mlir::OpTrait::OneOperand, mlir::OpTrait::OneResult, mlir::OpTrait::SameOperandsAndResultType, mlir::OpTrait::IsInvolution>
	{
		public:
		using Op::Op;
		using Adaptor = NegateOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder &builder, mlir::OperationState& state, mlir::Value operand);
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value operand();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::EqOp
	//===----------------------------------------------------------------------===//

	class EqOp;

	class EqOpAdaptor : public OpAdaptor<EqOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class EqOp : public mlir::Op<EqOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, OperandsAreSignlessIntegerOrFloatLike, mlir::OpTrait::SameTypeOperands, mlir::OpTrait::ResultsAreBoolLike, mlir::OpTrait::IsCommutative>
	{
		public:
		using Op::Op;
		using Adaptor = EqOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value lhs, mlir::Value rhs);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		mlir::LogicalResult verify();
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::NotEqOp
	//===----------------------------------------------------------------------===//

	class NotEqOp;

	class NotEqOpAdaptor : public OpAdaptor<NotEqOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class NotEqOp : public mlir::Op<NotEqOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, OperandsAreSignlessIntegerOrFloatLike, mlir::OpTrait::SameTypeOperands, mlir::OpTrait::ResultsAreBoolLike, mlir::OpTrait::IsCommutative>
	{
		public:
		using Op::Op;
		using Adaptor = NotEqOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value lhs, mlir::Value rhs);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		mlir::LogicalResult verify();
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::GtOp
	//===----------------------------------------------------------------------===//

	class GtOp;

	class GtOpAdaptor : public OpAdaptor<GtOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class GtOp : public mlir::Op<GtOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, OperandsAreSignlessIntegerOrFloatLike, mlir::OpTrait::SameTypeOperands, mlir::OpTrait::ResultsAreBoolLike>
	{
		public:
		using Op::Op;
		using Adaptor = GtOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value lhs, mlir::Value rhs);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		mlir::LogicalResult verify();
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::GteOp
	//===----------------------------------------------------------------------===//

	class GteOp;

	class GteOpAdaptor : public OpAdaptor<GteOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class GteOp : public mlir::Op<GteOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, OperandsAreSignlessIntegerOrFloatLike, mlir::OpTrait::SameTypeOperands, mlir::OpTrait::ResultsAreBoolLike>
	{
		public:
		using Op::Op;
		using Adaptor = GteOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value lhs, mlir::Value rhs);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		mlir::LogicalResult verify();
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::LtOp
	//===----------------------------------------------------------------------===//

	class LtOp;

	class LtOpAdaptor : public OpAdaptor<LtOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class LtOp : public mlir::Op<LtOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, OperandsAreSignlessIntegerOrFloatLike, mlir::OpTrait::SameTypeOperands, mlir::OpTrait::ResultsAreBoolLike>
	{
		public:
		using Op::Op;
		using Adaptor = LtOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value lhs, mlir::Value rhs);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		mlir::LogicalResult verify();
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::LteOp
	//===----------------------------------------------------------------------===//

	class LteOp;

	class LteOpAdaptor : public OpAdaptor<LteOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class LteOp : public mlir::Op<LteOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, OperandsAreSignlessIntegerOrFloatLike, mlir::OpTrait::SameTypeOperands, mlir::OpTrait::ResultsAreBoolLike>
	{
		public:
		using Op::Op;
		using Adaptor = LteOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Value lhs, mlir::Value rhs);
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		mlir::LogicalResult verify();
		void print(mlir::OpAsmPrinter& printer);

		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::AddOp
	//===----------------------------------------------------------------------===//

	class AddOp;

	class AddOpAdaptor : public OpAdaptor<AddOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class AddOp : public mlir::Op<AddOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, mlir::OpTrait::IsCommutative>
	{
		public:
		using Op::Op;
		using Adaptor = AddOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		void print(mlir::OpAsmPrinter& printer);

		mlir::Type resultType();
		mlir::Value lhs();
		mlir::Value rhs();
	};

	//===----------------------------------------------------------------------===//
	// Modelica::SubOp
	//===----------------------------------------------------------------------===//

	class SubOp;

	class SubOpAdaptor : public OpAdaptor<SubOp>
	{
		public:
		using OpAdaptor::OpAdaptor;

		mlir::Value lhs();
		mlir::Value rhs();
	};

	class SubOp : public mlir::Op<SubOp, mlir::OpTrait::NOperands<2>::Impl, mlir::OpTrait::OneResult, mlir::OpTrait::IsCommutative>
	{
		public:
		using Op::Op;
		using Adaptor = SubOpAdaptor;

		static llvm::StringRef getOperationName();
		static void build(mlir::OpBuilder& builder, mlir::OperationState& state, mlir::Type resultType, mlir::Value lhs, mlir::Value rhs);
		void print(mlir::OpAsmPrinter& printer);

		mlir::Type resultType();
		mlir::Value lhs();
		mlir::Value rhs();
	};
}