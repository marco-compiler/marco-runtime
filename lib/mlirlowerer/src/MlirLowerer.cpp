#include <llvm/ADT/SmallVector.h>
#include <mlir/Dialect/LLVMIR/LLVMDialect.h>
#include <mlir/Dialect/SCF/SCF.h>
#include <mlir/Dialect/StandardOps/IR/Ops.h>
#include <mlir/Pass/PassManager.h>
#include <modelica/mlirlowerer/LLVMLoweringPass.hpp>
#include <modelica/mlirlowerer/MlirLowerer.hpp>
#include <modelica/mlirlowerer/ModelicaDialect.hpp>
#include <modelica/mlirlowerer/ModelicaToStandard.hpp>

using namespace llvm;
using namespace mlir;
using namespace modelica;
using namespace std;

mlir::LogicalResult modelica::convertToLLVMDialect(mlir::MLIRContext* context, mlir::ModuleOp module)
{
	mlir::PassManager passManager(context);
	passManager.addPass(std::make_unique<ModelicaToStandardLoweringPass>());
	passManager.addPass(std::make_unique<LLVMLoweringPass>());
	return passManager.run(module);
}

mlir::Value Reference::operator*()
{
	if (isPtr)
		return builder->create<LoadOp>(builder->getUnknownLoc(), value, indexes);

	return value;
}

mlir::Value Reference::getValue() const
{
	return value;
}

bool Reference::isPointer() const
{
	return isPtr;
}

mlir::ValueRange Reference::getIndexes() const
{
	return indexes;
}

MlirLowerer::MlirLowerer(mlir::MLIRContext& context, bool x64)
		: builder(&context)
{
	// Check that the required dialects have been previously registered
	context.loadDialect<ModelicaDialect>();
	context.loadDialect<StandardOpsDialect>();
	context.loadDialect<scf::SCFDialect>();
	context.loadDialect<LLVM::LLVMDialect>();

	if (x64)
	{
		integerType = builder.getI64Type();
		floatType = builder.getF64Type();
	}
	else
	{
		integerType = builder.getI32Type();
		floatType = builder.getF32Type();
	}
}

mlir::OpBuilder& MlirLowerer::getOpBuilder()
{
	return builder;
}

mlir::Location MlirLowerer::loc(SourcePosition location)
{
	return builder.getFileLineColLoc(
			builder.getIdentifier(*location.file),
			location.line,
			location.column);
}

mlir::Value MlirLowerer::cast(mlir::Value value, mlir::Type destination)
{
	auto source = value.getType();

	if (source == destination)
		return value;

	mlir::Type sourceBase = source;

	if (sourceBase.isSignlessInteger())
	{
		if (destination.isF32() || destination.isF64())
			return builder.create<SIToFPOp>(builder.getUnknownLoc(), value, destination);
	}
	else if (source.isF32() || source.isF64())
	{
		if (destination.isSignlessInteger())
			return builder.create<FPToSIOp>(builder.getUnknownLoc(), value, destination);
	}

	assert(false && "Unsupported type conversion");
	return nullptr;
}

mlir::ModuleOp MlirLowerer::lower(llvm::ArrayRef<const modelica::ClassContainer> classes)
{
	mlir::ModuleOp module = ModuleOp::create(builder.getUnknownLoc());

	for (const auto& cls : classes)
	{
		auto* op = cls.visit([&](auto& obj) -> mlir::Operation* { return lower(obj); });

		if (op != nullptr)
			module.push_back(op);
	}

	return module;
}

mlir::Operation* MlirLowerer::lower(const modelica::Class& cls)
{
	return nullptr;
}

mlir::FuncOp MlirLowerer::lower(const modelica::Function& foo)
{
	// Create a scope in the symbol table to hold variable declarations.
	ScopedHashTableScope<StringRef, mlir::Value> varScope(symbolTable);

	auto location = loc(foo.getSourcePosition());

	SmallVector<llvm::StringRef, 3> argNames;
	SmallVector<mlir::Type, 3> argTypes;
	SmallVector<llvm::StringRef, 3> returnNames;
	SmallVector<mlir::Type, 3> returnTypes;

	for (const auto& member : foo.getResults())
	{
		if (member->getType().isScalar())
		{
			returnNames.emplace_back(member->getName());
			returnTypes.emplace_back(lower(member->getType()));
		}
		else
		{
			// If a result is an array, then it has to be placed among the arguments.
			// Otherwise, the callee should alloc the memory on the heap and thus
			// create a potential memory leak.

			argNames.emplace_back(member->getName());
			argTypes.emplace_back(lower(member->getType()));
		}
	}

	for (const auto& member : foo.getArgs())
	{
		argNames.emplace_back(member->getName());
		argTypes.emplace_back(lower(member->getType()));
	}

	auto functionType = builder.getFunctionType(argTypes, returnTypes);
	auto function = mlir::FuncOp::create(location, foo.getName(), functionType);

	// Start the body of the function.
	// In MLIR the entry block of the function is special: it must have the same
	// argument list as the function itself.
	auto &entryBlock = *function.addEntryBlock();

	// Declare all the function arguments in the symbol table
	for (const auto& pair : llvm::zip(argNames, entryBlock.getArguments())) {
		const auto& name = get<0>(pair);
		const auto& value = get<1>(pair);
		symbolTable.insert(name, value);
	}

	// Set the insertion point in the builder to the beginning of the function
	// body, it will be used throughout the codegen to create operations in this
	// function.
	builder.setInsertionPointToStart(&entryBlock);

	// Initialize members
	for (const auto& member : foo.getMembers())
		lower(member);

	// Emit the body of the function
	lower(foo.getAlgorithms()[0]);

	// Return statement
	std::vector<mlir::Value> results;

	for (const auto& name : returnNames)
	{
		auto ptr = symbolTable.lookup(name);
		auto value = builder.create<LoadOp>(builder.getUnknownLoc(), ptr);
		results.push_back(value);
	}

	builder.create<ReturnOp>(builder.getUnknownLoc(), results);
	return function;
}

mlir::Type MlirLowerer::lower(const modelica::Type& type)
{
	auto visitor = [&](auto& obj) -> mlir::Type
	{
		auto baseType = lower(obj);

		if (!type.isScalar())
		{
			const auto& dimensions = type.getDimensions();
			SmallVector<long, 3> shape(dimensions.begin(), dimensions.end());
			return mlir::MemRefType::get(shape, baseType);
		}

		return baseType;
	};

	return type.visit(visitor);
}

mlir::Type MlirLowerer::lower(const modelica::BuiltInType& type)
{
	switch (type)
	{
		case BuiltInType::None:
			return builder.getNoneType();
		case BuiltInType::Integer:
			return integerType;
		case BuiltInType::Float:
			return floatType;
		case BuiltInType::Boolean:
			return builder.getI1Type();
		default:
			assert(false && "Unexpected type");
	}
}

mlir::Type MlirLowerer::lower(const modelica::UserDefinedType& type)
{
	SmallVector<mlir::Type, 3> types;

	for (auto& subType : type)
		types.push_back(lower(subType));

	return builder.getTupleType(move(types));
}

void MlirLowerer::lower(const modelica::Member& member)
{
	if (member.getType().isScalar() || !member.isInput())
	{
		auto type = lower(member.getType());
		mlir::Value var;

		if (type.isa<ShapedType>())
		{
			auto shape = type.cast<ShapedType>().getShape();
			auto baseType = type.cast<ShapedType>().getElementType();
			var = builder.create<AllocaOp>(builder.getUnknownLoc(), MemRefType::get(shape, baseType));
			symbolTable.insert(member.getName(), var);
		}
		else
		{
			var = builder.create<AllocaOp>(builder.getUnknownLoc(), MemRefType::get({}, type));

			// Input variables already have an associated value. Copy it to the stack.
			if (member.isInput())
			{
				auto value = symbolTable.lookup(member.getName());
				builder.create<StoreOp>(builder.getUnknownLoc(), value, var);
			}

			symbolTable.insert(member.getName(), var);

			if (member.hasInitializer())
			{
				auto reference = lower<modelica::Expression>(member.getInitializer())[0];
				builder.create<StoreOp>(builder.getUnknownLoc(), *reference, var);
			}
		}
	}
}

void MlirLowerer::lower(const modelica::Algorithm& algorithm)
{
	for (const auto& statement : algorithm) {
		lower(statement);
	}
}

void MlirLowerer::lower(const modelica::Statement& statement)
{
	statement.visit([&](auto& obj) { lower(obj); });
}

void MlirLowerer::lower(const modelica::AssignmentStatement& statement)
{
	auto destinations = statement.getDestinations();
	auto values = lower<modelica::Expression>(statement.getExpression());
	assert(values.size() == destinations.size() && "Unequal number of destinations and results");

	for (auto pair : zip(destinations, values))
	{
		auto& ptr = lower<modelica::Expression>(get<0>(pair))[0];
		auto value = *get<1>(pair);
		builder.create<StoreOp>(builder.getUnknownLoc(), value, ptr.getValue(), ptr.getIndexes());
	}
}

void MlirLowerer::lower(const modelica::IfStatement& statement)
{
	ScopedHashTableScope<StringRef, mlir::Value> varScope(symbolTable);

	auto insertionPoint = builder.saveInsertionPoint();
	size_t blocks = statement.size();

	for (size_t i = 0; i < blocks; i++)
	{
		const auto& conditionalBlock = statement[i];
		auto condition = lower<modelica::Expression>(conditionalBlock.getCondition())[0];

		// The last conditional block can be at most an originally equivalent
		// "else" block, and thus doesn't need a lowered else block.
		bool elseBlock = i < blocks - 1;

		auto ifOp = builder.create<scf::IfOp>(builder.getUnknownLoc(), *condition, elseBlock);
		builder.setInsertionPointToStart(&ifOp.thenRegion().front());

		for (const auto& stmnt : conditionalBlock)
			lower(stmnt);

		if (elseBlock)
			builder.setInsertionPointToStart(&ifOp.elseRegion().front());
	}

	builder.restoreInsertionPoint(insertionPoint);
}

void MlirLowerer::lower(const modelica::ForStatement& statement)
{

}

void MlirLowerer::lower(const modelica::WhileStatement& statement)
{
	SmallVector<mlir::Type, 3> resultTypes;
	SmallVector<mlir::Value, 3> operands;
	auto whileOp = builder.create<scf::WhileOp>(builder.getUnknownLoc(), resultTypes, operands);

	auto insertionPoint = builder.saveInsertionPoint();

	builder.setInsertionPointToStart(&whileOp.before().front());
	auto condition = *lower<modelica::Expression>(statement.getCondition())[0];
	SmallVector<mlir::Value, 3> args;
	builder.create<scf::ConditionOp>(builder.getUnknownLoc(), condition, args);

	builder.setInsertionPointToStart(&whileOp.after().front());

	for (const auto& stmnt : statement)
		lower(stmnt);

	builder.restoreInsertionPoint(insertionPoint);
}

void MlirLowerer::lower(const modelica::WhenStatement& statement)
{

}

void MlirLowerer::lower(const modelica::BreakStatement& statement)
{

}

void MlirLowerer::lower(const modelica::ReturnStatement& statement)
{

}

template<>
MlirLowerer::Container<Reference> MlirLowerer::lower<Expression>(const modelica::Expression& expression)
{
	return expression.visit([&](auto& obj) {
		using type = decltype(obj);
		using deref = typename std::remove_reference<type>::type;
		using deconst = typename std::remove_const<deref>::type;
		return lower<deconst>(expression);
	});
}

template<>
MlirLowerer::Container<Reference> MlirLowerer::lower<modelica::Operation>(const modelica::Expression& expression)
{
	assert(expression.isA<modelica::Operation>());
	const auto& operation = expression.get<modelica::Operation>();
	auto kind = operation.getKind();
	mlir::Type resultType = lower(expression.getType());
	mlir::Location location = loc(expression.getLocation());

	if (kind == OperationKind::negate)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 1);
		mlir::Value result = builder.create<modelica::NegateOp>(location, args[0]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::add)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<modelica::AddOp>(location, args);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::subtract)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<modelica::SubOp>(location, args);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::multiply)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<modelica::MulOp>(location, args);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::divide)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<modelica::DivOp>(location, args);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::powerOf)
	{
		// TODO
		/*
		mlir::Value base = cast(*lower<modelica::Expression>(operation[0])[0], floatType);
		mlir::Value exponent = cast(*lower<modelica::Expression>(operation[1])[0], floatType);
		mlir::Value result = builder.create<>(location, base, exponent);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
		 */
		return { Reference(builder, nullptr, false) };
	}

	if (kind == OperationKind::equal)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 2);
		mlir::Value result = builder.create<EqOp>(location, args[0], args[1]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::different)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 2);
		mlir::Value result = builder.create<NotEqOp>(location, args[0], args[1]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::greater)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 2);
		mlir::Value result = builder.create<GtOp>(location, args[0], args[1]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::greaterEqual)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 2);
		mlir::Value result = builder.create<GteOp>(location, args[0], args[1]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::less)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 2);
		mlir::Value result = builder.create<LtOp>(location, args[0], args[1]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::lessEqual)
	{
		auto args = lowerOperationArgs(operation);
		assert(args.size() == 2);
		mlir::Value result = builder.create<LteOp>(location, args[0], args[1]);
		result = cast(result, resultType);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::ifelse)
	{
		mlir::Value condition = *lower<modelica::Expression>(operation[0])[0];
		mlir::Value trueValue = cast(*lower<modelica::Expression>(operation[1])[0], resultType);
		mlir::Value falseValue = cast(*lower<modelica::Expression>(operation[2])[0], resultType);
		mlir::Value result = builder.create<SelectOp>(location, condition, trueValue, falseValue);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::land)
	{
		mlir::Value lhs = *lower<modelica::Expression>(operation[0])[0];
		mlir::Value rhs = *lower<modelica::Expression>(operation[1])[0];

		unsigned int lhsBitWidth = lhs.getType().getIntOrFloatBitWidth();
		unsigned int rhsBitWidth = lhs.getType().getIntOrFloatBitWidth();

		if (lhsBitWidth < rhsBitWidth)
			lhs = builder.create<SignExtendIOp>(location, lhs, rhs.getType());
		else if (lhsBitWidth > rhsBitWidth)
			rhs = builder.create<SignExtendIOp>(location, rhs, lhs.getType());

		mlir::Value result = builder.create<AndOp>(location, lhs, rhs);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::lor)
	{
		mlir::Value lhs = *lower<modelica::Expression>(operation[0])[0];
		mlir::Value rhs = *lower<modelica::Expression>(operation[1])[0];

		unsigned int lhsBitWidth = lhs.getType().getIntOrFloatBitWidth();
		unsigned int rhsBitWidth = lhs.getType().getIntOrFloatBitWidth();

		if (lhsBitWidth < rhsBitWidth)
			lhs = builder.create<SignExtendIOp>(location, lhs, rhs.getType());
		else if (lhsBitWidth > rhsBitWidth)
			rhs = builder.create<SignExtendIOp>(location, rhs, lhs.getType());

		mlir::Value result = builder.create<OrOp>(location, lhs, rhs);
		return { Reference(builder, result, false) };
	}

	if (kind == OperationKind::subscription)
	{
		auto var = lower<modelica::Expression>(operation[0])[0];
		assert(var.isPointer());
		SmallVector<mlir::Value, 3> indexes;

		for (size_t i = 1; i < operation.argumentsCount(); i++)
		{
			auto subscript = *lower<modelica::Expression>(operation[i])[0];
			auto index = builder.create<IndexCastOp>(builder.getUnknownLoc(), subscript, builder.getIndexType());
			indexes.push_back(index);
		}

		return { Reference(builder, var.getValue(), true, indexes) };
	}

	if (kind == OperationKind::memberLookup)
	{
		// TODO
		return { Reference(builder, nullptr, false) };
	}

	assert(false && "Unexpected operation");
	return { Reference(builder, nullptr, false) };
}

template<>
MlirLowerer::Container<Reference> MlirLowerer::lower<modelica::Constant>(const modelica::Expression& expression)
{
	assert(expression.isA<modelica::Constant>());
	const auto& constant = expression.get<modelica::Constant>();

	auto value = builder.create<ConstantOp>(
			loc(expression.getLocation()),
			constantToType(constant),
			constant.visit([&](const auto& obj) { return getAttribute(obj); }));

	return { Reference(builder, value, false) };
}

template<>
MlirLowerer::Container<Reference> MlirLowerer::lower<modelica::ReferenceAccess>(const modelica::Expression& expression)
{
	assert(expression.isA<modelica::ReferenceAccess>());
	const auto& reference = expression.get<modelica::ReferenceAccess>();
	auto ptr = symbolTable.lookup(reference.getName());
	return { Reference(builder, ptr, true) };
}

template<>
MlirLowerer::Container<Reference> MlirLowerer::lower<modelica::Call>(const modelica::Expression& expression)
{
	assert(expression.isA<modelica::Call>());
	const auto& call = expression.get<modelica::Call>();
	const auto& function = call.getFunction();

	SmallVector<mlir::Value, 3> args;

	for (const auto& arg : call)
	{
		auto& reference = lower<modelica::Expression>(arg)[0];
		args.push_back(*reference);
	}

	auto op = builder.create<CallOp>(
			loc(expression.getLocation()),
			function.get<ReferenceAccess>().getName(),
			lower(function.getType()),
			args);

	Container<Reference> results;

	for (auto result : op.getResults())
		results.emplace_back(builder, result, false);

	return results;
}

template<>
MlirLowerer::Container<Reference> MlirLowerer::lower<modelica::Tuple>(const modelica::Expression& expression)
{
	assert(expression.isA<modelica::Tuple>());
	const auto& tuple = expression.get<modelica::Tuple>();
	Container<Reference> result;

	for (auto& exp : tuple)
	{
		auto values = lower<modelica::Expression>(expression);

		// The only way to have multiple returns is to call a function, but this
		// is forbidden in a tuple declaration. In fact, a tuple is just a
		// container of references.
		assert(values.size() == 1);
		result.emplace_back(values[0]);
	}

	return result;
}

MlirLowerer::Container<mlir::Value> MlirLowerer::lowerOperationArgs(const modelica::Operation& operation)
{
	SmallVector<mlir::Value, 3> args;
	SmallVector<mlir::Value, 3> castedArgs;
	bool containsFloat = false;

	for (const auto& arg : operation)
	{
		// For now, we only support operation between built-in types.
		// In future, this should be extended to support struct types.
		assert(arg.getType().isA<modelica::BuiltInType>());

		if (arg.getType().get<BuiltInType>() == BuiltInType::Float)
			containsFloat = true;

		args.push_back(*lower<modelica::Expression>(arg)[0]);
	}

	// Convert the arguments to a common type.
	// If any of the arguments is a float, also the others must be floats, in
	// order to preserve correctness. If a value is a boolean, it is first
	// extended to an integer in first place, and then, if needed, to a float.

	for (const auto& arg : args)
	{
		mlir::Value castedArg = arg;

		if (castedArg.getType().isInteger(1))
			castedArg = builder.create<SignExtendIOp>(arg.getLoc(), castedArg, integerType);

		if (containsFloat && castedArg.getType().isSignlessInteger())
			castedArg = builder.create<SIToFPOp>(arg.getLoc(), castedArg, floatType);

		castedArgs.push_back(castedArg);
	}

	return castedArgs;
}
