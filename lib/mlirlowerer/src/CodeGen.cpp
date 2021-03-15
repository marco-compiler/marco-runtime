#include <llvm/ADT/SmallVector.h>
#include <mlir/Conversion/Passes.h>
#include <mlir/Dialect/Affine/IR/AffineOps.h>
#include <mlir/Dialect/LLVMIR/LLVMDialect.h>
#include <mlir/Dialect/Linalg/Passes.h>
#include <mlir/Dialect/StandardOps/IR/Ops.h>
#include <mlir/Dialect/StandardOps/Transforms/Passes.h>
#include <mlir/Dialect/Vector/VectorOps.h>
#include <mlir/IR/Verifier.h>
#include <mlir/Pass/PassManager.h>
#include <mlir/Transforms/Passes.h>
#include <modelica/frontend/AST.h>
#include <modelica/mlirlowerer/CodeGen.h>
#include <modelica/mlirlowerer/ModelicaDialect.h>
#include <modelica/mlirlowerer/passes/LowerToLLVM.h>
#include <modelica/mlirlowerer/passes/ModelicaConversionPass.h>

using namespace modelica;
using namespace std;

mlir::LogicalResult modelica::convertToLLVMDialect(mlir::MLIRContext* context, mlir::ModuleOp module, ModelicaOptions options)
{
	mlir::PassManager passManager(context);

	ModelicaConversionOptions modelicaToLLVMOptions;
	passManager.addPass(createModelicaConversionPass(modelicaToLLVMOptions));

	//passManager.addNestedPass<mlir::FuncOp>(mlir::createConvertSCFToOpenMPPass());
	passManager.addPass(mlir::createLowerToCFGPass());
	passManager.addPass(createLLVMLoweringPass());

	return passManager.run(module);
}

Reference::Reference() : builder(nullptr), value(nullptr), reader(nullptr)
{
}

Reference::Reference(ModelicaBuilder* builder,
										 mlir::Value value,
										 bool initialized,
										 std::function<mlir::Value(ModelicaBuilder*, mlir::Value)> reader)
		: builder(builder),
			value(std::move(value)),
			initialized(initialized),
			reader(std::move(reader))
{
}

mlir::Value Reference::operator*()
{
	return reader(builder, value);
}

mlir::Value Reference::getReference() const
{
	return value;
}

bool Reference::isInitialized() const
{
	return initialized;
}

Reference Reference::ssa(ModelicaBuilder* builder, mlir::Value value)
{
	return Reference(
			builder, value, true,
			[](mlir::OpBuilder* builder, mlir::Value value) -> mlir::Value
			{
				return value;
			});
}

Reference Reference::memory(ModelicaBuilder* builder, mlir::Value value, bool initialized)
{
	return Reference(
			builder, value, initialized,
			[](mlir::OpBuilder* builder, mlir::Value value) -> mlir::Value
			{
				return builder->create<LoadOp>(value.getLoc(), value);
			});
}

MLIRLowerer::MLIRLowerer(mlir::MLIRContext& context, ModelicaOptions options)
		: builder(&context, options.getBitWidth()), options(move(options))
{
	context.loadDialect<ModelicaDialect>();
	context.loadDialect<mlir::StandardOpsDialect>();
}

mlir::Location MLIRLowerer::loc(SourcePosition location)
{
	return mlir::FileLineColLoc::get(
			builder.getIdentifier(*location.file),
			location.line,
			location.column);
}

llvm::Optional<mlir::ModuleOp> MLIRLowerer::lower(llvm::ArrayRef<const ClassContainer> classes)
{
	mlir::ModuleOp module = mlir::ModuleOp::create(builder.getUnknownLoc());

	for (const auto& cls : classes)
	{
		auto* op = cls.visit([&](auto& obj) -> mlir::Operation* { return lower(obj); });

		if (op != nullptr)
			module.push_back(op);
	}

	if (failed(mlir::verify(module)))
		return llvm::None;

	return module;
}

mlir::Operation* MLIRLowerer::lower(const Class& cls)
{
	return nullptr;
}

mlir::FuncOp MLIRLowerer::lower(const Function& foo)
{
	// Create a scope in the symbol table to hold variable declarations.
	llvm::ScopedHashTableScope<mlir::StringRef, Reference> varScope(symbolTable);

	auto location = loc(foo.getLocation());

	llvm::SmallVector<llvm::StringRef, 3> argNames;
	llvm::SmallVector<mlir::Type, 3> argTypes;

	for (const auto& member : foo.getArgs())
	{
		argNames.emplace_back(member->getName());
		argTypes.emplace_back(lower(member->getType()));
	}

	llvm::SmallVector<llvm::StringRef, 3> returnNames;
	llvm::SmallVector<mlir::Type, 3> returnTypes;
	auto outputMembers = foo.getResults();

	for (const auto& member : outputMembers)
	{
		mlir::Type type = lower(member->getType());

		if (auto pointerType = type.dyn_cast<PointerType>(); pointerType && pointerType.hasConstantShape())
		{
			argNames.emplace_back(member->getName());
			argTypes.emplace_back(type);
			continue;
		}

		if (type.isa<PointerType>())
			type = builder.getPointerType(true, type.cast<PointerType>().getElementType(), type.cast<PointerType>().getShape());

		returnNames.emplace_back(member->getName());
		returnTypes.emplace_back(type);
	}

	auto functionType = builder.getFunctionType(argTypes, returnTypes);
	auto function = mlir::FuncOp::create(location, foo.getName(), functionType);

	// If the function doesn't have a body, it means it is just a declaration
	if (foo.getAlgorithms().empty())
		return function;

	// Start the body of the function.
	// In MLIR the entry block of the function is special: it must have the same
	// argument list as the function itself.
	auto &entryBlock = *function.addEntryBlock();

	// Declare all the function arguments in the symbol table
	for (const auto& pair : llvm::zip(argNames, entryBlock.getArguments())) {
		const auto& name = get<0>(pair);
		const auto& value = get<1>(pair);
		symbolTable.insert(name, Reference::ssa(&builder, value));
	}

	// Set the insertion point in the builder to the beginning of the function
	// body, it will be used throughout the codegen to create operations in this
	// function.
	builder.setInsertionPointToStart(&entryBlock);

	// Initialize members
	for (const auto& member : foo.getMembers())
		lower(*member);

	// Emit the body of the function
	const auto& algorithm = foo.getAlgorithms()[0];

	// Create the variable to be checked for an early return
	auto algorithmLocation = loc(algorithm->getLocation());
	mlir::Value returnCondition = builder.create<AllocaOp>(algorithmLocation, builder.getBooleanType());
	mlir::Value falseValue = builder.create<ConstantOp>(algorithmLocation, builder.getBooleanAttribute(false));
	builder.create<StoreOp>(algorithmLocation, falseValue, returnCondition);
	symbolTable.insert(algorithm->getReturnCheckName(), Reference::memory(&builder, returnCondition, true));

	// Lower the statements
	lower(*foo.getAlgorithms()[0]);

	// Return statement
	llvm::SmallVector<mlir::Value, 3> results;

	for (const auto& name : returnNames)
	{
		auto ptr = symbolTable.lookup(name);
		results.push_back(*ptr);
	}

	builder.create<mlir::ReturnOp>(location, results);
	return function;
}

mlir::Type MLIRLowerer::lower(const Type& type)
{
	auto visitor = [&](auto& obj) -> mlir::Type
	{
		auto baseType = lower(obj);

		if (!type.isScalar())
		{
			const auto& dimensions = type.getDimensions();
			llvm::SmallVector<long, 3> shape;

			for (const auto& dimension : type.getDimensions())
			{
				if (dimension.isDynamic())
					shape.emplace_back(-1);
				else
					shape.emplace_back(dimension.getNumericSize());
			}

			return builder.getPointerType(false, baseType, shape);
		}

		return baseType;
	};

	return type.visit(visitor);
}

mlir::Type MLIRLowerer::lower(const BuiltInType& type)
{
	switch (type)
	{
		case BuiltInType::None:
			return builder.getNoneType();
		case BuiltInType::Integer:
			return builder.getIntegerType();
		case BuiltInType::Float:
			return builder.getRealType();
		case BuiltInType::Boolean:
			return builder.getBooleanType();
		default:
			assert(false && "Unexpected type");
			return builder.getNoneType();
	}
}

mlir::Type MLIRLowerer::lower(const UserDefinedType& type)
{
	llvm::SmallVector<mlir::Type, 3> types;

	for (auto& subType : type)
		types.push_back(lower(subType));

	return builder.getTupleType(move(types));
}

/**
 * Lower a member of a function.
 * If the size of the element can be determined statically, then it is
 * allocated on the stack. If not, it will be allocated when it will be
 * initialized. Output arrays with constant shape are considered to be
 * already allocated by the caller and passed as argument to the function.
 * Input members are ignored because they are supposed to be unmodifiable
 * as per the Modelica standard.
 *
 * 	                 INPUT        OUTPUT        PROTECTED
 * scalar              -          stack           stack
 * array ranked        -            -             stack
 * array unranked      -          heap            stack
 */
void MLIRLowerer::lower(const Member& member)
{
	auto location = loc(member.getLocation());

	// Input values are supposed to be read-only by the Modelica standard,
	// thus they don't need to be allocated on the stack for modifications.
	if (member.isInput())
		return;

	mlir::Type type = lower(member.getType());

	if (!member.isOutput() || !type.isa<PointerType>() || !type.cast<PointerType>().hasConstantShape())
	{
		if (type.isa<PointerType>())
		{
			auto pointerType = type.cast<PointerType>();

			if (member.isOutput())
				type = builder.getPointerType(true, pointerType.getElementType(), pointerType.getShape());
		}

		mlir::Value ptr = builder.create<AllocaOp>(location, type);
		bool initialized = false;

		if (type.isa<PointerType>())
		{
			auto pointerType = type.cast<PointerType>();
			llvm::SmallVector<mlir::Value, 3> sizes;

			for (const auto& dimension : member.getType().getDimensions())
				if (dimension.hasExpression())
				{
					mlir::Value size = *lower<Expression>(dimension.getExpression())[0];
					size = builder.create<CastOp>(location, size, builder.getIndexType());
					sizes.push_back(size);
				}

			if (sizes.size() == pointerType.getDynamicDimensions())
			{
				// All the dynamic dimensions have an expression to determine their values.
				// So we can instantiate the array.

				if (pointerType.isOnHeap())
				{
					mlir::Value var = builder.create<AllocOp>(location, pointerType.getElementType(), pointerType.getShape(), sizes);
					builder.create<StoreOp>(location, var, ptr);
				}
				else
				{
					mlir::Value var = builder.create<AllocaOp>(location, pointerType.getElementType(), pointerType.getShape(), sizes);
					builder.create<StoreOp>(location, var, ptr);
				}

				initialized = true;
			}
		}

		symbolTable.insert(member.getName(), Reference::memory(&builder, ptr, initialized));
	}

	/*
	if (member.hasInitializer())
	{
		mlir::Value reference = symbolTable.lookup(member.getName()).getReference();
		mlir::Value value = lower<modelica::Expression>(member.getInitializer())[0].getReference();
		builder.create<AssignmentOp>(loc(member.getInitializer().getLocation()), value, reference);
	}
	 */
}

void MLIRLowerer::lower(const Algorithm& algorithm)
{
	for (const auto& statement : algorithm)
		lower(*statement);
}

void MLIRLowerer::lower(const Statement& statement)
{
	statement.visit([&](auto& obj) { lower(obj); });
}

void MLIRLowerer::lower(const AssignmentStatement& statement)
{
	auto location = loc(statement.getLocation());
	auto destinations = statement.getDestinations();
	auto values = lower<Expression>(statement.getExpression());
	assert(values.size() == destinations.size() && "Unequal number of destinations and results");

	for (auto pair : zip(destinations, values))
	{
		auto destination = lower<Expression>(get<0>(pair))[0];
		auto value = get<1>(pair);

		auto destinationPointer = destination.getReference().getType().cast<PointerType>();

		if (destinationPointer.getElementType().isa<PointerType>())
		{
			if (destination.isInitialized())
				builder.create<AssignmentOp>(location, *value, *destination);
			else
			{
				auto pointer = destinationPointer.getElementType().cast<PointerType>();

				// Copy source on stack
				// Save the descriptor of the new copy into the destination using StoreOp

				mlir::Value copy = builder.create<ArrayCopyOp>(location, *value, pointer.isOnHeap());
				builder.create<StoreOp>(location, copy, destination.getReference());
			}
		}
		else
		{
			builder.create<AssignmentOp>(location, *value, destination.getReference());
		}
	}
}

void MLIRLowerer::lower(const IfStatement& statement)
{
	// Each conditional blocks creates an If operation, but we need to keep
	// track of the first one in order to restore the insertion point right
	// after that when we have finished to lower all the blocks.
	mlir::Operation* firstOp = nullptr;

	size_t blocks = statement.size();

	for (size_t i = 0; i < blocks; i++)
	{
		llvm::ScopedHashTableScope<mlir::StringRef, Reference> varScope(symbolTable);
		const auto& conditionalBlock = statement[i];
		auto condition = lower<Expression>(conditionalBlock.getCondition())[0];

		// The last conditional block can be at most an originally equivalent
		// "else" block, and thus doesn't need a lowered else block.
		bool elseBlock = i < blocks - 1;

		auto ifOp = builder.create<IfOp>(loc(statement.getLocation()), *condition, elseBlock);

		if (firstOp == nullptr)
			firstOp = ifOp;

		// "Then" block
		builder.setInsertionPointToStart(&ifOp.thenRegion().front());

		for (const auto& stmnt : conditionalBlock)
			lower(*stmnt);

		builder.create<YieldOp>(loc(statement.getLocation()));

		if (i > 0)
		{
			builder.setInsertionPointAfter(ifOp);
			builder.create<YieldOp>(loc(statement.getLocation()));
		}

		// The next conditional blocks will be placed as new If operations
		// nested inside the "else" block.
		if (elseBlock)
			builder.setInsertionPointToStart(&ifOp.elseRegion().front());
	}

	builder.setInsertionPointAfter(firstOp);
}

void MLIRLowerer::lower(const ForStatement& statement)
{
	llvm::ScopedHashTableScope<mlir::StringRef, Reference> varScope(symbolTable);

	auto location = loc(statement.getLocation());

	// Variable to be set when calling "break"
	mlir::Value breakCondition = builder.create<AllocaOp>(location, builder.getBooleanType());
	mlir::Value falseValue = builder.create<ConstantOp>(location, builder.getBooleanAttribute(false));
	builder.create<StoreOp>(location, falseValue, breakCondition);
	symbolTable.insert(statement.getBreakCheckName(),
			Reference::memory(&builder, breakCondition, true));

	// Variable to be set when calling "return"
	mlir::Value returnCondition = symbolTable.lookup(statement.getReturnCheckName()).getReference();

	const auto& induction = statement.getInduction();

	mlir::Value lowerBound = *lower<Expression>(induction.getBegin())[0];
	lowerBound = builder.create<CastOp>(lowerBound.getLoc(), lowerBound, builder.getIndexType());

	auto forOp = builder.create<ForOp>(location, breakCondition, returnCondition, lowerBound);

	{
		// Check the loop condition
		llvm::ScopedHashTableScope<mlir::StringRef, Reference> scope(symbolTable);
		symbolTable.insert(induction.getName(), Reference::ssa(&builder, forOp.condition().front().getArgument(0)));

		builder.setInsertionPointToStart(&forOp.condition().front());

		mlir::Value upperBound = *lower<Expression>(induction.getEnd())[0];
		upperBound = builder.create<CastOp>(lowerBound.getLoc(), upperBound, builder.getIndexType());

		mlir::Value condition = builder.create<LteOp>(location, builder.getBooleanType(), forOp.condition().front().getArgument(0), upperBound);
		builder.create<ConditionOp>(location, condition, *symbolTable.lookup(induction.getName()));
	}

	{
		// Body
		llvm::ScopedHashTableScope<mlir::StringRef, Reference> scope(symbolTable);
		symbolTable.insert(induction.getName(), Reference::ssa(&builder, forOp.body().front().getArgument(0)));

		builder.setInsertionPointToStart(&forOp.body().front());

		for (const auto& stmnt : statement)
			lower(*stmnt);

		builder.create<YieldOp>(location, *symbolTable.lookup(induction.getName()));
	}

	{
		// Step
		llvm::ScopedHashTableScope<mlir::StringRef, Reference> scope(symbolTable);
		symbolTable.insert(induction.getName(), Reference::ssa(&builder, forOp.step().front().getArgument(0)));

		builder.setInsertionPointToStart(&forOp.step().front());

		mlir::Value step = builder.create<mlir::ConstantOp>(location, builder.getIndexAttribute(1));
		mlir::Value incremented = builder.create<mlir::AddIOp>(location, *symbolTable.lookup(induction.getName()), step);
		builder.create<YieldOp>(location, incremented);
	}

	builder.setInsertionPointAfter(forOp);
}

void MLIRLowerer::lower(const WhileStatement& statement)
{
	auto location = loc(statement.getLocation());

	// Variable to be set when calling "break"
	mlir::Value breakCondition = builder.create<AllocaOp>(location, builder.getBooleanType());
	mlir::Value falseValue = builder.create<ConstantOp>(location, builder.getBooleanAttribute(false));
	builder.create<StoreOp>(location, falseValue, breakCondition);
	symbolTable.insert(statement.getBreakCheckName(),
			Reference::memory(&builder, breakCondition, true));

	// Variable to be set when calling "return"
	mlir::Value returnCondition = symbolTable.lookup(statement.getReturnCheckName()).getReference();

	// Create the operation
	auto whileOp = builder.create<WhileOp>(location, breakCondition, returnCondition);

	{
		// Condition
		llvm::ScopedHashTableScope<mlir::StringRef, Reference> varScope(symbolTable);
		mlir::Block* conditionBlock = &whileOp.condition().front();
		builder.setInsertionPointToStart(conditionBlock);
		const auto& condition = statement.getCondition();

		builder.create<ConditionOp>(
				loc(condition.getLocation()),
				*lower<Expression>(condition)[0]);
	}

	{
		// Body
		llvm::ScopedHashTableScope<mlir::StringRef, Reference> varScope(symbolTable);
		builder.setInsertionPointToStart(&whileOp.body().front());

		for (const auto& stmnt : statement)
			lower(*stmnt);

		builder.create<YieldOp>(location);
	}

	// Keep populating after the while operation
	builder.setInsertionPointAfter(whileOp);
}

void MLIRLowerer::lower(const WhenStatement& statement)
{

}

void MLIRLowerer::lower(const BreakStatement& statement)
{
	assert(false && "Break statement encountered. BreakRemovingPass may have not been run before lowering the AST.");
}

void MLIRLowerer::lower(const ReturnStatement& statement)
{
	assert(false && "Return statement encountered. ReturnRemovingPass may have not been run before lowering the AST.");
}

template<>
MLIRLowerer::Container<Reference> MLIRLowerer::lower<Expression>(const Expression& expression)
{
	return expression.visit([&](auto& obj) {
		using type = decltype(obj);
		using deref = typename std::remove_reference<type>::type;
		using deconst = typename std::remove_const<deref>::type;
		return lower<deconst>(expression);
	});
}

template<>
MLIRLowerer::Container<Reference> MLIRLowerer::lower<Operation>(const Expression& expression)
{
	assert(expression.isA<Operation>());
	const auto& operation = expression.get<Operation>();
	auto kind = operation.getKind();
	mlir::Type resultType = lower(expression.getType());
	mlir::Location location = loc(expression.getLocation());

	if (kind == OperationKind::negate)
	{
		auto arg = lower<Expression>(operation[0])[0].getReference();
		mlir::Value result = builder.create<NotOp>(location, arg);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::add)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = foldBinaryOperation(
				args,
				[&](mlir::Value lhs, mlir::Value rhs) -> mlir::Value
				{
					return builder.create<AddOp>(location, resultType, lhs, rhs);
				});

		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::subtract)
	{
		auto args = lowerOperationArgs(operation);

		if (args.size() == 1)
		{
			// Special case for sign change (i.e "-x").
			// TODO
			// In future, when all the project will rely on MLIR, a different
			// operation in the frontend should be created for this purpose.

			mlir::Value result = builder.create<NegateOp>(location, resultType, args[0]);
			return { Reference::ssa(&builder, result) };
		}
		else
		{
			mlir::Value result = foldBinaryOperation(
					args,
					[&](mlir::Value lhs, mlir::Value rhs) -> mlir::Value
					{
						return builder.create<SubOp>(location, resultType, lhs, rhs);
					});

			return { Reference::ssa(&builder, result) };
		}
	}

	if (kind == OperationKind::multiply)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = foldBinaryOperation(
				args,
				[&](mlir::Value lhs, mlir::Value rhs) -> mlir::Value
				{
					return builder.create<MulOp>(location, resultType, lhs, rhs);
				});

		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::divide)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = foldBinaryOperation(
				args,
				[&](mlir::Value lhs, mlir::Value rhs) -> mlir::Value
				{
					return builder.create<DivOp>(location, resultType, lhs, rhs);
				});

		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::powerOf)
	{
		mlir::Value base = *lower<Expression>(operation[0])[0];
		mlir::Value exponent = *lower<Expression>(operation[1])[0];

		if (base.getType().isa<PointerType>())
		{
			exponent = builder.create<CastOp>(base.getLoc(), exponent, builder.getIntegerType());
		}
		else
		{
			base = builder.create<CastOp>(base.getLoc(), base, builder.getRealType());
			exponent = builder.create<CastOp>(base.getLoc(), exponent, builder.getRealType());
		}

		mlir::Value result = builder.create<PowOp>(location, resultType, base, exponent);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::equal)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<EqOp>(location, resultType, args[0], args[1]);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::different)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<NotEqOp>(location, resultType, args[0], args[1]);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::greater)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<GtOp>(location, resultType, args[0], args[1]);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::greaterEqual)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<GteOp>(location, resultType, args[0], args[1]);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::less)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<LtOp>(location, resultType, args[0], args[1]);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::lessEqual)
	{
		auto args = lowerOperationArgs(operation);
		mlir::Value result = builder.create<LteOp>(location, resultType, args[0], args[1]);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::ifelse)
	{
		mlir::Value condition = *lower<Expression>(operation[0])[0];

		mlir::Value trueValue = *lower<Expression>(operation[1])[0];
		trueValue = builder.create<CastOp>(trueValue.getLoc(), trueValue, resultType);

		mlir::Value falseValue = *lower<Expression>(operation[2])[0];
		falseValue = builder.create<CastOp>(falseValue.getLoc(), falseValue, resultType);

		mlir::Value result = builder.create<mlir::SelectOp>(location, condition, trueValue, falseValue);
		result = builder.create<CastOp>(result.getLoc(), result, resultType);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::land)
	{
		mlir::Value lhs = *lower<Expression>(operation[0])[0];
		mlir::Value rhs = *lower<Expression>(operation[1])[0];

		mlir::Value result = builder.create<AndOp>(location, resultType, lhs, rhs);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::lor)
	{
		mlir::Value lhs = *lower<Expression>(operation[0])[0];
		mlir::Value rhs = *lower<Expression>(operation[1])[0];

		mlir::Value result = builder.create<OrOp>(location, resultType, lhs, rhs);
		return { Reference::ssa(&builder, result) };
	}

	if (kind == OperationKind::subscription)
	{
		auto buffer = *lower<Expression>(operation[0])[0];
		assert(buffer.getType().isa<PointerType>());

		llvm::SmallVector<mlir::Value, 3> indexes;

		for (size_t i = 1; i < operation.argumentsCount(); i++)
		{
			auto subscript = *lower<Expression>(operation[i])[0];
			mlir::Value index = builder.create<CastOp>(subscript.getLoc(), subscript, builder.getIndexType());
			indexes.push_back(index);
		}

		mlir::Value result = builder.create<SubscriptionOp>(location, buffer, indexes);
		return { Reference::memory(&builder, result, true) };
	}

	if (kind == OperationKind::memberLookup)
	{
		// TODO
		return { Reference::ssa(&builder, nullptr) };
	}

	assert(false && "Unexpected operation");
	return { Reference::ssa(&builder, nullptr) };
}

template<>
MLIRLowerer::Container<Reference> MLIRLowerer::lower<Constant>(const Expression& expression)
{
	assert(expression.isA<Constant>());
	const auto& constant = expression.get<Constant>();

	auto value = builder.create<ConstantOp>(
			loc(expression.getLocation()),
			constant.visit([&](const auto& obj) { return getAttribute(obj); }));

	return { Reference::ssa(&builder, value) };
}

template<>
MLIRLowerer::Container<Reference> MLIRLowerer::lower<ReferenceAccess>(const Expression& expression)
{
	assert(expression.isA<ReferenceAccess>());
	const auto& reference = expression.get<ReferenceAccess>();
	return { symbolTable.lookup(reference.getName()) };
}

template<>
MLIRLowerer::Container<Reference> MLIRLowerer::lower<Call>(const Expression& expression)
{
	assert(expression.isA<Call>());
	const auto& call = expression.get<Call>();
	const auto& function = call.getFunction();
	mlir::Location location = loc(expression.getLocation());

	llvm::SmallVector<mlir::Value, 3> args;

	for (const auto& arg : call)
	{
		auto& reference = lower<Expression>(arg)[0];
		args.push_back(*reference);
	}

	auto& functionName = function.get<ReferenceAccess>().getName();

	Container<Reference> results;

	if (functionName == "ndims")
	{
		assert(args.size() == 1);
		mlir::Type resultType = lower(expression.getType());
		mlir::Value result = builder.create<NDimsOp>(location, resultType, args[0]);
		results.emplace_back(Reference::ssa(&builder, result));
	}
	else if (functionName == "size")
	{
		assert(args.size() == 1 || args.size() == 2);
		mlir::Type resultType = lower(expression.getType());

		if (args.size() == 1)
		{
			mlir::Value result = builder.create<SizeOp>(location, resultType, args[0]);
			results.emplace_back(Reference::ssa(&builder, result));
		}
		else if (args.size() == 2)
		{
			mlir::Value oneValue = builder.create<ConstantOp>(location, builder.getIntegerAttribute(1));
			mlir::Value index = builder.create<SubOp>(location, builder.getIndexType(), args[1], oneValue);
			mlir::Value result = builder.create<SizeOp>(location, resultType, args[0], index);
			results.emplace_back(Reference::ssa(&builder, result));
		}
	}
	else
	{
		auto op = builder.create<CallOp>(
				loc(expression.getLocation()),
				function.get<ReferenceAccess>().getName(),
				lower(function.getType()),
				args);

		for (auto result : op.getResults())
			results.emplace_back(Reference::ssa(&builder, result));
	}

	return results;
}

template<>
MLIRLowerer::Container<Reference> MLIRLowerer::lower<Tuple>(const Expression& expression)
{
	assert(expression.isA<Tuple>());
	const auto& tuple = expression.get<Tuple>();
	Container<Reference> result;

	for (auto& exp : tuple)
	{
		auto values = lower<Expression>(expression);

		// The only way to have multiple returns is to call a function, but this
		// is forbidden in a tuple declaration. In fact, a tuple is just a
		// container of references.
		assert(values.size() == 1);
		result.emplace_back(values[0]);
	}

	return result;
}

mlir::Value MLIRLowerer::foldBinaryOperation(llvm::ArrayRef<mlir::Value> args, std::function<mlir::Value(mlir::Value, mlir::Value)> callback)
{
	assert(args.size() >= 2);
	mlir::Value result = callback(args[0], args[1]);

	for (size_t i = 2, e = args.size(); i < e; ++i)
		result = callback(result, args[i]);

	return result;
}

MLIRLowerer::Container<mlir::Value> MLIRLowerer::lowerOperationArgs(const Operation& operation)
{
	llvm::SmallVector<mlir::Value, 3> args;

	for (const auto& arg : operation)
		args.push_back(*lower<Expression>(arg)[0]);

	return args;
}