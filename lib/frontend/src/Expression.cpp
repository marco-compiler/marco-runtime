#include <modelica/frontend/Expression.hpp>
#include <modelica/utils/IRange.hpp>

using namespace llvm;
using namespace modelica;
using namespace std;

struct LValueVisitor
{
	bool operator()(const Operation& obj) const { return obj.isLValue(); }
	bool operator()(const Constant& obj) const { return false; }
	bool operator()(const ReferenceAccess& obj) const { return true; }
	bool operator()(const Call& obj) const { return false; }
	bool operator()(const Tuple& obj) const { return false; }
};

Expression::Expression(Type type, Constant constant)
		: content(move(constant)), type(move(type))
{
}

Expression::Expression(Type type, ReferenceAccess access)
		: content(move(access)), type(move(type))
{
}

Expression::Expression(Type type, Call call)
		: content(move(call)), type(move(type))
{
}

Expression::Expression(Type type, Tuple tuple)
		: content(move(tuple)), type(move(type))
{
}

Expression::Expression(Type type, OperationKind kind, Operation::Container args)
		: content(Operation(kind, move(args))), type(move(type))
{
}

bool Expression::operator==(const Expression& other) const
{
	return type == other.type && content == other.content;
}

bool Expression::operator!=(const Expression& other) const
{
	return !(*this == other);
}

void Expression::dump() const { dump(outs(), 0); }

void Expression::dump(raw_ostream& os, size_t indents) const
{
	os.indent(indents);
	os << "type: ";
	getType().dump(os);
	os << "\n";

	visit([&](const auto& exp) { exp.dump(os, indents + 1); });
}

bool Expression::isLValue() const { return visit(LValueVisitor()); }

Type& Expression::getType() { return type; }

const Type& Expression::getType() const { return type; }

void Expression::setType(Type tp) { type = move(tp); }

Expression modelica::makeCall(Expression fun, llvm::ArrayRef<Expression> args)
{
	return Expression(Type::unknown(), Call(move(fun), move(args)));
}
