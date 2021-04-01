#include <algorithm>
#include <modelica/frontend/AST.h>

using namespace modelica;
using namespace frontend;

template<typename T>
using Container = Function::Container<T>;

Function::Function(
		SourcePosition location,
		std::string name,
		bool pure,
		llvm::ArrayRef<Member> members,
		llvm::ArrayRef<Algorithm> algorithms,
		Annotation annotation)
		: location(std::move(location)),
			name(std::move(name)),
			pure(pure),
			annotation(std::move(annotation)),
			type(Type::unknown())
{
	assert(!this->name.empty());

	for (const auto& member : members)
		this->members.emplace_back(std::make_shared<Member>(member));

	for (const auto& algorithm : algorithms)
		this->algorithms.emplace_back(std::make_shared<Algorithm>(algorithm));
}

Member& Function::operator[](llvm::StringRef str)
{
	for (auto& member : members)
		if (member->getName() == str)
			return *member;

	assert(false && "Not found");
}

const Member& Function::operator[](llvm::StringRef str) const
{
	for (const auto& member : members)
		if (member->getName() == str)
			return *member;

	assert(false && "Not found");
}

void Function::dump() const { dump(llvm::outs(), 0); }

void Function::dump(llvm::raw_ostream& os, size_t indents) const
{
	os.indent(indents);
	os << "function " << name << "\n";

	for (const auto& member : members)
		member->dump(os, indents + 1);

	for (const auto& algorithm : algorithms)
		algorithm->dump(os, indents + 1);
}

SourcePosition Function::getLocation() const { return location; }

std::string& Function::getName() { return name; }

const std::string& Function::getName() const { return name; }

bool Function::isPure() const
{
	return pure;
}

Container<Member>& Function::getMembers() { return members; }

const Container<Member>& Function::getMembers() const{ return members; }

Container<Member> Function::getArgs() const
{
	Container<Member> args;

	for (const auto& member : members)
		if (member->isInput())
			args.push_back(member);

	return args;
}

Container<Member> Function::getResults() const
{
	Container<Member> results;

	for (const auto& member : members)
		if (member->isOutput())
			results.push_back(member);

	return results;
}

Container<Member> Function::getProtectedMembers() const
{
	Container<Member> results;

	for (const auto& member : members)
		if (!member->isInput() && !member->isOutput())
			results.push_back(member);

	return results;
}

void Function::addMember(Member member) {
	members.emplace_back(std::make_shared<Member>(std::move(member)));
}

Container<Algorithm>& Function::getAlgorithms() { return algorithms; }

const Container<Algorithm>& Function::getAlgorithms() const
{
	return algorithms;
}

Annotation Function::getAnnotation() const
{
	return annotation;
}

Type& Function::getType() { return type; }

const Type& Function::getType() const { return type; }

void Function::setType(Type t) { type = std::move(t); }
