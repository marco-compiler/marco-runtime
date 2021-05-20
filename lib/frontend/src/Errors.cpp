#include <modelica/frontend/Errors.h>

using namespace modelica;
using namespace modelica::frontend;

ParserErrorCategory ParserErrorCategory::category;
char UnexpectedToken::ID;
char UnexpectedIdentifier::ID;
char NotImplemented::ID;
char ChoiseNotFound::ID;
char IncompatibleType::ID;
char BranchesTypeDoNotMatch::ID;
char EmptyList::ID;
char BadSemantic::ID;

std::error_condition modelica::frontend::make_error_condition(ParserErrorCode errc)
{
	return std::error_condition(
			static_cast<int>(errc), ParserErrorCategory::category);
}

/**
 * This is required by std::error, just add a line every time you need to
 * create a new error type.
 */
[[nodiscard]] std::error_condition ParserErrorCategory::default_error_condition(
		int ev) const noexcept
{
	if (ev == 0)
		return std::error_condition(ParserErrorCode::success);
	if (ev == 1)
		return std::error_condition(ParserErrorCode::not_implemented);
	if (ev == 2)
		return std::error_condition(ParserErrorCode::unexpected_token);
	if (ev == 3)
		return std::error_condition(ParserErrorCode::choise_not_found);
	if (ev == 4)
		return std::error_condition(ParserErrorCode::branches_types_do_not_match);
	if (ev == 5)
		return std::error_condition(ParserErrorCode::incompatible_type);
	if (ev == 6)
		return std::error_condition(ParserErrorCode::empty_list);
	if (ev == 7)
		return std::error_condition(ParserErrorCode::bad_semantic);

	return std::error_condition(ParserErrorCode::unexpected_token);
}

[[nodiscard]] bool ParserErrorCategory::equivalent(
		const std::error_code& code, int condition) const noexcept
{
	bool equal = *this == code.category();
	auto v = default_error_condition(code.value()).value();
	equal = equal && static_cast<int>(v) == condition;
	return equal;
}

/**
 * Decides the messaged based upon the type.
 * This is done for compatibility with std::error, but when writing
 * tools code you should report error with ExitOnError and that will use
 * the string provided by the class extending ErrorInfo.
 */
[[nodiscard]] std::string ParserErrorCategory::message(int ev) const noexcept
{
	switch (ev)
	{
		case (0):
			return "Success";
		case (1):
			return "Not Implemented";
		case (2):
			return "Unexpected Token";
		case (3):
			return "Choise Not Found";
		case (4):
			return "If else branches type do not match";
		case (5):
			return "Expression type is incompatible";
		case (6):
			return "List was empty when expected not";
		default:
			return "Unknown Error";
	}
}

UnexpectedToken::UnexpectedToken(SourceRange location, Token token)
		: location(std::move(location)),
			token(token)
{
}

SourceRange UnexpectedToken::getLocation() const
{
	return location;
}

void UnexpectedToken::printMessage(llvm::raw_ostream& os) const
{
	os << "unexpected token [";
	os.changeColor(llvm::raw_ostream::SAVEDCOLOR, true);
	os << token;
	os << "]";
}

void UnexpectedToken::log(llvm::raw_ostream& os) const
{
	print(os);
}

UnexpectedIdentifier::UnexpectedIdentifier(SourceRange location,
																					 llvm::StringRef identifier,
																					 llvm::StringRef expected)
		: location(std::move(location)),
			identifier(identifier.str()),
			expected(expected.str())
{
}

SourceRange UnexpectedIdentifier::getLocation() const
{
	return location;
}

void UnexpectedIdentifier::printMessage(llvm::raw_ostream& os) const
{
	os << "unexpected identifier \"";
	os.changeColor(llvm::raw_ostream::SAVEDCOLOR, true);
	os << identifier;
	os.resetColor();
	os << "\" (expected: \"";
	os.changeColor(llvm::raw_ostream::SAVEDCOLOR, true);
	os << expected;
	os.resetColor();
	os << "\")";
}

void UnexpectedIdentifier::log(llvm::raw_ostream& os) const
{
	print(os);
}

BadSemantic::BadSemantic(SourceRange location, llvm::StringRef message)
		: location(std::move(location)),
			message(message.str())
{
}

SourceRange BadSemantic::getLocation() const
{
	return location;
}

void BadSemantic::printMessage(llvm::raw_ostream& os) const
{
	os << message;
}

void BadSemantic::log(llvm::raw_ostream& os) const
{
	print(os);
}
