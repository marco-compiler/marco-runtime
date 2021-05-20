#pragma once

#include <llvm/ADT/SmallVector.h>
#include <llvm/Support/Error.h>
#include <memory>
#include <modelica/frontend/Errors.h>
#include <modelica/frontend/LexerStateMachine.h>
#include <modelica/utils/Lexer.hpp>
#include <optional>

#include "AST.h"

namespace modelica::frontend
{
	class Argument;
	class ClassModification;
	class ElementModification;
	class ElementRedeclaration;
	class ElementReplaceable;
	class Expression;
	class Statement;
	class AssignmentStatement;
	class IfStatement;
	class ForStatement;
	class Modification;
	class WhileStatement;
	class WhenStatement;

	/**
	 * The parser encapsulates the lexer but not the memory where the string we
	 * are reading is held. It exposes parts of the grammatical rules that are
	 * available in the grammar (can be found at page ~ 265 of the 3.4 doc).
	 */
	class Parser
	{
		public:

		template<typename T>
		class ValueWrapper
		{
			public:
			ValueWrapper(SourceRange location, T value)
					: location(std::move(location)), value(value)
			{
			}

			SourceRange getLocation()
			{
				return location;
			}

			T& getValue()
			{
				return value;
			}

			const T& getValue() const
			{
				return value;
			}

			private:
			SourceRange location;
			T value;
		};

		Parser(llvm::StringRef fileName, const char* source);
		Parser(const std::string& source);
		Parser(const char* source);

		/**
		 * Return the current position in the source stream.
		 */
		[[nodiscard]] SourceRange getPosition() const;

		[[nodiscard]] Token getCurrentToken() const;

		llvm::Expected<ValueWrapper<std::string>> identifier();

		llvm::Expected<std::unique_ptr<Class>> classDefinition();

		llvm::Expected<std::unique_ptr<Expression>> primary();
		llvm::Expected<std::unique_ptr<Expression>> factor();
		llvm::Expected<llvm::Optional<std::unique_ptr<Expression>>> termModification();
		llvm::Expected<std::unique_ptr<Expression>> term();
		llvm::Expected<Type> typeSpecifier();

		llvm::Expected<std::unique_ptr<Expression>> arithmeticExpression();

		llvm::Error arrayDimensions(
				llvm::SmallVectorImpl<ArrayDimension>& dimensions);

		llvm::Error elementList(
				llvm::SmallVectorImpl<std::unique_ptr<Member>>& members,
				bool publicSection = true);

		llvm::Expected<TypePrefix> typePrefix();
		llvm::Expected<std::unique_ptr<Member>> element(bool publicSection = true);
		llvm::Optional<OperationKind> relationalOperator();

		llvm::Expected<std::unique_ptr<Expression>> logicalTerm();
		llvm::Expected<std::unique_ptr<Expression>> logicalExpression();

		llvm::Expected<std::unique_ptr<Induction>> induction();

		llvm::Error equationSection(
				llvm::SmallVectorImpl<std::unique_ptr<Equation>>& equations,
				llvm::SmallVectorImpl<std::unique_ptr<ForEquation>>& forEquations);

		llvm::Expected<std::unique_ptr<Equation>> equation();

		llvm::Error forEquation(
				llvm::SmallVectorImpl<std::unique_ptr<ForEquation>>& equations,
				int nestingLevel);

		llvm::Error forEquationBody(
				llvm::SmallVectorImpl<std::unique_ptr<ForEquation>>& equations,
				int nestingLevel);

		llvm::Expected<std::unique_ptr<Expression>> expression();
		llvm::Expected<std::unique_ptr<Expression>> logicalFactor();
		llvm::Expected<std::unique_ptr<Expression>> relation();
		llvm::Expected<std::unique_ptr<Expression>> componentReference();

		llvm::Error functionCallArguments(llvm::SmallVectorImpl<std::unique_ptr<Expression>>& args);

		llvm::Expected<std::unique_ptr<Algorithm>> algorithmSection();
		llvm::Expected<std::unique_ptr<Statement>> statement();
		llvm::Expected<std::unique_ptr<Statement>> assignmentStatement();
		llvm::Expected<std::unique_ptr<Statement>> ifStatement();
		llvm::Expected<std::unique_ptr<Statement>> forStatement();
		llvm::Expected<std::unique_ptr<Statement>> whileStatement();
		llvm::Expected<std::unique_ptr<Statement>> whenStatement();
		llvm::Expected<std::unique_ptr<Statement>> breakStatement();
		llvm::Expected<std::unique_ptr<Statement>> returnStatement();

		llvm::Error outputExpressionList(llvm::SmallVectorImpl<std::unique_ptr<Expression>>& expressions);

		llvm::Error arraySubscript(llvm::SmallVectorImpl<std::unique_ptr<Expression>>& subscripts);

		llvm::Expected<std::unique_ptr<Annotation>> annotation();
		llvm::Expected<std::unique_ptr<Modification>> modification();
		llvm::Expected<std::unique_ptr<ClassModification>> classModification();
		llvm::Expected<std::unique_ptr<Argument>> argument();
		llvm::Expected<std::unique_ptr<Argument>> elementModification(bool each, bool final);
		llvm::Expected<std::unique_ptr<Argument>> elementRedeclaration();
		llvm::Expected<std::unique_ptr<Argument>> elementReplaceable(bool each, bool final);

		private:
		/**
		 * Read the next token.
		 */
		void next();

		/**
		 * Regular accept: if the current token is t then the next one will be read
		 * and true will be returned, else false.
		 */
		bool accept(Token t);

		/**
		 * fancy overloads if you know at compile time
		 * which token you want.
		 */
		template<Token t>
		bool accept()
		{
			if (current == t)
			{
				next();
				return true;
			}

			return false;
		}

		llvm::Expected<bool> expect(Token t);

		void updateTokenSourceRange();

		const std::string filename;
		Lexer<ModelicaStateMachine> lexer;
		Token current;
		SourceRange tokenRange;
	};
}
