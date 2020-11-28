#pragma once

#include <llvm/ADT/ArrayRef.h>
#include <llvm/ADT/SmallVector.h>
#include <llvm/Support/raw_ostream.h>
#include <modelica/frontend/Algorithm.hpp>
#include <modelica/frontend/Class.hpp>
#include <modelica/frontend/Equation.hpp>
#include <modelica/frontend/ForEquation.hpp>
#include <modelica/frontend/Function.hpp>
#include <modelica/frontend/Member.hpp>
#include <string>

namespace modelica
{
	enum class ClassType
	{
		Function,
		Model
	};

	class ClassContainer
	{
		public:
		explicit ClassContainer(Function function);
		explicit ClassContainer(Class model);

		void dump() const;
		void dump(llvm::raw_ostream& os, size_t indents = 0) const;

		template<typename T>
		[[nodiscard]] bool isA() const
		{
			return std::holds_alternative<T>(content);
		}

		template<typename T>
		[[nodiscard]] const T& get() const
		{
			assert(isA<T>());
			return std::get<T>(content);
		}

		template<class Visitor>
		auto visit(Visitor&& vis)
		{
			return std::visit(std::forward<Visitor>(vis), content);
		}

		template<class Visitor>
		auto visit(Visitor&& vis) const
		{
			return std::visit(std::forward<Visitor>(vis), content);
		}

		private:
		std::variant<Function, Class> content;
	};
}	 // namespace modelica
