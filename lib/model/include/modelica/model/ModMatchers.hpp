#pragma once
#include "llvm/ADT/SmallPtrSet.h"
#include "llvm/ADT/SmallVector.h"
#include "modelica/model/ModEquation.hpp"
#include "modelica/model/ModExp.hpp"

namespace modelica
{
	class ReferenceMatcher
	{
		public:
		ReferenceMatcher() = default;
		ReferenceMatcher(const ModEquation& eq) { visit(eq); }
		void visitAt(const ModExp& exp)
		{
			const ModExp* childs = &exp.getLeftHand();
			toIgnore.insert(childs);

			while (childs->isOperation<ModExpKind::at>())
			{
				childs = &childs->getLeftHand();
				toIgnore.insert(childs);
			}

			if (childs->isReference())
				vars.emplace_back(&exp);
		}

		void afterVisit(const ModExp&) {}

		void visit(const ModExp& exp)
		{
			if (toIgnore.find(&exp) != toIgnore.end())
				return;
			if (exp.isOperation<ModExpKind::at>())
			{
				visitAt(exp);
				return;
			}

			if (exp.isReference())
			{
				toIgnore.insert(&exp);
				vars.emplace_back(&exp);
			}
		}

		void visit(const ModEquation& equation)
		{
			modelica::visit(equation.getLeft(), *this);
			modelica::visit(equation.getRight(), *this);
		}

		[[nodiscard]] auto begin() const { return vars.begin(); }
		[[nodiscard]] auto end() const { return vars.end(); }
		[[nodiscard]] auto begin() { return vars.begin(); }
		[[nodiscard]] auto end() { return vars.end(); }
		[[nodiscard]] size_t size() const { return vars.size(); }
		[[nodiscard]] const ModExp& at(size_t index) const { return *vars[index]; }
		[[nodiscard]] const ModExp& operator[](size_t index) const
		{
			return at(index);
		}

		private:
		llvm::SmallPtrSet<const ModExp*, 4> toIgnore;
		llvm::SmallVector<const ModExp*, 3> vars;
	};

};	// namespace modelica
