#pragma once

#include <cstddef>

#include "llvm/ADT/SmallVector.h"
#include "modelica/model/ModExp.hpp"
namespace modelica
{
	class ModExpPath
	{
		public:
		ModExpPath(const ModExp& exp, llvm::SmallVector<size_t, 3> path, bool left)
				: path(std::move(path)), exp(&exp), left(left)
		{
		}

		[[nodiscard]] const ModExp& getExp() const { return *exp; }
		[[nodiscard]] bool isOnEquationLeftHand() const { return left; }
		[[nodiscard]] size_t depth() const { return path.size(); }

		[[nodiscard]] auto begin() const { return path.begin(); }
		[[nodiscard]] auto end() const { return path.end(); }

		[[nodiscard]] ModExp& reach(ModExp& exp) const
		{
			ModExp* e = &exp;
			for (auto i : path)
				e = &e->getChild(i);

			return *e;
		}
		[[nodiscard]] const ModExp& reach(const ModExp& exp) const
		{
			const ModExp* e = &exp;
			for (auto i : path)
				e = &e->getChild(i);

			return *e;
		}

		private:
		llvm::SmallVector<size_t, 3> path;
		const ModExp* exp;
		bool left;
	};
}	 // namespace modelica
