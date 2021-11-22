#ifndef MARCO_MATCHING_MCIM_H
#define MARCO_MATCHING_MCIM_H

#include <boost/numeric/ublas/matrix.hpp>

#include "AccessFunction.h"
#include "MCIS.h"
#include "Range.h"

namespace marco::matching
{
	using Matrix = boost::numeric::ublas::matrix<int>;

	class MCIMElement
	{
		public:
		MCIMElement(long delta, MCIS k);

		private:
		long delta;
		MCIS k;
	};

	/**
	 * Multidimensional Compressed Index Map (MCIM).
	 *
	 * It replaces the multidimensional incidence matrices in order to achieve
	 * O(1) scaling.
	 */
	class MCIM
	{
		public:
		MCIM();

		private:
	};
}

#endif	// MARCO_MATCHING_MCIM_H