#include <modelica/runtime/Runtime.h>
#include <numeric>

/**
 * Set all the elements of an array to a given value.
 *
 * @tparam T 		 data type
 * @param array  array to be populated
 * @param value  value to be set
 */
template<typename T>
inline void fill(UnsizedArrayDescriptor<T> array, T value)
{
	for (auto& element : array)
		element = value;
}

RUNTIME_FUNC_DEF(fill, void, ARRAY(bool), bool)
RUNTIME_FUNC_DEF(fill, void, ARRAY(int), int)
RUNTIME_FUNC_DEF(fill, void, ARRAY(long), long)
RUNTIME_FUNC_DEF(fill, void, ARRAY(float), float)
RUNTIME_FUNC_DEF(fill, void, ARRAY(double), double)

/**
 * Set a multi-dimensional array to an identity like matrix.
 *
 * @tparam T 	   data type
 * @param array  array to be populated
 */
template<typename T>
inline void identity(UnsizedArrayDescriptor<T> array)
{
	// Check that the array is square-like (all the dimensions have the same
	// size). Note that the implementation is generalized to n-D dimensions,
	// while the "identity" Modelica function is defined only for 2-D arrays.
	// Still, the implementation complexity would be the same.

	assert(array.hasSameSizes());

	// Directly use the iterators, as we need to determine the current indexes
	// so that we can place a 1 if the access is on the matrix diagonal.

	for (auto it = array.begin(), end = array.end(); it != end; ++it)
	{
		auto indexes = it.getCurrentIndexes();
		assert(!indexes.empty());

		bool isIdentityAccess = llvm::all_of(indexes, [&indexes](const auto& i) {
			return i == indexes[0];
		});

		*it = isIdentityAccess ? 1 : 0;
	}
}

RUNTIME_FUNC_DEF(identity, void, ARRAY(bool))
RUNTIME_FUNC_DEF(identity, void, ARRAY(int))
RUNTIME_FUNC_DEF(identity, void, ARRAY(long))
RUNTIME_FUNC_DEF(identity, void, ARRAY(float))
RUNTIME_FUNC_DEF(identity, void, ARRAY(double))

/**
 * Place some values on the diagonal of a matrix, and set all the other
 * elements to zero.
 *
 * @tparam T 					destination matrix type
 * @tparam U 					source values type
 * @param destination destination matrix
 * @param values 			source values
 */
template<typename T, typename U>
inline void diagonal(UnsizedArrayDescriptor<T> destination, UnsizedArrayDescriptor<U> values)
{
	// Check that the array is square-like (all the dimensions have the same
	// size). Note that the implementation is generalized to n-D dimensions,
	// while the "identity" Modelica function is defined only for 2-D arrays.
	// Still, the implementation complexity would be the same.

	assert(destination.hasSameSizes());

	// Check that the sizes of the matrix dimensions match with the amount of
	// values to be set.

	assert(destination.getRank() > 0);
	assert(values.getRank() == 1);
	assert(destination.getDimensionSize(0) == values.getDimensionSize(0));

	// Directly use the iterators, as we need to determine the current indexes
	// so that we can place a 1 if the access is on the matrix diagonal.

	for (auto it = destination.begin(), end = destination.end(); it != end; ++it)
	{
		auto indexes = it.getCurrentIndexes();
		assert(!indexes.empty());

		bool isIdentityAccess = llvm::all_of(indexes, [&indexes](const auto& i) {
			return i == indexes[0];
		});

		*it = isIdentityAccess ? values.get(indexes[0]) : 0;
	}
}

RUNTIME_FUNC_DEF(diagonal, void, ARRAY(bool), ARRAY(bool))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(bool), ARRAY(int))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(bool), ARRAY(long))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(bool), ARRAY(float))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(bool), ARRAY(double))

RUNTIME_FUNC_DEF(diagonal, void, ARRAY(int), ARRAY(bool))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(int), ARRAY(int))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(int), ARRAY(long))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(int), ARRAY(float))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(int), ARRAY(double))

RUNTIME_FUNC_DEF(diagonal, void, ARRAY(long), ARRAY(bool))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(long), ARRAY(int))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(long), ARRAY(long))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(long), ARRAY(float))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(long), ARRAY(double))

RUNTIME_FUNC_DEF(diagonal, void, ARRAY(float), ARRAY(bool))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(float), ARRAY(int))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(float), ARRAY(long))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(float), ARRAY(float))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(float), ARRAY(double))

RUNTIME_FUNC_DEF(diagonal, void, ARRAY(double), ARRAY(bool))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(double), ARRAY(int))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(double), ARRAY(long))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(double), ARRAY(float))
RUNTIME_FUNC_DEF(diagonal, void, ARRAY(double), ARRAY(double))

/**
 * Set all the elements of an array to zero.
 *
 * @tparam T data type
 * @param array   array to be populated
 */
template<typename T>
inline void zeros(UnsizedArrayDescriptor<T> array)
{
	for (auto& element : array)
		element = 0;
}

RUNTIME_FUNC_DEF(zeros, void, ARRAY(bool))
RUNTIME_FUNC_DEF(zeros, void, ARRAY(int))
RUNTIME_FUNC_DEF(zeros, void, ARRAY(long))
RUNTIME_FUNC_DEF(zeros, void, ARRAY(float))
RUNTIME_FUNC_DEF(zeros, void, ARRAY(double))

/**
 * Set all the elements of an array to ones.
 *
 * @tparam T data type
 * @param array   array to be populated
 */
template<typename T>
inline void ones(UnsizedArrayDescriptor<T> array)
{
	for (auto& element : array)
		element = 1;
}

RUNTIME_FUNC_DEF(ones, void, ARRAY(bool))
RUNTIME_FUNC_DEF(ones, void, ARRAY(int))
RUNTIME_FUNC_DEF(ones, void, ARRAY(long))
RUNTIME_FUNC_DEF(ones, void, ARRAY(float))
RUNTIME_FUNC_DEF(ones, void, ARRAY(double))

/**
 * Populate a 1-D array with equally spaced elements.
 *
 * @tparam T 		 data type
 * @param array  array to be populated
 * @param start  start value
 * @param end 	 end value
 */
template<typename T>
inline void linspace(UnsizedArrayDescriptor<T> array, double start, double end)
{
	assert(array.getRank() == 1);

	size_t n = array.getDimensionSize(0);
	double step = (end - start) / ((double) n - 1);

	for (size_t i = 0; i < n; ++i)
		array.get(i) = start + i * step;
}

RUNTIME_FUNC_DEF(linspace, void, ARRAY(bool), float, float)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(bool), double, double)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(int), float, float)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(int), double, double)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(long), float, float)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(long), double, double)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(float), float, float)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(float), double, double)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(double), float, float)
RUNTIME_FUNC_DEF(linspace, void, ARRAY(double), double, double)

template<typename T>
inline T min(UnsizedArrayDescriptor<T> array)
{
	return *std::min_element(array.begin(), array.end());
}

RUNTIME_FUNC_DEF(min, bool, ARRAY(bool))
RUNTIME_FUNC_DEF(min, int, ARRAY(int))
RUNTIME_FUNC_DEF(min, long, ARRAY(long))
RUNTIME_FUNC_DEF(min, float, ARRAY(float))
RUNTIME_FUNC_DEF(min, double, ARRAY(double))

template<typename T>
inline T min(T x, T y)
{
	return std::min(x, y);
}

RUNTIME_FUNC_DEF(min, bool, bool, bool)
RUNTIME_FUNC_DEF(min, int, int, int)
RUNTIME_FUNC_DEF(min, long, long, long)
RUNTIME_FUNC_DEF(min, float, float, float)
RUNTIME_FUNC_DEF(min, double, double, double)

template<typename T>
inline T max(UnsizedArrayDescriptor<T> array)
{
	return *std::max_element(array.begin(), array.end());
}

RUNTIME_FUNC_DEF(max, bool, ARRAY(bool))
RUNTIME_FUNC_DEF(max, int, ARRAY(int))
RUNTIME_FUNC_DEF(max, long, ARRAY(long))
RUNTIME_FUNC_DEF(max, float, ARRAY(float))
RUNTIME_FUNC_DEF(max, double, ARRAY(double))

template<typename T>
inline T max(T x, T y)
{
	return std::max(x, y);
}

RUNTIME_FUNC_DEF(max, bool, bool, bool)
RUNTIME_FUNC_DEF(max, int, int, int)
RUNTIME_FUNC_DEF(max, long, long, long)
RUNTIME_FUNC_DEF(max, float, float, float)
RUNTIME_FUNC_DEF(max, double, double, double)

template<typename T>
inline T sum(UnsizedArrayDescriptor<T> array)
{
	return std::accumulate(array.begin(), array.end(), 0, std::plus<T>());
}

RUNTIME_FUNC_DEF(sum, bool, ARRAY(bool))
RUNTIME_FUNC_DEF(sum, int, ARRAY(int))
RUNTIME_FUNC_DEF(sum, long, ARRAY(long))
RUNTIME_FUNC_DEF(sum, float, ARRAY(float))
RUNTIME_FUNC_DEF(sum, double, ARRAY(double))

template<typename T>
inline T product(UnsizedArrayDescriptor<T> array)
{
	return std::accumulate(array.begin(), array.end(), 1, std::multiplies<T>());
}

RUNTIME_FUNC_DEF(product, bool, ARRAY(bool))
RUNTIME_FUNC_DEF(product, int, ARRAY(int))
RUNTIME_FUNC_DEF(product, long, ARRAY(long))
RUNTIME_FUNC_DEF(product, float, ARRAY(float))
RUNTIME_FUNC_DEF(product, double, ARRAY(double))