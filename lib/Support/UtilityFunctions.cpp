#include "marco/Runtime/Support/UtilityFunctions.h"
#include <cstring>
#include <iostream>
#include <numeric>

//===----------------------------------------------------------------------===//
// clone
//===----------------------------------------------------------------------===//

namespace {
/// Clone an array into another one.
///
/// @tparam T 					destination array type
/// @tparam U 					source array type
/// @param destination  destination array
/// @param values 			source values
template <typename T, typename U>
void clone_void(UnrankedMemRefType<T> *destination,
                UnrankedMemRefType<U> *source) {
  DynamicMemRefType<U> dynamicSource(*source);
  DynamicMemRefType<T> dynamicDestination(*destination);

  // Check that the two arrays have the same number of elements
  [[maybe_unused]] int64_t sourceFlatSize = std::accumulate(
      dynamicSource.sizes, dynamicSource.sizes + dynamicSource.rank,
      static_cast<int64_t>(1), std::multiplies<int64_t>());

  [[maybe_unused]] int64_t destinationFlatSize =
      std::accumulate(dynamicDestination.sizes,
                      dynamicDestination.sizes + dynamicDestination.rank,
                      static_cast<int64_t>(1), std::multiplies<int64_t>());

  assert(sourceFlatSize == destinationFlatSize);

  auto sourceIt = std::begin(dynamicSource);
  auto destinationIt = std::begin(dynamicDestination);

  for (size_t i = 0, e = sourceFlatSize; i < e; ++i) {
    *destinationIt = *sourceIt;

    ++sourceIt;
    ++destinationIt;
  }
}
} // namespace

RUNTIME_FUNC_DEF(clone, void, ARRAY(bool), ARRAY(bool))
RUNTIME_FUNC_DEF(clone, void, ARRAY(bool), ARRAY(int32_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(bool), ARRAY(int64_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(bool), ARRAY(float))
RUNTIME_FUNC_DEF(clone, void, ARRAY(bool), ARRAY(double))

RUNTIME_FUNC_DEF(clone, void, ARRAY(int32_t), ARRAY(bool))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int32_t), ARRAY(int32_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int32_t), ARRAY(int64_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int32_t), ARRAY(float))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int32_t), ARRAY(double))

RUNTIME_FUNC_DEF(clone, void, ARRAY(int64_t), ARRAY(bool))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int64_t), ARRAY(int32_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int64_t), ARRAY(int64_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int64_t), ARRAY(float))
RUNTIME_FUNC_DEF(clone, void, ARRAY(int64_t), ARRAY(double))

RUNTIME_FUNC_DEF(clone, void, ARRAY(float), ARRAY(bool))
RUNTIME_FUNC_DEF(clone, void, ARRAY(float), ARRAY(int32_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(float), ARRAY(int64_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(float), ARRAY(float))
RUNTIME_FUNC_DEF(clone, void, ARRAY(float), ARRAY(double))

RUNTIME_FUNC_DEF(clone, void, ARRAY(double), ARRAY(bool))
RUNTIME_FUNC_DEF(clone, void, ARRAY(double), ARRAY(int32_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(double), ARRAY(int64_t))
RUNTIME_FUNC_DEF(clone, void, ARRAY(double), ARRAY(float))
RUNTIME_FUNC_DEF(clone, void, ARRAY(double), ARRAY(double))

extern "C" void memrefCopy(int64_t elemSize, UnrankedMemRefType<char> *srcArg,
                           UnrankedMemRefType<char> *dstArg) {
  DynamicMemRefType<char> src(*srcArg);
  DynamicMemRefType<char> dst(*dstArg);

  int64_t rank = src.rank;

  // Handle empty shapes -> nothing to copy.
  for (int rankp = 0; rankp < rank; ++rankp)
    if (src.sizes[rankp] == 0)
      return;

  char *srcPtr = src.data + src.offset * elemSize;
  char *dstPtr = dst.data + dst.offset * elemSize;

  if (rank == 0) {
    memcpy(dstPtr, srcPtr, elemSize);
    return;
  }

  int64_t *indices = static_cast<int64_t *>(alloca(sizeof(int64_t) * rank));
  int64_t *srcStrides = static_cast<int64_t *>(alloca(sizeof(int64_t) * rank));
  int64_t *dstStrides = static_cast<int64_t *>(alloca(sizeof(int64_t) * rank));

  // Initialize index and scale strides.
  for (int rankp = 0; rankp < rank; ++rankp) {
    indices[rankp] = 0;
    srcStrides[rankp] = src.strides[rankp] * elemSize;
    dstStrides[rankp] = dst.strides[rankp] * elemSize;
  }

  int64_t readIndex = 0, writeIndex = 0;
  for (;;) {
    // Copy over the element, byte by byte.
    memcpy(dstPtr + writeIndex, srcPtr + readIndex, elemSize);
    // Advance index and read position.
    for (int64_t axis = rank - 1; axis >= 0; --axis) {
      // Advance at current axis.
      auto newIndex = ++indices[axis];
      readIndex += srcStrides[axis];
      writeIndex += dstStrides[axis];
      // If this is a valid index, we have our next index, so continue copying.
      if (src.sizes[axis] != newIndex)
        break;
      // We reached the end of this axis. If this is axis 0, we are done.
      if (axis == 0)
        return;
      // Else, reset to 0 and undo the advancement of the linear index that
      // this axis had. Then continue with the axis one outer.
      indices[axis] = 0;
      readIndex -= src.sizes[axis] * srcStrides[axis];
      writeIndex -= dst.sizes[axis] * dstStrides[axis];
    }
  }
}
