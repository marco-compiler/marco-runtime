#ifndef MARCO_RUNTIME_MODELING_MULTIDIMENSIONALRANGE_H
#define MARCO_RUNTIME_MODELING_MULTIDIMENSIONALRANGE_H

#include "marco/Runtime/Modeling/Range.h"
#include <functional>
#include <vector>

namespace marco::runtime {
using MultidimensionalRange = std::vector<Range>;

class MultidimensionalRangeIterator {
public:
  using iterator_category = std::forward_iterator_tag;
  using value_type = const int64_t *;
  using difference_type = std::ptrdiff_t;
  using pointer = const int64_t **;
  using reference = const int64_t *&;

  static MultidimensionalRangeIterator
  begin(const MultidimensionalRange &range);
  static MultidimensionalRangeIterator end(const MultidimensionalRange &range);

  bool operator==(const MultidimensionalRangeIterator &it) const;

  bool operator!=(const MultidimensionalRangeIterator &it) const;

  MultidimensionalRangeIterator &operator++();

  MultidimensionalRangeIterator operator++(int);

  const int64_t *operator*() const;

private:
  MultidimensionalRangeIterator(
      const MultidimensionalRange &range,
      std::function<RangeIterator(const Range &)> initFunction);

  void fetchNext();

private:
  std::vector<RangeIterator> beginIterators;
  std::vector<RangeIterator> currentIterators;
  std::vector<RangeIterator> endIterators;
  std::vector<int64_t> indices;
};

uint64_t getFlatSize(const MultidimensionalRange &ranges);

uint64_t getFlatIndex(const std::vector<int64_t> &indices,
                      const MultidimensionalRange &ranges);

void getIndicesFromFlatIndex(uint64_t flatIndex, std::vector<int64_t> &result,
                             const MultidimensionalRange &ranges);

void getBeginIndices(std::vector<int64_t> &indices,
                     const MultidimensionalRange &ranges);

void getEndIndices(std::vector<int64_t> &indices,
                   const MultidimensionalRange &ranges);

/// Given an array of indices and the dimensions of an equation, increase the
/// indices within the induction bounds of the equation. Return false if the
/// indices exceed the equation bounds, which means the computation has
/// finished, true otherwise.
bool advanceEquationIndices(int64_t *indices,
                            const MultidimensionalRange &ranges);

bool advanceEquationIndices(std::vector<int64_t> &indices,
                            const MultidimensionalRange &ranges);

bool advanceEquationIndicesUntil(std::vector<int64_t> &indices,
                                 const MultidimensionalRange &ranges,
                                 const std::vector<int64_t> &end);
} // namespace marco::runtime

#endif // MARCO_RUNTIME_MODELING_MULTIDIMENSIONALRANGE_H
