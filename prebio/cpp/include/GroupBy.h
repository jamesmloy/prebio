#ifndef _GROUP_BY_H_
#define _GROUP_BY_H_

#include <utility>
#include <vector>

namespace __groupByImpl
{
  struct DefaultEqFunc
  {
    template <typename T>
    auto operator()(T const &a, T const &b) const
    {
      return a == b;
    }
  };
} // namespace __groupByImpl


template <typename It, typename PropFunc = __groupByImpl::DefaultEqFunc>
auto
groupByFromSortedRange(It beg, It fin,
  PropFunc &&func = __groupByImpl::DefaultEqFunc()) -> decltype(auto)
{
  using RangePair = std::pair<It, It>;
  using RangeList = std::vector<RangePair>;

  // empty list if empty
  if (beg == fin)
    return RangeList{};

  auto begRng = beg;
  auto endRng = beg;

  RangeList theRanges;

  while (begRng != fin)
  {
    if (endRng != fin && func(*begRng, *endRng))
      ++endRng;
    else
    {
      theRanges.emplace_back(begRng, endRng);
      begRng = endRng;
    }
  }

  return theRanges;
}


template <typename ListType, typename PropFunc = __groupByImpl::DefaultEqFunc>
auto
groupByFromSorted(ListType const &theList,
  PropFunc &&func = __groupByImpl::DefaultEqFunc()) -> decltype(auto)
{
  return groupByFromSortedRange(
    std::begin(theList), std::end(theList), std::forward<PropFunc &&>(func));
}

#endif