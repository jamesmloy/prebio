#ifndef _MARK_UTILS_H_
#define _MARK_UTILS_H_

#include "gemmi/neighbor.hpp"

#include <iostream>
#include <utility>
#include <vector>

namespace ApBio
{
  namespace __MarkUtilsImpl
  {
    using Mark = gemmi::NeighborSearch::Mark;

    using MarkPtrVec = std::vector<Mark *>;
    using MarkPtrVecIt = MarkPtrVec::iterator;
    using MarkPtrVecConstIt = MarkPtrVec::const_iterator;

    using ConstMarkPtrVec = const MarkPtrVec;
    using ConstMarkPtrVecIt = ConstMarkPtrVec::iterator;
    using ConstMarkPtrVecConstIt = ConstMarkPtrVec::const_iterator;

    using MarkConstPtrVec = std::vector<Mark const *>;
    using MarkConstPtrVecIt = MarkConstPtrVec::iterator;
    using MarkConstPtrVecConstIt = MarkConstPtrVec::const_iterator;

    using ConstMarkConstPtrVec = const MarkConstPtrVec;
    using ConstMarkConstPtrVecIt = ConstMarkConstPtrVec::iterator;
    using ConstMarkConstPtrVecConstIt = ConstMarkConstPtrVec::const_iterator;

    using ConstMarkConstPtrVecRng
      = std::pair<ConstMarkConstPtrVecConstIt, ConstMarkConstPtrVecConstIt>;

    using ConstMarkConstPtrVecRngList = std::vector<ConstMarkConstPtrVecRng>;

    using ConstMarkPtrRng
      = std::pair<Mark const *const *const, Mark const *const *const>;
    using ConstMarkPtrRngList = std::vector<ConstMarkPtrRng>;

    template <typename PropFunc>
    __MarkUtilsImpl::ConstMarkPtrRngList __groupBy(
      __MarkUtilsImpl::Mark const *const *const beg,
      __MarkUtilsImpl::Mark const *const *const fin, PropFunc &&func)
    {
      using namespace __MarkUtilsImpl;

      // empty list if empty
      if (beg == fin)
        return {};

      auto begRng = beg;
      auto endRng = beg;

      ConstMarkPtrRngList theRanges;

      while (begRng != fin)
      {
        if (endRng != fin && func(**begRng) == func(**endRng))
          ++endRng;
        else
        {
          theRanges.emplace_back(begRng, endRng);
          begRng = endRng;
        }
      }

      return theRanges;
    }
  } // namespace __MarkUtilsImpl

  struct MarkSorterOld
  {
    inline bool operator()(gemmi::NeighborSearch::Mark const *const a,
      gemmi::NeighborSearch::Mark const *const b) const
    {
      return this->operator()(*a, *b);
    }

    inline bool operator()(
      gemmi::NeighborSearch::Mark const &a, gemmi::NeighborSearch::Mark const &b) const
    {
      if (a.chain_idx == b.chain_idx)
      {
        if (a.residue_idx == b.residue_idx)
        {
          return a.atom_idx < b.atom_idx;
        }
        else
          return a.residue_idx < b.residue_idx;
      }
      else
        return a.chain_idx < b.chain_idx;
    }
  };


  template <typename PropFunc>
  __MarkUtilsImpl::ConstMarkPtrRngList groupBy(
    std::vector<__MarkUtilsImpl::Mark const *> const &vec, PropFunc &&func)
  {
    using namespace __MarkUtilsImpl;

    auto const beg = &vec.front();
    auto const fin = beg + vec.size();
    return __groupBy(beg, fin, func);
  }


  template <typename PropFunc>
  __MarkUtilsImpl::ConstMarkPtrRngList groupBy(
    std::vector<__MarkUtilsImpl::Mark *> const &vec, PropFunc &&func)
  {
    using namespace __MarkUtilsImpl;

    auto const beg = &vec.front();
    auto const fin = beg + vec.size();
    return __groupBy(beg, fin, func);
  }
} // namespace ApBio

#endif