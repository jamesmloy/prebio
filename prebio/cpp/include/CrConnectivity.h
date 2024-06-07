#ifndef _CRCONNECTIVITY_H_
#define _CRCONNECTIVITY_H_

#include "PreBioExport.h"

#include "blaze/math/DynamicVector.h"

#include <utility>

namespace ApBio
{

  class PREBIO_EXPORT CrConnectivity
  {
  public:

    using IndexArray = blaze::DynamicVector<int>;

    CrConnectivity(int nRows, IndexArray offsets, IndexArray entries);

    ~CrConnectivity() {}

    CrConnectivity(CrConnectivity &&o) = default;
    CrConnectivity(CrConnectivity &o) = delete;

    CrConnectivity& operator=(CrConnectivity &&o) = default;
    CrConnectivity& operator=(CrConnectivity const &o) = delete;

    inline int const& neib(int const i, int const off) const
    {
      return _entries[_offsets[i] + off];
    }

    inline int count(int const i) const
    {
      return _offsets[i + 1] - _offsets[i];
    }

    inline int nRows() const noexcept
    {
      return _nRows;
    }

    using ConstIterator = IndexArray::ConstIterator;
    using ConstItPair = std::pair<ConstIterator, ConstIterator>;

    inline ConstItPair neibs(int const i) const
    {
      auto start = begin(_offsets) + i;
      auto finish = begin(_offsets) + i + 1;
      return std::make_pair(start, finish);
    }

  private:
    int _nRows;
    IndexArray _offsets;
    IndexArray _entries;
  };

}

#endif