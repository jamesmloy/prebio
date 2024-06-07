  
#include "CrConnectivity.h"

namespace ApBio
{

  CrConnectivity::
  CrConnectivity(int nRows, CrConnectivity::IndexArray offsets, CrConnectivity::IndexArray entries)
    : _nRows(nRows)
    , _offsets(std::move(offsets))
    , _entries(std::move(entries))
  {
  }

}