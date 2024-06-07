#ifndef _CONNECTIVITY_BUILDER_H_
#define _CONNECTIVITY_BUILDER_H_

#include "CrConnectivity.h"
#include "PreBioExport.h"

#include <memory>

namespace ApBio
{
  class PREBIO_EXPORT ConnectivityBuilder
  {
  public:
    using IndexArray = CrConnectivity::IndexArray;

    ConnectivityBuilder(int const nRows, int const nCols);

    void addCount(int const i, int const count = 1);
    void finalizeCount();
    void addEntry(int const i, int const j);

    std::unique_ptr<CrConnectivity> makeConnectivity();

  private:
    int const _nRows;
    int const _nCols;
    IndexArray _offsets;
    IndexArray _entries;
  };
}

#endif