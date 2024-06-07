#include "ConnectivityBuilder.h"

#include <stdexcept>
#include <algorithm>
#include <sstream>
#include <iostream>

namespace ApBio
{

  ConnectivityBuilder::
  ConnectivityBuilder(int const nRows, int const nCols)
    : _nRows(nRows)
    , _nCols(nCols)
    , _offsets(nRows + 1, 0)
    , _entries()
  {}


  void
  ConnectivityBuilder::
  addCount(int const i, int const count)
  {
    if (!(i < _nRows))
    {
      std::stringstream msg;
      msg << "Adding count for row " << i
          << " which is larger than num rows "
          << _nRows;
      throw std::runtime_error(msg.str());
    }
    _offsets[i] += count;
  }


  void
  ConnectivityBuilder::
  finalizeCount()
  {
    int nnz = 0;
    for (int i = 0; i < _nRows + 1; ++i)
    {
      int const cnt = _offsets[i];
      _offsets[i] = nnz;
      nnz += cnt;
    }

    _entries.resize(nnz + 1, false);
  }


  void
  ConnectivityBuilder::
  addEntry(int const i, int const j)
  {
    if (!(i < _nRows))
    {
      std::stringstream msg;
      msg << "Adding entry for row " << i
          << " which is larger than num rows "
          << _nRows;
      throw std::runtime_error(msg.str());
    }

    if (!(j < _nCols))
    {
      std::stringstream msg;
      msg << "Adding entry for column " << j
          << " which is larger than num cols "
          << _nCols;
      throw std::runtime_error(msg.str());
    }

    _entries[_offsets[i]++] = j;
  }


  std::unique_ptr<CrConnectivity>
  ConnectivityBuilder::
  makeConnectivity()
  {
    int nnz = 0;
    for (int i = 0; i < _nRows + 1; ++i)
    {
      int const nnzOld = _offsets[i];
      _offsets[i] = nnz;
      nnz = nnzOld;
    }

    return std::make_unique<CrConnectivity>(_nRows, std::move(_offsets), std::move(_entries));
  }

}