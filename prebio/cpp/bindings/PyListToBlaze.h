#ifndef _PYLIST_TO_BLAZE_H_
#define _PYLIST_TO_BLAZE_H_

#include "pybind11/pybind11.h"
#include "blaze/math/StaticVector.h"

#include <sstream>

template <typename T, int N>
blaze::StaticVector<T, N>
copyPyListToStaticVector(pybind11::list theList)
{
  if (theList.size() != N)
  {
    std::stringstream msg;
    msg << "Trying to convert a list to a static vector length ";
    msg << N << ", ";
    msg << "but the list size is " << theList.size();
    throw std::runtime_error(msg.str());
  }
  blaze::StaticVector<T, N> theVector;
  for (int i = 0; i != N; ++i)
  {
    theVector[i] = theList[i].cast<T>();
  }
  // for (auto listItem : theList) { theVector.push_back(listItem.cast<T>()); }
  return theVector;
}

#endif