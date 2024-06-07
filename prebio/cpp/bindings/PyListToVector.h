#ifndef _PYLIST_TO_VECTOR_H_
#define _PYLIST_TO_VECTOR_H_

#include "pybind11/pybind11.h"

template <typename T>
std::vector<T>
copyPyListToVector(pybind11::list theList)
{
  std::vector<T> theVector;
  theVector.reserve(theList.size());

  for (auto listItem : theList) { theVector.push_back(listItem.cast<T>()); }
  return theVector;
}

#endif