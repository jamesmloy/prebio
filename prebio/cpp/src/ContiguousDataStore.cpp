#include "ContiguousDataStore.h"


ContiguousDataStore::ContiguousDataStore(unsigned int const numElems)
  : _numElems(numElems)
{
}


ContiguousDataStore::~ContiguousDataStore() {}


ContiguousDataStore::ContiguousDataStore(ContiguousDataStore &&other)
  : _numElems(other.getElemCount())
  , _dataMap(std::move(other._dataMap))
{
}


ContiguousDataStore&
ContiguousDataStore::operator=(ContiguousDataStore &&other)
{
  _numElems = other.getElemCount();
  _dataMap.clear();
  _dataMap = std::move(other._dataMap);
  return *this;
}


unsigned int
ContiguousDataStore::getElemCount() const
{
  return _numElems;
}