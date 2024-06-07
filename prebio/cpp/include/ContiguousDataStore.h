#ifndef _CONTIGUOUS_DATA_STORE_H_
#define _CONTIGUOUS_DATA_STORE_H_

#include <unordered_map>
#include <memory>
#include <iostream>
#include <sstream>

#include "DataKey.h"
#include "PreBioExport.h"

#include "blaze/math/DynamicVector.h"


class PREBIO_EXPORT ContiguousDataStore
{
  struct DataHandle
  {
    virtual ~DataHandle() {}
  };

  template <typename NumType>
  struct SpecificDataHandle : DataHandle
  {
    using Data = blaze::DynamicVector<NumType>;
    Data _data;

    SpecificDataHandle(Data &&data)
      : _data(data)
    {
    }

    Data &getData() noexcept { return _data; }
    Data const &getData() const noexcept { return _data; }
  };

  using DhPtr = std::unique_ptr<DataHandle>;
  using DataMap = std::unordered_map<DataKey, DhPtr>;

public:
  template <typename T>
  using VecType = blaze::DynamicVector<T>;

  ContiguousDataStore(unsigned int const numElems);

  ~ContiguousDataStore();

  /// cannot copy construct
  ContiguousDataStore(ContiguousDataStore const &other) = delete;
  /// cannot copy assign
  ContiguousDataStore& operator=(ContiguousDataStore const &other) = delete;

  /// move ctor
  ContiguousDataStore(ContiguousDataStore &&other);
  /// move assign
  ContiguousDataStore& operator=(ContiguousDataStore &&other);

  unsigned int getElemCount() const;

  template <typename T>
  VecType<T> &getData(std::string const &tag)
  {
    auto const key = createDataKey<T>(tag);
    auto it = _dataMap.find(key);

    if (it == end(_dataMap))
    {
      std::stringstream msg;
      msg << "Attempting to access unallocated storage: ";
      msg << key._tag << ", ";
      msg << key._typeTag << std::endl;
      throw std::runtime_error(msg.str());
    }

    auto &dh = *it->second;
    auto &sdh = static_cast<SpecificDataHandle<T> &>(dh);
    return sdh.getData();
  }

  template <typename T>
  VecType<T> const &getData(std::string const &tag) const
  {
    auto const key = createDataKey<T>(tag);
    auto it = _dataMap.find(key);

    if (it == end(_dataMap))
    {
      std::stringstream msg;
      msg << "Attempting to access unallocated storage: ";
      msg << key._tag << ", ";
      msg << key._typeTag << std::endl;
      throw std::runtime_error(msg.str());
    }

    auto &dh = *it->second;
    auto const &sdh = static_cast<SpecificDataHandle<T> const &>(dh);
    return sdh.getData();
  }

  template <typename T>
  bool hasData(std::string const &tag) const
  {
    auto const key = createDataKey<T>(tag);
    auto it = _dataMap.find(key);

    return it != end(_dataMap);
  }

  template <typename T>
  void addData(std::string const &tag)
  {
    VecType<T> data(_numElems, T(0));

    auto key = createDataKey<T>(tag);
    auto sdh = std::make_unique<SpecificDataHandle<T>>(std::move(data));
    _dataMap.emplace(key, std::move(sdh));
  }

  inline void clearData() noexcept
  {
    _dataMap.clear();
  }

private:
  unsigned int _numElems;
  DataMap _dataMap;
};

#endif