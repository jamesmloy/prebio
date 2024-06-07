#ifndef _DATA_KEY_H_
#define _DATA_KEY_H_

#include <string>

#include "PreBioExport.h"

#include "NumericTraits.h"

struct PREBIO_EXPORT DataKey
{
  std::string const _tag;
  std::string const _typeTag;
  size_t const _hash;

private:
  DataKey(std::string tag, std::string typeTag);

  template <typename T>
  friend DataKey createDataKey(std::string const& tag);
};


namespace std
{
  template <>
  struct hash<DataKey>
  {
    size_t operator()(DataKey const &k) const
    {
      return k._hash;
    }
  };
}


inline
bool operator==(DataKey const &a, DataKey const &b)
{
  return a._hash == b._hash;
}


inline
bool operator!=(DataKey const &a, DataKey const &b)
{
  return !(a == b);
}


template <typename T>
DataKey createDataKey(std::string const &tag)
{
  return {tag, NumericTraits<T>::name()};
}


#endif