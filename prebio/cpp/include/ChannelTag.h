#ifndef _CHANNEL_TAG_H_
#define _CHANNEL_TAG_H_

#include "DataKey.h"
#include "PreBioExport.h"

class PREBIO_EXPORT ChannelTag
{
public:

  virtual ~ChannelTag();

  inline std::string const &channelName() const;

  inline DataKey const &getKey() const;

protected:

  /// cannot instantiate this directly, it
  /// is only created by derived classes
  inline ChannelTag(DataKey key);

private:
  DataKey _key;
};


/// pass through class for convenience
template <typename T>
class TypedChannelTag : public ChannelTag
{
public:
  using DType = T;
protected:
  TypedChannelTag(std::string const &name);
};


ChannelTag::ChannelTag(DataKey key)
  : _key(key)
{
}


std::string const &
ChannelTag::channelName() const
{
  return _key._tag;
}


DataKey const &
ChannelTag::getKey() const
{
  return _key;
}


template <typename T>
TypedChannelTag<T>::TypedChannelTag(std::string const &name)
  : ChannelTag(createDataKey<T>(name))
{}

#endif