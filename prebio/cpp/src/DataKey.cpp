#include "DataKey.h"

#include <iostream>


DataKey::DataKey(std::string tag, std::string typeTag)
  : _tag(std::move(tag))
  , _typeTag(std::move(typeTag))
  , _hash(std::hash<std::string>()(_tag + _typeTag))
{
}