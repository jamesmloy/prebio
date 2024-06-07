#ifndef _SINGLE_ELEMENT_CHANNEL_H_
#define _SINGLE_ELEMENT_CHANNEL_H_

#include "ChannelTag.h"
#include "PreBioExport.h"

#include "Atom.h"


template <typename T>
class PREBIO_EXPORT SingleElementChannel : public TypedChannelTag<T>
{
public:
  using Base = TypedChannelTag<T>;
  using DType = T;

  explicit SingleElementChannel(std::string elementSymbol)
    : Base("SingleElement-" + elementSymbol)
    , _elementSymbol(std::move(elementSymbol))
  {
  }

  std::vector<std::string> getLabelTags() const
  {
    return {_elementSymbol};
  }

  bool hasData(MolecularObjs::Atom const &a) const
  {
    return DType(a.element() == _elementSymbol);
  }

  DType extract(MolecularObjs::Atom const &a) const
  {
    return DType(a.element() == _elementSymbol);
  }

private:
  std::string const _elementSymbol;
};


#endif