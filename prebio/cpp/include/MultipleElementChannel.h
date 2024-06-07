#ifndef _MULTIPLE_ELEMENT_CHANNEL_H_
#define _MULTIPLE_ELEMENT_CHANNEL_H_

#include "ChannelTag.h"
#include "PreBioExport.h"

#include "AtomicCollection.h"

#include <algorithm>
#include <sstream>
#include <string>


template <typename T>
class PREBIO_EXPORT MultipleElementChannel : public TypedChannelTag<T>
{

  std::string makeStorageTag(std::vector<std::string> const &elements) const
  {
    std::stringstream tag;
    tag << "MultipleElement-";
    std::for_each(
      begin(elements), end(elements), [&tag](auto const &e) { tag << e; });

    return tag.str();
  }

public:
  using Base = TypedChannelTag<T>;
  using DType = T;

  explicit MultipleElementChannel(std::vector<std::string> elements)
    : Base(makeStorageTag(elements))
    , _elements(std::move(elements))
  {
  }

  std::vector<std::string> getLabelTags() const { return _elements; }

  bool hasData(MolecularObjs::Atom const &a) const { return extract(a); }

  DType extract(MolecularObjs::Atom const &a) const
  {
    auto const it = std::find_if(begin(_elements), end(_elements),
      [&a](auto const &el) { return a.element() == el; });

    return it == end(_elements) ? DType(0) : DType(1);
  }

private:
  std::vector<std::string> _elements;
};

#endif