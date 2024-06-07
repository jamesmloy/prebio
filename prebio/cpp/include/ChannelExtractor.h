#ifndef _CHANNEL_EXTRACTOR_H_
#define _CHANNEL_EXTRACTOR_H_

#include "DiscretizedSpace.h"


template <typename Evidence_t>
class ChannelExtractor
{
public:
  ChannelExtractor();
  virtual ~ChannelExtractor();

  virtual void extractChannel(DiscretizedSpace &box, Evidence_t const &evidence)
    = 0;
};


template <typename Evidence_t>
ChannelExtractor<Evidence_t>::ChannelExtractor()
{}


template <typename Evidence_t>
ChannelExtractor<Evidence_t>::~ChannelExtractor()
{}

#endif