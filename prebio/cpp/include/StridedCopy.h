#ifndef _STRIDED_COPY_H_
#define _STRIDED_COPY_H_

#include "NumericTraits.h"

#include "blaze/math/DynamicVector.h"
#include "blaze/math/StaticVector.h"

#include <type_traits>
#include <vector>
#include <iostream>

namespace __StridedCopyDetail
{
  template <typename TargetDType, typename T>
  void stridedCopy_scalar(typename std::vector<TargetDType>::iterator const beg,
    blaze::DynamicVector<T> const &input, std::true_type channelFirstTrue)
  {
    // TODO: this function can be optimized
    unsigned int const numElements = input.size();

    auto it = beg;
    for (unsigned int inIdx = 0; inIdx != numElements; ++inIdx)
      *(it++) = TargetDType(input[inIdx]);
  }
} // namespace __StridedCopyDetail


template <typename TargetDType, typename T, size_t N>
void
stridedCopy(typename std::vector<TargetDType>::iterator const beg,
  blaze::DynamicVector<blaze::StaticVector<T, N>> const &input,
  std::true_type channelFirstTrue)
{
  unsigned int const numElements = input.size();

  auto it = beg;
  for (unsigned int chan = 0; chan != N; ++chan)
  {
    for (unsigned int bIdx = 0; bIdx != numElements; ++bIdx)
      *(it++) = TargetDType(input[bIdx][chan]);
  }
}


template <typename TargetDType>
void
stridedCopy(typename std::vector<TargetDType>::iterator const beg,
  blaze::DynamicVector<double> const &input, std::true_type channelFirstTrue)
{
  using namespace __StridedCopyDetail;
  stridedCopy_scalar<TargetDType, double>(beg, input, channelFirstTrue);
}


template <typename TargetDType>
void
stridedCopy(typename std::vector<TargetDType>::iterator const beg,
  blaze::DynamicVector<float> const &input, std::true_type channelFirstTrue)
{
  using namespace __StridedCopyDetail;
  stridedCopy_scalar<TargetDType, float>(beg, input, channelFirstTrue);
}


template <typename TargetDType, typename T, size_t N>
void
channelLastCopy(typename std::vector<TargetDType> &data,
  blaze::DynamicVector<blaze::StaticVector<T, N>> const &input, unsigned int const chanIdx,
  unsigned int const numChan)
{
  throw std::runtime_error("Channel last only implemented for scalar channels");
}


template <typename TargetDType, typename SourceDType>
void
channelLastCopy(typename std::vector<TargetDType> &data,
  blaze::DynamicVector<SourceDType> const &input, unsigned int const chanIdx,
  unsigned int const numChan)
{
  for (unsigned int i = 0, e = chanIdx; e < data.size(); e += numChan, ++i)
    data[e] = TargetDType(input[i]);
}

#endif