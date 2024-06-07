#include "TensorFlowConverter.h"

#include "PreBioExport.h"


namespace
{
  template <typename TargetDType>
  unsigned int getTotalDataSize(
    typename TensorFlowConverter<TargetDType>::FlattenerVec const &flatteners)
  {
    unsigned int eltCnt = 0;
    for (auto const &f : flatteners) eltCnt += f->elementSize();
    return eltCnt;
  }
} // namespace


template <typename TargetDType>
TensorFlowConverter<TargetDType>::TensorFlowConverter(bool channelFirst)
  : _channelFirst(channelFirst)
{
}


/// copy ctor
template <typename TargetDType>
TensorFlowConverter<TargetDType>::TensorFlowConverter(
  TensorFlowConverter const &o)
  : _channelFirst(o._channelFirst)
{
  _flatteners.reserve(o._flatteners.size());

  for (auto const &f : o._flatteners) _flatteners.push_back(f->clone());
}


/// copy assign
template <typename TargetDType>
TensorFlowConverter<TargetDType> &
TensorFlowConverter<TargetDType>::operator=(TensorFlowConverter const &o)
{
  TensorFlowConverter tmp(o);
  _channelFirst = o._channelFirst;
  _flatteners.swap(tmp._flatteners);
  return *this;
}


/// dtor
template <typename TargetDType>
TensorFlowConverter<TargetDType>::~TensorFlowConverter<TargetDType>()
{
}


template <typename TargetDType>
auto
TensorFlowConverter<TargetDType>::convert(DiscretizedSpace const &ds)
  -> ConvertResult
{
  unsigned int const nVoxels = ds.getDataStore().getElemCount();
  unsigned int const totalEltCount = getTotalDataSize<TargetDType>(_flatteners);
  unsigned int const totalDataCount = totalEltCount * nVoxels;

  std::vector<TargetDType> allData(totalDataCount, TargetDType(0));
  unsigned int chanTot = 0;

  if (_channelFirst)
  {
    unsigned int idxStart = 0;
    for (auto &f : _flatteners)
    {
      f->addChannelFirstData(allData, idxStart, ds);
      idxStart += nVoxels * f->elementSize();
      chanTot += f->elementSize();
    }
  }
  else
  {
    unsigned int idxStart = 0;
    chanTot = _flatteners.size();
    for (auto &f : _flatteners)
    {
      f->addChannelLastData(allData, idxStart, chanTot, ds);
      ++idxStart;
    }
  }

  return {std::move(allData), _channelFirst, chanTot, ds.nx(), ds.ny(), ds.nz()};
}


template class PREBIO_EXPORT TensorFlowConverter<double>;
template class PREBIO_EXPORT TensorFlowConverter<float>;