#ifndef _TRUNC_SMOOTH_EXTRACTOR_H_
#define _TRUNC_SMOOTH_EXTRACTOR_H_

#include "SymmetricGaussSmoother.h"

#include "Gauss3dFunction.h"

#include "Atom.h"
#include "AtomicRadius.h"

#include "blaze/math/DynamicVector.h"
#include "blaze/math/StaticVector.h"

#include <chrono>
#include <cmath>


template <typename Scalar_t>
auto makeRadiusFunction(Scalar_t *sigma) -> std::function<Scalar_t(MolecularObjs::Atom const &)>
{
  if (sigma == nullptr)
  {
    return [](MolecularObjs::Atom const &atom) -> Scalar_t {
      return MolecularObjs::atomicRadius(atom);
    };
  }
  else
  {
    Scalar_t const sig_val = *sigma;
    return [sig_val](MolecularObjs::Atom const &atom) -> Scalar_t { return sig_val; };
  }
}

/// Primary template
template <typename... Args>
class TruncSmoothExtractor;


/// Specialization for generic cases. Channels can implement
/// their own specialization
template <typename ChannelTag_t, typename Evidence_t>
class TruncSmoothExtractor<ChannelTag_t, Evidence_t>
  : public ChannelExtractor<Evidence_t>
{
public:
  using ChanDType = typename ChannelTag_t::DType;
  using Scalar_t = typename NumericTraits<ChanDType>::Scalar_t;
  using RadFunc = typename std::function<Scalar_t(MolecularObjs::Atom const &)>;

  template <typename... ChannelArgs>
  TruncSmoothExtractor(Scalar_t trunc, Scalar_t *sigma, ChannelArgs &&...args)
    : _trunc(trunc)
    , _radiusFunction(makeRadiusFunction(sigma))
    , _channelTag(std::forward<ChannelArgs>(args)...)
  {
  }

  ~TruncSmoothExtractor()
  {
    // if (_totalTime.count() > 0.0)
    // {
    //   std::cout << "TRUNCATED SMOOTHING TIMING "
    //               "DATA:\n"
    //             << _channelTag.channelName() << "\n"
    //             << "============================================\n"
    //             << "Shell: " << _shellTime.count() << "("
    //             << _shellTime.count() / _totalTime.count() << ")"
    //             << "\n"
    //             << "Copy: " << _copyTime.count() << "("
    //             << _copyTime.count() / _totalTime.count() << ")"
    //             << "\n"
    //             << "Smooth: " << _smoothTime.count() << "("
    //             << _smoothTime.count() / _totalTime.count() << ")"
    //             << "\n"
    //             << "ReCopy: " << _recopyTime.count() << "("
    //             << _recopyTime.count() / _totalTime.count() << ")"
    //             << "\n"
    //             << "Total: " << _totalTime.count() << "\n"
    //             << "============================================\n"
    //             << std::endl;
    // }
  }

  void extractChannel(
    DiscretizedSpace &box, Evidence_t const &evidence) override
  {
    extractTruncSmoothChannel(box, _channelTag, _radiusFunction, evidence,
      _trunc, _totalTime, _shellTime, _copyTime, _smoothTime, _recopyTime);
  }

private:
  Scalar_t _trunc;
  RadFunc _radiusFunction;
  ChannelTag_t _channelTag;
  std::chrono::duration<double> _totalTime{0};
  std::chrono::duration<double> _shellTime{0};
  std::chrono::duration<double> _copyTime{0};
  std::chrono::duration<double> _smoothTime{0};
  std::chrono::duration<double> _recopyTime{0};
};


template <typename PointArray, typename IdxArray>
PointArray
copySubsetPoints(PointArray const &inPoints, IdxArray const &ids)
{
  PointArray outPoints(ids.size());
  for (unsigned int i = 0; i != ids.size(); ++i)
    outPoints[i] = inPoints[ids[i]];
  return outPoints;
}


template <typename Scalar_t>
unsigned int
calculateShellLayers(
  Scalar_t const &dx, Scalar_t const &trunc, MolecularObjs::Atom const &atom)
{
  auto const &rad = MolecularObjs::atomicRadius(atom);
  return (unsigned int)round(rad * trunc / dx + 0.5);
}


template <typename Channel_t, typename RadiusFunction, typename Evidence,
  typename Scalar_t>
void
extractTruncSmoothChannel(DiscretizedSpace &box, Channel_t &channel,
  RadiusFunction radiusFunction, Evidence const &evidence, Scalar_t const trunc,
  std::chrono::duration<double> &_totalTime,
  std::chrono::duration<double> &_shellTime,
  std::chrono::duration<double> &_copyTime,
  std::chrono::duration<double> &_smoothTime,
  std::chrono::duration<double> &_recopyTime)
{
  using ChanDType = typename Channel_t::DType;
  using SmootherFactory = SymmetricGaussSmoother<Scalar_t, Scalar_t>;
  auto const tStart = std::chrono::high_resolution_clock::now();

  ContiguousDataStore &ds = box.getDataStore();
  ds.addData<ChanDType>(channel.channelName());

  auto &channelData = ds.getData<ChanDType>(channel.channelName());
  auto const &coords = box.centroids();

  auto const &atoms = evidence.atoms();
  SmootherFactory factory;
  Scalar_t const dx = box.getSideLengths()[0] / Scalar_t(box.nx());
  unsigned int cnt = 0;

  for (auto const &a : atoms)
  {
    if (channel.hasData(a))
    {
      auto const tInnerStart = std::chrono::high_resolution_clock::now();
      unsigned int const shellSize = calculateShellLayers(dx, trunc, a);
      auto const cellIds = box.getShellLab(a.position(), shellSize);
      auto const tShell = std::chrono::high_resolution_clock::now();

      _shellTime += std::chrono::duration<double>(tShell - tInnerStart);

      if (cellIds.size() != 0)
      {
        ++cnt;
        auto smoother = factory.makeArraySmoother(a);

        auto const points = copySubsetPoints(coords, cellIds);
        auto const tCopy = std::chrono::high_resolution_clock::now();
        auto const vals = smoother(points);
        auto const tSmooth = std::chrono::high_resolution_clock::now();
        auto const coef = factory.normalizingCoef(a) * channel.extract(a);

        for (unsigned int i = 0; i != cellIds.size(); ++i)
          channelData[cellIds[i]] += vals[i] * coef;
        auto const tRecopy = std::chrono::high_resolution_clock::now();

        _copyTime += std::chrono::duration<double>(tCopy - tShell);
        _smoothTime += std::chrono::duration<double>(tSmooth - tCopy);
        _recopyTime += std::chrono::duration<double>(tRecopy - tSmooth);
      }
    }
  }

  auto const tEnd = std::chrono::high_resolution_clock::now();
  _totalTime += std::chrono::duration<double>(tEnd - tStart);
}

#endif
