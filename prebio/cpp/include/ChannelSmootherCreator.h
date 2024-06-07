#ifndef _CHANNEL_SMOOTHER_CREATOR_H_
#define _CHANNEL_SMOOTHER_CREATOR_H_

#include "Atom.h"

#include "blaze/math/StaticVector.h"
#include "blaze/math/DynamicVector.h"

#include <functional>

template <typename Out_t, typename Point_t>
class ChannelSmootherCreator
{
public:
  using Point = blaze::StaticVector<Point_t, 3>;
  using Point_f = blaze::StaticVector<float, 3>;
  using PointArray = blaze::DynamicVector<Point_f>;
  using OutArray = blaze::DynamicVector<Out_t>;
  using PointSmoother = std::function<Out_t(Point const &)>;
  using ArraySmoother = std::function<OutArray(PointArray const &)>;

  virtual ~ChannelSmootherCreator() {};

  virtual PointSmoother makeSmoother(MolecularObjs::Atom const& atom) = 0;
  virtual ArraySmoother makeArraySmoother(MolecularObjs::Atom const& atom) = 0;
  virtual Out_t normalizingCoef(MolecularObjs::Atom const& atom) = 0;
};

#endif