#ifndef _SYMMETRIC_GAUSS_SMOOTHER_H_
#define _SYMMETRIC_GAUSS_SMOOTHER_H_

#include "ChannelSmootherCreator.h"

#include "Constants.h"
#include "Gauss3dFunction.h"

#include "AtomicRadius.h"


template <typename Out_t, typename Point_t>
class SymmetricGaussSmoother : public ChannelSmootherCreator<Out_t, Point_t>
{
public:
  using Parent = ChannelSmootherCreator<Out_t, Point_t>;
  using typename Parent::PointArray;
  using typename Parent::OutArray;
  using typename Parent::PointSmoother;
  using typename Parent::ArraySmoother;

  PointSmoother makeSmoother(MolecularObjs::Atom const &atom) override;
  ArraySmoother makeArraySmoother(MolecularObjs::Atom const& atom) override;
  Out_t normalizingCoef(MolecularObjs::Atom const& atom) override;
};


template <typename Out_t, typename Point_t>
auto
SymmetricGaussSmoother<Out_t, Point_t>::makeSmoother(
  MolecularObjs::Atom const &atom) -> PointSmoother
{
  return gauss3dFunction<Out_t, Point_t>(
    atom.position(), MolecularObjs::atomicRadius(atom) * 0.5);
}


namespace __SgsUtils
{
  template <typename OutArray, typename PointArray, typename Point,
    typename Out_t>
  inline OutArray initExp(
    PointArray const &a, Point const &mean, Out_t const &innerCoef)
  {
    OutArray val(a.size(), Out_t(0));
    for (int i = 0; i != (int)a.size(); ++i)
      val[i] = blaze::sqrNorm(a[i] - mean) * innerCoef;
    return val;
  }
} // namespace __SgsUtils


template <typename Out_t, typename Point_t>
Out_t
SymmetricGaussSmoother<Out_t, Point_t>::normalizingCoef(
  MolecularObjs::Atom const &atom)
{
  using namespace Constants;
  auto const stDev = MolecularObjs::atomicRadius(atom) * 0.5;
  Out_t const twoPi = Out_t(2) * Out_t(pi);
  Out_t const det = blaze::pow(stDev, Out_t(6));
  Out_t const piCoef = twoPi * twoPi * twoPi;
  Out_t const coef = blaze::pow(det * piCoef, Out_t(-0.5));

  return coef;
}


template <typename Out_t, typename Point_t>
auto
SymmetricGaussSmoother<Out_t, Point_t>::makeArraySmoother(
  MolecularObjs::Atom const &atom) -> ArraySmoother
{
  using namespace Constants;
  auto const stDev = MolecularObjs::atomicRadius(atom) * 0.5;
  auto const &mean = atom.position();
  return [stDev, mean](PointArray const &a) -> OutArray {

    Out_t const innerCoef = Out_t(-0.5) / (stDev * stDev);
    auto const toExp = __SgsUtils::initExp<OutArray>(a, mean, innerCoef);

    return blaze::exp(toExp);
  };
}

#endif