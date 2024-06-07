#ifndef _GAUSS_3D_FUNCTION_H_
#define _GAUSS_3D_FUNCTION_H_

#include "Constants.h"

#include "blaze/math/StaticVector.h"

#include <functional>

/**
 * return a symmetric 3d gaussian function
 */

template <typename T, typename Point_t>
decltype(auto)
gauss3dFunction(blaze::StaticVector<Point_t, 3> const& mean, Point_t const& stDev)
{
  using PointVec = blaze::StaticVector<Point_t, 3>;
  return [=] (PointVec const &x) -> T
  {
    using namespace Constants;
    T const twoPi = T(2) * T(pi);
    T const det = blaze::pow(stDev, T(6));
    T const piCoef = twoPi * twoPi * twoPi;
    T const coef = blaze::pow(det*piCoef, T(-0.5));

    return coef * blaze::exp(T(-0.5) * blaze::sqrNorm(x - mean) / (stDev * stDev));
  };
}



#endif