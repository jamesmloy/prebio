#ifndef _COORDINATE_TRANSFORMATION_H_
#define _COORDINATE_TRANSFORMATION_H_

#include <iostream>

#include "blaze/math/StaticMatrix.h"
#include "blaze/math/StaticVector.h"

template <typename AType, int N>
class CoordinateTransformation_t
{
public:
  using MatType = blaze::StaticMatrix<AType, N, N>;
  using PosType = blaze::StaticVector<AType, N>;

  CoordinateTransformation_t(MatType mat, PosType zero)
    : _mat(mat)
    , _zero(zero)
  {
  }

  PosType transform(PosType const &x) const { return _mat * x - _zero; }

  CoordinateTransformation_t inverse() const
  {
    try
    {
      MatType const theInv = blaze::inv(_mat);
      return {theInv, -theInv * _zero};
    }
    catch(const std::exception& e)
    {
      std::cout << "Could not invert matrix\n" << _mat << "\n";
      throw;
    }
  }

private:
  MatType _mat;
  PosType _zero;
};


using CoordinateTransformation = CoordinateTransformation_t<float, 3>;

template <typename AType = float>
inline
CoordinateTransformation_t<AType, 3>
identityTransformation()
{
  return CoordinateTransformation_t<AType, 3>({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, {0, 0, 0});
}


inline CoordinateTransformation
noShiftTransform(CoordinateTransformation::MatType const &mat)
{
  return CoordinateTransformation(mat, {0, 0, 0});
}


inline CoordinateTransformation
translateOnlyTransform(CoordinateTransformation::PosType const &pos)
{
  return CoordinateTransformation({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}, pos);
}

template <typename DType, int N>
CoordinateTransformation_t<DType, N>
fromCenterAndBasis(typename CoordinateTransformation_t<DType, N>::PosType const &pos,
                   typename CoordinateTransformation_t<DType, N>::MatType const &mat)
{
  return CoordinateTransformation_t<DType, N>(mat, mat * pos);
}


#endif