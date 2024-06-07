#ifndef _CENTERED_BOX_BUILDER_H_
#define _CENTERED_BOX_BUILDER_H_

#include <vector>

#include "BoxBuilder.h"

#include "blaze/math/dense/StaticMatrix.h"
#include "blaze/math/dense/StaticVector.h"

class PREBIO_EXPORT CenteredBoxBuilder : public BoxBuilder
{
public:
  using MatType = blaze::StaticMatrix<double, 3, 3>;
  using PosType = blaze::StaticVector<double, 3>;
  using CoordTrans = typename BoxBuilder::CoordTrans;

  CenteredBoxBuilder(MatType const &mat, unsigned int nx, unsigned int ny,
    unsigned int nz, float lx, float ly, float lz);

  CenteredBoxBuilder(unsigned int nx, unsigned int ny, unsigned int nz,
    float lx, float ly, float lz);

  std::vector<DiscretizedSpace> buildBoxes(
    CoordTrans const &locAndOrientation) override;

private:
  MatType const _mat;
  unsigned int _nx, _ny, _nz;
  float _lx, _ly, _lz;
};


#endif