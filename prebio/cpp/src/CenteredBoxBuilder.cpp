#include "CenteredBoxBuilder.h"

#include "ResidueUtils.h"

using namespace MolecularObjs;


CenteredBoxBuilder::CenteredBoxBuilder(MatType const &mat, unsigned int nx,
  unsigned int ny, unsigned int nz, float lx, float ly, float lz)
  : _mat(mat)
  , _nx(nx)
  , _ny(ny)
  , _nz(nz)
  , _lx(lx)
  , _ly(ly)
  , _lz(lz)
{
}

CenteredBoxBuilder::CenteredBoxBuilder(unsigned int nx, unsigned int ny,
  unsigned int nz, float lx, float ly, float lz)
  : CenteredBoxBuilder(
    MatType({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}), nx, ny, nz, lx, ly, lz)
{
}


std::vector<DiscretizedSpace>
CenteredBoxBuilder::buildBoxes(CoordTrans const &locAndOrientation)
{
  std::vector<DiscretizedSpace> boxes;
  boxes.emplace_back(_nx, _ny, _nz, _lx, _ly, _lz, locAndOrientation);

  return boxes;
}