#ifndef _RESIDUE_UTILS_H_
#define _RESIDUE_UTILS_H_

#include <tuple>

#include "Residue.h"

#include "CoordinateTransformation.h"

#include "blaze/math/dense/StaticVector.h"

namespace MolecularObjs
{
  namespace __detail
  {
    template <typename OutNumType, typename InNumType>
    blaze::StaticMatrix<OutNumType, 3, 3> toMatrix(
      blaze::StaticVector<InNumType, 3> const &a,
      blaze::StaticVector<InNumType, 3> const &b,
      blaze::StaticVector<InNumType, 3> const &c)
    {
      return blaze::StaticMatrix<OutNumType, 3, 3> {
        {a[0], a[1], a[2]}, {b[0], b[1], b[2]}, {c[0], c[1], c[2]}};
    }
  } // namespace __detail

  /// return val: (<C - CA>, <N - CA>) vectors
  template <typename NumType>
  std::tuple<blaze::StaticVector<NumType, 3>, blaze::StaticVector<NumType, 3>>
  extractPlaneVectors(Residue const &res)
  {
    using Vec_t = blaze::StaticVector<NumType, 3>;
    auto const &ac = res.getBackboneAtom("CA");
    auto const &cbb = res.getBackboneAtom("C");
    auto const &nbb = res.getBackboneAtom("N");

    return {Vec_t {cbb.position() - ac.position()},
      Vec_t {nbb.position() - ac.position()}};
  }


  template <typename NumType>
  CoordinateTransformation_t<NumType, 3> calculateTransformation(
    Residue const &res)
  {
    using namespace __detail;
    using Vec_t = blaze::StaticVector<NumType, 3>;

    auto&& [a, b] = MolecularObjs::extractPlaneVectors<NumType>(res);
    Vec_t const y = blaze::normalize(a);
    Vec_t z = blaze::normalize(blaze::cross(y, b));

    auto const &acPos = res.getBackboneAtom("CA").position();

    if (res.getInfo().residueName() != "GLY")
    {
      auto const &bc = res.getSideChainAtom("CB");
      if (blaze::dot(bc.position() - acPos, y) < NumType(0))
        z.scale(NumType(-1));
    }

    Vec_t const x = blaze::cross(y, z);

    return {toMatrix<NumType>(x, y, z), acPos};
  }
} // namespace MolecularObjs

#endif