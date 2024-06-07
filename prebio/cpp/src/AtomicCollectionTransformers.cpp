#include "AtomicCollectionTransformers.h"

#include "gemmi/model.hpp"

#include <exception>
#include <unordered_set>

namespace
{
  template <typename DType>
  inline blaze::StaticVector<DType, 3> atomPosition(gemmi::Atom const &atom)
  {
    return blaze::StaticVector<DType, 3>(
      {(DType)atom.pos.x, (DType)atom.pos.y, (DType)atom.pos.z});
  }
} // namespace

namespace ApBio
{
  bool isAminoAcid(gemmi::Residue const &r)
  {
    static std::unordered_set<std::string> const aas({
      "ALA",
      "ARG",
      "ASN",
      "ASP",
      "CYS",
      "GLN",
      "GLU",
      "GLY",
      "HIS",
      "ILE",
      "LEU",
      "LYS",
      "MET",
      "PHE",
      "PRO",
      "SER",
      "THR",
      "TRP",
      "TYR",
      "VAL",
    });

    return aas.count(r.name) != 0;
  }

  template <typename DType>
  blaze::StaticMatrix<DType, 3, 3> getBackboneOrientation(
    gemmi::Residue const &residue)
  {
    using Vec = blaze::StaticVector<DType, 3>;

    auto const *caPtr = residue.get_ca();
    auto const *cPtr = residue.get_c();
    auto const *nPtr = residue.get_n();

    if (!caPtr) throw IncompleteResidue("Could not find alpha carbon.");
    if (!cPtr) throw IncompleteResidue("Could not find carbon.");
    if (!nPtr) throw IncompleteResidue("Could not find nitrogen.");

    auto const caPos = atomPosition<DType>(*caPtr);
    auto const cPos = atomPosition<DType>(*cPtr);
    auto const nPos = atomPosition<DType>(*nPtr);

    auto const y = Vec(blaze::normalize(cPos - caPos));
    auto const caToN = Vec(nPos - caPos);

    auto z = Vec(blaze::normalize(blaze::cross(y, caToN)));

    if (residue.name != "GLY")
    {
      auto const *cbPtr = residue.find_atom("CB", '*', gemmi::El::C);

      if (!cbPtr) throw IncompleteResidue("Could not find beta carbon");

      auto const cbPos = atomPosition<DType>(*cbPtr);
      auto const caToCb = Vec(cbPos - caPos);
      if (blaze::dot(caToCb, y) < DType(0))
        z.scale(DType(-1));
    }

    auto const x = blaze::normalize(blaze::cross(y, z));

    return blaze::StaticMatrix<DType, 3, 3>({
      {x[0], x[1], x[2]},
      {y[0], y[1], y[2]},
      {z[0], z[1], z[2]},
    });
  }


  template <typename DType>
  blaze::StaticVector<DType, 3> getCenterOfMass(gemmi::Residue const &residue)
  {
    DType totalWeight(0);
    blaze::StaticVector<DType, 3> weightedCoords({0, 0, 0});

    for (auto const &atom : residue.atoms)
    {
      totalWeight += DType(atom.element.weight());
      weightedCoords += atomPosition<DType>(atom).scale(atom.element.weight());
    }

    return weightedCoords.scale(DType(1) / totalWeight);
  }

  /**
   * EXPLICIT INSTANTIATIONS
   */

  template blaze::StaticMatrix<float, 3, 3> PREBIO_EXPORT getBackboneOrientation<float>(
    gemmi::Residue const &residue);

  template blaze::StaticMatrix<double, 3, 3> PREBIO_EXPORT getBackboneOrientation<double>(
    gemmi::Residue const &residue);

  template blaze::StaticVector<float, 3> PREBIO_EXPORT getCenterOfMass<float>(
    gemmi::Residue const &residue);

  template blaze::StaticVector<double, 3> PREBIO_EXPORT getCenterOfMass<double>(
    gemmi::Residue const &residue);
} // namespace ApBio
