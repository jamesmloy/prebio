#ifndef _ATOMIC_COLLECTION_TRANSFORMERS_H_
#define _ATOMIC_COLLECTION_TRANSFORMERS_H_

#include "PreBioExport.h"

#include "json.h"

#include "CoordinateTransformation.h"

#include "gemmi/model.hpp"

#include <string>
#include <vector>
#include <exception>

namespace MolecularObjs
{
  class AtomicCollection;
}

namespace gemmi
{
  struct Residue;
  namespace cif
  {
    struct Document;
  }
} // namespace gemmi

namespace ApBio
{
  class IncompleteResidue : public std::runtime_error
  {
  public:

    IncompleteResidue(std::string const& what_arg)
      : std::runtime_error(what_arg)
    {}

    IncompleteResidue(char const* what_arg)
      : std::runtime_error(what_arg)
    {}
  };

  struct CraSorter
  {
    inline bool operator()(gemmi::CRA const& aCra, gemmi::CRA const& bCra) const
    {
      if (aCra.chain->name == bCra.chain->name)
      {
        if (aCra.residue->seqid.num.value == bCra.residue->seqid.num.value)
        {
          return aCra.atom->serial < bCra.atom->serial;
        }
        else
          return aCra.residue->seqid.num.value < bCra.residue->seqid.num.value;
      }
      else
        return aCra.chain->name < bCra.chain->name;
    }
  };

  MolecularObjs::AtomicCollection PREBIO_EXPORT transformProteinChainResidue(gemmi::cif::Document &doc,
    Json const &snapshot, double filterRadius, std::vector<std::string> const &additionalParams,
    bool const includeSelf, bool const includeHetatm, bool const includeWater);

  MolecularObjs::AtomicCollection PREBIO_EXPORT transformWholeProtein(gemmi::cif::Document &doc,
    Json const &snapshot, std::vector<std::string> const &additionalParams,
    bool const includeHetatm, bool const includeWater);

  template <typename DType>
  blaze::StaticMatrix<DType, 3, 3> getBackboneOrientation(gemmi::Residue const &residue);

  template <typename DType>
  blaze::StaticVector<DType, 3> getCenterOfMass(gemmi::Residue const &residue);

  bool isAminoAcid(gemmi::Residue const &r);
}

#endif