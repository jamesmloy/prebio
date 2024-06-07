#include "GemmiConvertUtils.h"

#include <algorithm>
#include <sstream>

namespace
{
  template <typename AType>
  AType getParamAs(std::string const &label,
    std::vector<std::string> const &keys, std::vector<std::string> const &vals,
    AType orElse = 0)
  {
    auto const it = std::find_if(
      begin(keys), end(keys), [&label](auto const &s) { return s == label; });

    AType param(orElse);
    if (it != end(keys))
    {
      auto const idx = std::distance(begin(keys), it);
      std::stringstream ss;
      ss << vals[idx];
      ss >> param;
    }
    return param;
  }

} // namespace

namespace ApBio
{
  MolecularObjs::Atom makeAtomFrom(gemmi::CRA const &cra,
    std::vector<std::string> const &keys, std::vector<std::string> const &vals)
  {
    return makeAtomFrom(*cra.atom, *cra.residue, cra.chain->name, keys, vals);
  }

  MolecularObjs::Atom makeAtomFrom(gemmi::Atom const &inAtom,
    gemmi::Residue const &residue, std::string const &chainId,
    std::vector<std::string> const &keys, std::vector<std::string> const &vals)
  {
    auto const isBb = isBackbone(inAtom);
    std::vector<std::string> keysCopy(keys);
    std::vector<std::string> valsCopy(vals);
    return MolecularObjs::Atom(inAtom.name, inAtom.element.name(), residue.name,
      inAtom.serial, residue.seqid.num.value, chainId[0], isBb, !isBb,
      getParamAs<double>("_atom_site.partial_charge", keys, vals),
      getParamAs<double>("_atom_site.solvent_accessibility", keys, vals),
      {inAtom.pos.x, inAtom.pos.y, inAtom.pos.z}, std::move(keysCopy),
      std::move(valsCopy));
  }

  MolecularObjs::ResidueInfo makeResidueFrom(gemmi::CRA const &cra)
  {
    auto const &inResidue = *cra.residue;
    unsigned int const seqid = inResidue.seqid.num.value;
    return MolecularObjs::ResidueInfo(cra.chain->name, inResidue.name, seqid);
  }

  MolecularObjs::ResidueInfo makeResidueFrom(gemmi::Residue const &inResidue, std::string const& chainId)
  {
    unsigned int const seqid = inResidue.seqid.num.value;
    return MolecularObjs::ResidueInfo(chainId, inResidue.name, seqid);
  }
} // namespace ApBio