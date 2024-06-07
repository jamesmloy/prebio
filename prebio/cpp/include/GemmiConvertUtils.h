#ifndef _GEMMI_CONVERT_UTILS_H_
#define _GEMMI_CONVERT_UTILS_H_

#include "Atom.h"
#include "Residue.h"

#include "gemmi/model.hpp"

namespace ApBio
{
  MolecularObjs::Atom makeAtomFrom(gemmi::CRA const &cra,
    std::vector<std::string> const &keys, std::vector<std::string> const &vals);

  MolecularObjs::Atom makeAtomFrom(gemmi::Atom const &inAtom,
    gemmi::Residue const &residue, std::string const &chainId,
    std::vector<std::string> const &keys, std::vector<std::string> const &vals);

  MolecularObjs::ResidueInfo makeResidueFrom(
    gemmi::Residue const &inResidue, std::string const &chainId);
  MolecularObjs::ResidueInfo makeResidueFrom(gemmi::CRA const &cra);

  inline bool isBackbone(gemmi::Atom const &atom)
  {
    return atom.name == "C" || atom.name == "CA" || atom.name == "N"
      || atom.name == "O";
  }
} // namespace ApBio

#endif