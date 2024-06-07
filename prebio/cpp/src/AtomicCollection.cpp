#include "AtomicCollection.h"

#include <sstream>

namespace MolecularObjs
{
  AtomicCollection::Builder::Builder()
    : _coordSys(identityTransformation())
  {
  }

  void AtomicCollection::Builder::addResidue(std::vector<Atom> &&atoms,
    ResidueInfo resInfo, const unsigned int nInterior)
  {
    // add entry into interior map if there
    // are interior atoms
    if (nInterior > 0)
    {
      auto const insRet = _residueToInteriorAtomMap.insert({resInfo, {}});
      if (!insRet.second)
      {
        throw AlreadyInsertedResidue(resInfo);
      }

      auto &idxs = _residueToInteriorAtomMap[resInfo];
      idxs.reserve(nInterior);
      unsigned int currIdx = _interiorAtoms.size();
      for (unsigned int i = 0; i != nInterior; ++i)
      {
        _interiorAtoms.emplace_back(std::move(atoms[i]));
        idxs.push_back(currIdx++);
      }
    }

    // add entry if there are exterior atoms
    if (atoms.size() > nInterior)
    {
      auto const insRet = _residueToExteriorAtomMap.insert({resInfo, {}});
      if (!insRet.second)
      {
        throw AlreadyInsertedResidue(resInfo);;
      }

      auto &idxs = _residueToExteriorAtomMap[resInfo];
      idxs.reserve(atoms.size() - nInterior);
      unsigned int currIdx = _exteriorAtoms.size();
      for (unsigned int i = nInterior; i != atoms.size(); ++i)
      {
        _exteriorAtoms.emplace_back(std::move(atoms[i]));
        idxs.push_back(currIdx++);
      }
    }
  }


  void AtomicCollection::Builder::setCoordinateSystem(
    CoordinateTransformation const &cs)
  {
    _coordSys = cs;
  }


  AtomicCollection AtomicCollection::Builder::makeCollection()
  {
    unsigned int const nInterior = _interiorAtoms.size();
    unsigned int const nTotal = nInterior + _exteriorAtoms.size();

    // copy over all atoms to new vector
    std::vector<Atom> allAtoms;
    allAtoms.reserve(nTotal);

    std::for_each(begin(_interiorAtoms), end(_interiorAtoms),
      [&allAtoms](auto &atom) { allAtoms.emplace_back(std::move(atom)); });

    std::for_each(begin(_exteriorAtoms), end(_exteriorAtoms),
      [&allAtoms](auto &atom) { allAtoms.emplace_back(std::move(atom)); });

    // now combine the info maps into one
    ResidueInfoMap allResidues(std::move(_residueToInteriorAtomMap));

    for (auto &resPair : _residueToExteriorAtomMap)
    {
      // first, let's increment the indexing to be correct
      for (auto &idx : resPair.second) idx += nInterior;

      auto it = allResidues.find(resPair.first);

      // if there is no corresponding interior then we
      // just insert it
      if (it == allResidues.end())
      {
        allResidues.insert({resPair.first, std::move(resPair.second)});
      }
      else
      {
        it->second.insert(
          it->second.end(), resPair.second.begin(), resPair.second.end());
      }
    }

    return {std::move(allAtoms), nInterior, std::move(allResidues), _coordSys};
  }


  AtomicCollection::AtomicCollection(std::vector<Atom> &&atoms,
    unsigned int const nInterior, ResidueInfoMap &&residues,
    CoordinateTransformation const &cs)
    : _nInteriorAtoms(nInterior)
    , _atoms(atoms)
    , _residues(residues)
    , _coordSys(cs)
  {
  }

  AtomicCollection::~AtomicCollection() {}

  CoordinateTransformation const &
  AtomicCollection::coordinateSystem() const noexcept
  {
    return _coordSys;
  }

  std::vector<unsigned int> const &AtomicCollection::atomsAt(
    ResidueInfo const &res) const
  {
    auto it = _residues.find(res);
    if (it != _residues.end())
      return it->second;
    else
    {
      static std::vector<unsigned int> empty;
      return empty;
    }
  }


  void AtomicCollection::removeResidue(ResidueInfo const &residue)
  {
    auto const &indices = atomsAt(residue);
    if (!indices.empty())
    {
      unsigned int to_swap = _atoms.size() - 1;
      unsigned int const prevInteriorCount = _nInteriorAtoms;
      for (auto const &i : indices)
      {
        std::swap(_atoms[i], _atoms[to_swap--]);
        if (i < prevInteriorCount)
          --_nInteriorAtoms;
      }
      _atoms.resize(_atoms.size() - indices.size());
    }
  }
} // namespace MolecularObjs