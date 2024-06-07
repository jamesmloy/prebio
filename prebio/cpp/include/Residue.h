#ifndef _RESIDUE_H_
#define _RESIDUE_H_

#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include "Atom.h"


namespace MolecularObjs
{
  /**
   * ResidueInfo is a class that is comprised of the residue's name and sequence position.
   * It creates a unique hash from the combination of the name and sequence position. 
   * With these hashes, you can directly compare 2 objects of ResidueInfo and know if they
   * are the same residue are not. 
   * 
   * This class needs to be updated so the hash incorporates the chain ID and maybe 
   * even the PDB code that the residue came from. This will enable residue comparison
   * between chains and proteins. 
   */
  class ResidueInfo // This class should also include the chain. The information is incomplete.
  {
    friend std::hash<ResidueInfo>;

    inline std::size_t getHash() const;

  public:
    inline ResidueInfo(std::string chain_id, std::string residueName, int sequenceNumber);

    inline std::string const &chainID() const;
    inline std::string const &residueName() const;
    inline int const &sequenceNumber() const;

    inline bool operator==(ResidueInfo const &o) const;
    inline bool operator!=(ResidueInfo const &o) const;

  private:
    std::string _chainID;
    std::string _residueName;
    int _sequenceNumber;
    std::size_t _hash;
  };


  class Residue
  {
  public:
    using AtomVec = std::vector<Atom const *>;

    Residue(ResidueInfo &&info, AtomVec &&atoms)
      : _info(info)
      , _atoms(atoms)
    {
    }

    Residue(std::string chainID, std::string residueName, int sequenceNumber, AtomVec &&atoms)
      : Residue(
        {std::move(chainID), std::move(residueName), sequenceNumber}
      , std::forward<AtomVec>(atoms))
    {
    }

    ResidueInfo const &getInfo() const { return _info; }
    Atom const &getBackboneAtom(std::string const &aTag) const
    {
      auto ac
        = std::find_if(begin(_atoms), end(_atoms), [&aTag](auto const &a) {
            return (a->atomName() == aTag) && a->isBackbone();
          });

      if (ac == end(_atoms))
        throw std::runtime_error("could not find backbone atom " + aTag
          + " for residue:" + _info.residueName());

      return **ac;
    }

    Atom const &getSideChainAtom(std::string const &aTag) const
    {
      auto ac
        = std::find_if(begin(_atoms), end(_atoms), [&aTag](auto const &a) {
            return (a->atomName() == aTag) && !a->isBackbone();
          });

      if (ac == end(_atoms))
        throw std::runtime_error("could not find side chain atom " + aTag
          + " for residue:" + _info.residueName());

      return **ac;
    }


  private:
    ResidueInfo _info;
    std::vector<Atom const *> _atoms;
  };

  
// In the hash we should include the chain id. Bc lots of proteins have duplicate chains. 
// (homodimers, homotetramers, etc.)
  ResidueInfo::ResidueInfo(std::string chainID, std::string residueName, int sequenceNumber)
    :  _chainID(chainID) 
    ,  _residueName(std::move(residueName))
    , _sequenceNumber(sequenceNumber)
    , _hash(std::hash<std::string>()(
        _chainID + _residueName + std::to_string(_sequenceNumber)))
  {
  }

  std::size_t ResidueInfo::getHash() const { return _hash; }

  std::string const &ResidueInfo::chainID() const { return _chainID; }

  std::string const &ResidueInfo::residueName() const { return _residueName; }

  int const &ResidueInfo::sequenceNumber() const { return _sequenceNumber; }

  bool ResidueInfo::operator==(ResidueInfo const &o) const
  {
    return getHash() == o.getHash();
  }

  bool ResidueInfo::operator!=(ResidueInfo const &o) const
  {
    return !(*this == o);
  }

} // namespace MolecularObjs


namespace std
{
  template <>
  struct hash<MolecularObjs::ResidueInfo>
  {
    std::size_t operator()(MolecularObjs::ResidueInfo const &info) const
      noexcept
    {
      return info.getHash();
    }
  };


  template <>
  struct hash<MolecularObjs::Residue>
  {
    std::size_t operator()(MolecularObjs::Residue const &res) const noexcept
    {
      return hash<MolecularObjs::ResidueInfo>()(res.getInfo());
    }
  };
} // namespace std


#endif