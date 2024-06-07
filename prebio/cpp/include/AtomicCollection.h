#ifndef _ATOMIC_COLLECTION_H_
#define _ATOMIC_COLLECTION_H_

#include "Atom.h"
#include "PreBioExport.h"
#include "Residue.h"

#include "CoordinateTransformation.h"

#include <unordered_map>
#include <vector>
#include <exception>

/**
 * MAINLY A PLACEHOLDER CLASS FOR NOW
 */

namespace MolecularObjs
{

  class AlreadyInsertedResidue : public std::runtime_error
  {
    std::string generateMessage(ResidueInfo const& problemResidue) const
    {
      std::stringstream msg;
      msg << "Already inserted interior residue: " << problemResidue.chainID()
          << ", " << problemResidue.residueName() << " (" << problemResidue.sequenceNumber()
          << ").";
      return msg.str();
    }

  public:

    explicit AlreadyInsertedResidue(ResidueInfo const& problemResidue)
      : std::runtime_error(generateMessage(problemResidue))
    {}
  };

  class PREBIO_EXPORT AtomicCollection
  {
  public:
    using ResidueInfoMap
      = std::unordered_map<ResidueInfo, std::vector<unsigned int>>;

    AtomicCollection()
      : _nInteriorAtoms(0)
      , _coordSys(identityTransformation())
    {
    }

    class Builder
    {
    public:
      Builder();

      void addResidue(std::vector<Atom> &&atoms, ResidueInfo resInfo,
        unsigned int const nInterior);

      void setCoordinateSystem(CoordinateTransformation const &cs);

      AtomicCollection makeCollection();

    private:
      std::vector<Atom> _interiorAtoms;
      std::vector<Atom> _exteriorAtoms;
      ResidueInfoMap _residueToInteriorAtomMap;
      ResidueInfoMap _residueToExteriorAtomMap;
      CoordinateTransformation _coordSys;
    };

    ~AtomicCollection();
    AtomicCollection(AtomicCollection const &) = default;
    AtomicCollection(AtomicCollection &&) = default;
    AtomicCollection &operator=(AtomicCollection const &) = default;
    AtomicCollection &operator=(AtomicCollection &&) = default;

    unsigned int interiorCount() const noexcept { return _nInteriorAtoms; }

    std::vector<Atom> const &atoms() const noexcept { return _atoms; }

    std::vector<unsigned int> const& atomsAt(ResidueInfo const &res) const;

    CoordinateTransformation const &coordinateSystem() const noexcept;

    void removeResidue(ResidueInfo const &res);

  private:
    AtomicCollection(std::vector<Atom> &&atoms, unsigned int const nInterior,
      ResidueInfoMap &&residues, CoordinateTransformation const &cs);

    unsigned int _nInteriorAtoms;
    std::vector<Atom> _atoms;
    ResidueInfoMap _residues;
    CoordinateTransformation _coordSys;
  };


  inline CoordinateTransformation const &getCoordinateSystem(
    MolecularObjs::AtomicCollection const &ac)
  {
    return ac.coordinateSystem();
  }

} // namespace MolecularObjs


#endif