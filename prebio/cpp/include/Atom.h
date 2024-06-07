#ifndef _ATOM_H_
#define _ATOM_H_

#include "PreBioExport.h"

#include <any>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

#include "blaze/math/StaticVector.h"

namespace MolecularObjs
{

  class PREBIO_EXPORT Atom
  {
  public:
    using AType = double;
    using Coord_t = blaze::StaticVector<AType, 3>;

    /// main ctor
    Atom(std::string name, std::string elt, std::string rName, int sn, int rsn,
      char chainId, bool isBb, bool isSca, AType pqr, AType sasa, Coord_t pos,
      std::vector<std::string> &&keys, std::vector<std::string> &&values);

    /// no additional properties
    Atom(std::string name, std::string elt, std::string rName, int sn, int rsn,
      char chainId, bool isBb, bool isSca, AType pqr, AType sasa, Coord_t pos);

    explicit Atom();

    /// copy ctor
    Atom(Atom const &other);

    /// move ctor
    Atom(Atom &&other);

    /// copy assign
    Atom &operator=(Atom const &other);

    /// move assign
    Atom &operator=(Atom &&other);

    /// defualt dtor
    ~Atom() = default;

    std::string const &atomName() const { return _atomName; }

    std::string const &element() const { return _element; }

    std::string const &residueName() const { return _residueName; }

    int atomSerialNumber() const { return _atomSerialNumber; }

    int residueSequenceNumber() const { return _residueSequenceNumber; }

    char chainIdentifier() const { return _chainIdentifier; }

    bool isBackbone() const { return _isBackbone; }

    bool isSideChainAtom() const { return _isSideChainAtom; }

    AType pqr() const { return _pqr; }

    AType sasa() const { return _sasa; }

    Coord_t const &position() const { return _position; }

    std::string const &additionalParam(std::string const &key) const;

    std::vector<std::string> validColumnNames() const;

  private:
    std::string _atomName;
    std::string _element;
    std::string _residueName;
    int _atomSerialNumber;
    int _residueSequenceNumber;
    char _chainIdentifier;
    bool _isBackbone;
    bool _isSideChainAtom;
    AType _pqr;
    AType _sasa;
    Coord_t _position;
    std::unordered_map<std::string, std::any> _additionalParams;
  };


  PREBIO_EXPORT
  std::ostream &operator<<(std::ostream &s, Atom const &a);

  template <typename DType>
  DType getProperty(Atom const &a, std::string const &key)
  {
    auto const &prop = a.additionalParam(key);
    if (prop == "")
    {
      std::stringstream msg;
      msg << "Column with name " << key
          << " does not contain any data. Valid columns are:\n";

      for (auto const n: a.validColumnNames())
        msg << n << ", ";
      msg << a;

      throw std::runtime_error(msg.str());
    }

    std::stringstream ss(prop);
    DType asType;
    ss >> asType;
    return asType;
  }

} // namespace MolecularObjs

#endif
