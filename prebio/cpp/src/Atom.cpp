#include "Atom.h"

namespace
{
  static const std::string emptyParam("");
}

namespace MolecularObjs
{
  Atom::Atom(std::string name, std::string elt, std::string rName, int sn,
    int rsn, char chainId, bool isBb, bool isSca, AType pqr, AType sasa, Coord_t pos,
    std::vector<std::string> &&keys, std::vector<std::string> &&values)
    : _atomName(std::move(name))
    , _element(std::move(elt))
    , _residueName(std::move(rName))
    , _atomSerialNumber(sn)
    , _residueSequenceNumber(rsn)
    , _chainIdentifier(chainId)
    , _isBackbone(isBb)
    , _isSideChainAtom(isSca)
    , _pqr(pqr)
    , _sasa(sasa)
    , _position(pos)
    , _additionalParams()
  {
    if (keys.size() != values.size())
      throw std::runtime_error("keys and values must be the same size");
    
    std::vector<std::string> localKeys(keys);
    std::vector<std::string> localValues(values);

    for (unsigned int i = 0; i != localKeys.size(); ++i)
      _additionalParams.insert({std::move(localKeys[i]), std::move(localValues[i])});
  }

 Atom::Atom(std::string name, std::string elt, std::string rName, int sn,
    int rsn, char chainId, bool isBb, bool isSca, AType pqr, AType sasa, Coord_t pos)
    : Atom(name, elt, rName, sn, rsn, chainId, isBb, isSca, pqr, sasa, pos,
    {}, {})
  {
  }

  Atom::Atom()
    : Atom("", "", "", -1, -1, 0, true, false, 0, 0, {0, 0, 0},
    {}, {})
  {
  }


  Atom::Atom(Atom const &other)
    : Atom(other.atomName(), other.element(), other.residueName(),
      other.atomSerialNumber(), other.residueSequenceNumber(),
      other.chainIdentifier(), other.isBackbone(), other.isSideChainAtom(),
      other.pqr(), other.sasa(), other.position())
  {
    _additionalParams = other._additionalParams;
  }


  Atom::Atom(Atom &&other)
    : Atom(std::move(other._atomName), std::move(other._element),
      std::move(other._residueName), other.atomSerialNumber(),
      other.residueSequenceNumber(), other.chainIdentifier(),
      other.isBackbone(), other.isSideChainAtom(), other.pqr(),
      other.sasa(), other.position())
  {
    _additionalParams = std::move(other._additionalParams);
  }


  Atom& Atom::operator=(Atom const& other)
  {
    _atomName = other.atomName();
    _element = other.element();
    _residueName = other.residueName();
    _atomSerialNumber = other.atomSerialNumber();
    _residueSequenceNumber = other.residueSequenceNumber();
    _chainIdentifier = other.chainIdentifier();
    _isBackbone = other.isBackbone();
    _isSideChainAtom = other.isSideChainAtom();
    _pqr = other.pqr();
    _sasa = other.sasa();
    _position = other.position();
    _additionalParams = other._additionalParams;

    return *this;
  }


  Atom& Atom::operator=(Atom &&other)
  {
    _atomName = std::move(other.atomName());
    _element = std::move(other.element());
    _residueName = std::move(other.residueName());
    _atomSerialNumber = other.atomSerialNumber();
    _residueSequenceNumber = other.residueSequenceNumber();
    _chainIdentifier = other.chainIdentifier();
    _isBackbone = other.isBackbone();
    _isSideChainAtom = other.isSideChainAtom();
    _pqr = other.pqr();
    _sasa = other.sasa();
    _position = other.position();
    _additionalParams = std::move(other._additionalParams);

    return *this;
  }


  std::string const&
  Atom::additionalParam(std::string const& key) const
  {
    auto const it = _additionalParams.find(key);

    if (it == end(_additionalParams))
      return emptyParam;
    else
      return std::any_cast<std::string const&>(it->second);
  }


  std::vector<std::string> Atom::validColumnNames() const
  {
    std::vector<std::string> names;
    for (auto const &it : _additionalParams)
    {
      names.push_back(it.first);
    }
    return names;
  }


  std::ostream&
  operator<<(std::ostream &s, Atom const& a)
  {
    s << "ATOM NAME: " << a.atomName() << "\n"
      << "ELEMENT: " << a.element() << "\n"
      << "RESIDUE NAME: " << a.residueName() << "\n"
      << "ATOM SERIAL NUMBER: " << a.atomSerialNumber() << "\n"
      << "RESIDUE SEQUENCE NUMBER: " << a.residueSequenceNumber() << "\n"
      << "CHAIN IDENTIFIER: " << a.chainIdentifier() << "\n"
      << std::boolalpha
      << "IS BACKBONE: " << a.isBackbone() << "\n"
      << "IS SIDE CHAIN ATOM: " << a.isSideChainAtom() << "\n"
      << "PQR: " << a.pqr() << "\n"
      << "SASA: " << a.sasa() << "\n"
      << "POSITION:\n" << a.position() << "\n";

    s << "ADDITIONAL:\n";
    for (auto const & it : a.validColumnNames())
    {
      s << "\t" << it << ", " << a.additionalParam(it);
    }
    return s;
  }

} // namespace MolecularObjs