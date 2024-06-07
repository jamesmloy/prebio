#ifndef _PHYSICAL_PROPERTIES_CHANNEL_H_
#define _PHYSICAL_PROPERTIES_CHANNEL_H_

#include "ChannelTag.h"
#include "PreBioExport.h"

#include "Atom.h"

#include "blaze/math/StaticVector.h"

#include <unordered_map>
#include <string>
#include <sstream>


class PREBIO_EXPORT PhysicalPropertiesChannel
  : public TypedChannelTag<blaze::StaticVector<double, 2>>
{
public:
  using Base = TypedChannelTag<blaze::StaticVector<double, 2>>;
  using DType = typename Base::DType;

  explicit PhysicalPropertiesChannel();

  static std::vector<std::string> getLabelTags();

  bool hasData(MolecularObjs::Atom const &a) const;

  DType extract(MolecularObjs::Atom const &a) const;

private:
};


template <typename T>
class PREBIO_EXPORT PartialChargeChannel : public TypedChannelTag<T>
{
public:
  using Base = TypedChannelTag<T>;
  using DType = typename Base::DType;
  using ElementChargeMap = std::unordered_map<std::string, DType>;
  using FormalMap = std::unordered_map<std::string, ElementChargeMap>;

  explicit PartialChargeChannel()
    : Base("PartialChargeChannel")
    , _formalCharges()
  {
  }

  explicit PartialChargeChannel(FormalMap const& formalCharges, std::string const columnName)
    : Base("PartialChargeChannel")
    , _formalCharges(formalCharges)
    , _columnName(columnName)
  {
  }

  explicit PartialChargeChannel(std::string const columnName)
    : Base("PartialChargeChannel")
    , _formalCharges()
    , _columnName(columnName)
  {
  }

  static std::vector<std::string> getLabelTags() { return {"PQR"}; }

  bool hasData(MolecularObjs::Atom const &a) const { return true; }

  DType extract(MolecularObjs::Atom const &a) const;

private:

  DType getFromAtom(MolecularObjs::Atom const &a) const
  {
    using namespace MolecularObjs;
    return _columnName == "" ? a.pqr() : getProperty<DType>(a, _columnName);
  }

  FormalMap _formalCharges;
  std::string _columnName = "";
};


template <typename T>
class PREBIO_EXPORT SolventAccessibilityChannel : public TypedChannelTag<T>
{
public:
  using Base = TypedChannelTag<T>;
  using DType = typename Base::DType;

  explicit SolventAccessibilityChannel(std::string const columnName)
    : Base("SolventAccessibilityChannel")
    , _columnName(columnName)
  {
  }

  static std::vector<std::string> getLabelTags() { return {"SASA"}; }

  bool hasData(MolecularObjs::Atom const &a) const { return true; }

  DType extract(MolecularObjs::Atom const &a) const
  {
    using namespace MolecularObjs;
    return _columnName == "" ? a.sasa() : getProperty<DType>(a, _columnName);
  }

public:
  std::string _columnName = "";
};


template <typename T>
class PREBIO_EXPORT BetaFactorChannel : public TypedChannelTag<T>
{
public:
  using Base = TypedChannelTag<T>;
  using DType = typename Base::DType;

  explicit BetaFactorChannel()
    : Base("BetaFactorChannel")
  {
  }

  static std::vector<std::string> getLabelTags() { return {"BetaFactor"}; }

  bool hasData(MolecularObjs::Atom const &a) const { return true; }

  std::string additionalParam() const noexcept
  {
    return "_atom_site.B_iso_or_equiv";
  }

  DType extract(MolecularObjs::Atom const &a) const
  {
    auto const & bfString = a.additionalParam(additionalParam());
    return DType(std::stod(bfString));
  }
};


template <typename T>
auto
PartialChargeChannel<T>::extract(MolecularObjs::Atom const &a) const -> DType
{

  auto const componentIt = _formalCharges.find(a.residueName());
  if (componentIt != _formalCharges.end())
  {
    auto const elementIt = componentIt->second.find(a.element());
    if (elementIt != componentIt->second.end())
      return elementIt->second;
    else
      return getFromAtom(a);
  }
  else
    return getFromAtom(a);
}

#endif