#include "PhysicalPropertiesChannel.h"

#include "DiscretizedSpace.h"

#include "AtomicCollection.h"


PhysicalPropertiesChannel::PhysicalPropertiesChannel()
  : Base("PhysicalProperties")
{
}


std::vector<std::string>
PhysicalPropertiesChannel::getLabelTags()
{
  return {"PQR", "SASA"};
}


bool
PhysicalPropertiesChannel::hasData(MolecularObjs::Atom const &a) const
{
  return true;
}


PhysicalPropertiesChannel::DType
PhysicalPropertiesChannel::extract(MolecularObjs::Atom const &a) const
{
  return DType{a.pqr(), a.sasa()};
}
