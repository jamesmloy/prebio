#ifndef _ATOMIC_RADIUS_H_
#define _ATOMIC_RADIUS_H_

#include "Atom.h"
#include "PreBioExport.h"

namespace MolecularObjs
{
  PREBIO_EXPORT double atomicRadius(MolecularObjs::Atom const &atom);

  PREBIO_EXPORT std::unordered_map<std::string, double> const & radiiMap();
}

#endif