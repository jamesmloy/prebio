#ifndef _GEMMI_DOCUMENT_CREATOR_H_
#define _GEMMI_DOCUMENT_CREATOR_H_

#include "PreBioExport.h"

#include "PdbCode.h"

#include "gemmi/cifdoc.hpp"

#include <string>

namespace ApBio
{
  class PREBIO_EXPORT GemmiDocumentCreator
  {
  public:
    GemmiDocumentCreator() {}
    virtual ~GemmiDocumentCreator() {}

    virtual gemmi::cif::Document makeDocument(PdbCode const &code) = 0;
    virtual gemmi::cif::Document makeDocument(std::string const &filename) = 0;
  };
}

#endif