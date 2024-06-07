#ifndef _LOCAL_DOCUMENT_CREATOR_H_
#define _LOCAL_DOCUMENT_CREATOR_H_

#include "GemmiDocumentCreator.h"

#include "LruCache.h"

#include <string>

namespace ApBio
{
  class PREBIO_EXPORT LocalDocumentCreator : public GemmiDocumentCreator
  {
  public:
    LocalDocumentCreator(std::string folder, std::string prefix = "",
      std::string suffix = "", size_t cacheSize = 20);

    gemmi::cif::Document makeDocument(PdbCode const &code) override;
    gemmi::cif::Document makeDocument(std::string const &filename) override;

  private:
    inline std::string getFilePath(PdbCode const &code) const
    {
      return _folder + _prefix + code.code() + _suffix + ".cif";
    }

    inline std::string getFilePath(std::string const &filename) const
    {
      return _folder + filename;
    }

    std::string const &getCifAsString(PdbCode const &code);
    std::string readCifFile(PdbCode const &code) const;

    std::string _folder;
    std::string _prefix;
    std::string _suffix;
    LruCache<PdbCode, std::string> _cache{10};
  };
} // namespace ApBio

#endif