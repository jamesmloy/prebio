#include "LocalDocumentCreator.h"

#include "gemmi/cif.hpp"

#include <fstream>
#include <iostream>

namespace
{
  std::string readFile(std::ifstream &stream)
  {
    std::string out;
    std::string line;
    while (std::getline(stream, line))
    {
      line.append("\n");
      out.append(line.data());
    }
    return out;
  }
} // namespace

namespace ApBio
{
  LocalDocumentCreator::LocalDocumentCreator(std::string folder,
    std::string prefix, std::string suffix, size_t cacheSize)
    : _folder(folder)
    , _prefix(prefix)
    , _suffix(suffix)
    , _cache(cacheSize)
  {
    if (folder.back() != '/')
      _folder = folder + '/';
  }

  std::string LocalDocumentCreator::readCifFile(PdbCode const &code) const
  {
    std::ifstream cifFile(getFilePath(code));
    if (!cifFile.is_open())
      throw std::runtime_error("Failed to open file: " + getFilePath(code));

    return readFile(cifFile);
  }

  std::string const &LocalDocumentCreator::getCifAsString(PdbCode const &code)
  {
    if (!_cache.contains(code))
      _cache.insert(code, readCifFile(code));

    return _cache.get(code);
  }

  gemmi::cif::Document LocalDocumentCreator::makeDocument(PdbCode const &code)
  {
    return gemmi::cif::read_string(getCifAsString(code));
  }

  gemmi::cif::Document LocalDocumentCreator::makeDocument(
    std::string const &filename)
  {
    return gemmi::cif::read_file(getFilePath(filename));
  }
} // namespace ApBio