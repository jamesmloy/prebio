#ifndef _GEMMI_DOC_UTILS_H_
#define _GEMMI_DOC_UTILS_H_

#include "gemmi/cifdoc.hpp"

#include <sstream>
#include <algorithm>

namespace ApBio
{
  inline auto findLoop(gemmi::cif::Block &block, std::string const &tag)
    -> decltype(auto)
  {
    auto atomSiteTable = block.find_mmcif_category("_atom_site.");
    auto tags = atomSiteTable.tags();
    auto it = std::find(tags.begin(), tags.end(), tag);
    if (it == tags.end())
    {
      std::stringstream msg;
      msg << "Unable to find loop " << tag << ".  Tags available:\n";
      for (auto t: tags)
      {
        msg << t << "\n";
      }
      throw std::runtime_error(msg.str());
    }
    return atomSiteTable.column(std::distance(tags.begin(), it));
  }

  inline auto findLoop(gemmi::cif::Document &doc, std::string const &tag)
    -> decltype(auto)
  {
    return findLoop(doc.sole_block(), tag);
  }
}

#endif